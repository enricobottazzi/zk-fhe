use crate::poly::Poly;
use axiom_eth::{rlp::rlc::RlcChip, Field};
use halo2_base::{
    gates::GateChip,
    safe_types::{GateInstructions, RangeChip, RangeInstructions},
    utils::bigint_to_fe,
    AssignedValue, Context,
    QuantumCell::{Constant, Existing},
};
use num_bigint::BigInt;
use num_traits::Num;

/// Chip for polynomial operations in the circuit
#[derive(Clone, Debug)]
pub struct PolyChip<const DEG: usize, F: Field>
where
    [(); DEG + 1]: Sized,
{
    pub assigned_coefficients: [AssignedValue<F>; DEG + 1],
    pub max_num_bits: usize,
}

impl<const DEG: usize, F: Field> PolyChip<DEG, F>
where
    [(); DEG + 1]: Sized,
    [(); DEG * 2 + 1]: Sized,
{
    /// Assign the coefficients of the polynomial
    pub fn new(poly: Poly, ctx: &mut Context<F>) -> Self {
        // assert that the degree of the input polynomial is equal to DEG
        assert_eq!(poly.deg(), DEG);

        let mut assigned_coefficients = vec![];
        let mut max_num_bits = 0;

        for item in poly.coefficients.iter().take(DEG + 1) {
            let num_bits = item.bits();
            let val = bigint_to_fe(item);
            let assigned_coeff = ctx.load_witness(val);
            assigned_coefficients.push(assigned_coeff);
            if num_bits > max_num_bits {
                max_num_bits = num_bits;
            }
        }

        Self {
            assigned_coefficients: assigned_coefficients.try_into().unwrap(),
            max_num_bits: max_num_bits as usize,
        }
    }

    /// Expose the polynomial coefficients to public
    pub fn to_public(&self, make_public: &mut Vec<AssignedValue<F>>) {
        for &coeff in self.assigned_coefficients.iter() {
            make_public.push(coeff);
        }
    }

    /// Enforce that self * b = c
    /// This constraint leverages the Axiom's Challenge API
    /// The challenge API allows us to access a random challenge gamma inside the circuit.
    /// Any polynomial identity can be checked using gamma. e.g. p(x) = q(x) with high probability if p(gamma) = q(gamma) for a challenge gamma
    /// because a random gamma has vanishing probability of being a root of p(x) - q(x).
    /// Analogously, we can check the identity a(gamma) * b(gamma) = c(gamma) for a random gamma.
    ///
    /// Complexity:
    /// Computing the polynomial multiplication using the direct method would take O(N^2) constraints
    /// This algorithm takes O(N) constraints as it requires to:
    /// - Evaluate the polynomials a, b and c at gamma (3N constraints)
    /// - Enforce the identity a(gamma) * b(gamma) - c(gamma) = 0 (1 constraint)
    ///
    /// Note that the constraint will fail if the coefficients of the resulting polynomial c overflows the prime field p.
    ///
    /// # Assumptions
    /// * The coefficients of the polynomial `c` have not overflowed the prime field p during the assignment phase
    pub fn constrain_poly_mul(
        &self,
        b: PolyChip<DEG, F>,
        c: PolyChip<{ DEG * 2 }, F>,
        ctx_gate: &mut Context<F>,
        ctx_rlc: &mut Context<F>,
        rlc: &RlcChip<F>,
    ) {
        // Assert that `max_num_bits` of the coefficients of the polynomial `c` is less than the number of bits of the modulus of the circuit field
        let p = BigInt::from_str_radix(&F::MODULUS[2..], 16).unwrap();
        let p_bits = p.bits();
        let c_max_bits = c.max_num_bits;

        assert!(c_max_bits < p_bits.try_into().unwrap());

        // `compute_rlc` evaluates the polynomial at gamma and returns the evaluation
        let poly_a_trace = rlc.compute_rlc_fixed_len(ctx_rlc, self.assigned_coefficients);
        let poly_a_eval_assigned = poly_a_trace.rlc_val;

        let poly_b_trace = rlc.compute_rlc_fixed_len(ctx_rlc, b.assigned_coefficients);
        let poly_b_eval_assigned = poly_b_trace.rlc_val;

        let poly_c_trace = rlc.compute_rlc_fixed_len(ctx_rlc, c.assigned_coefficients);
        let poly_c_eval_assigned = poly_c_trace.rlc_val;

        // enforce gate a(gamma) * b(gamma) - c(gamma) = 0
        ctx_gate.assign_region(
            [
                Constant(F::from(0)),
                Existing(poly_a_eval_assigned),
                Existing(poly_b_eval_assigned),
                Existing(poly_c_eval_assigned),
            ],
            [0],
        );
    }

    /// Takes a polynomial represented by its coefficients in a vector and output a new polynomial reduced by applying modulo Q to each coefficient
    pub fn reduce_by_modulo_q<const Q: u64>(
        &self,
        ctx: &mut Context<F>,
        range: &RangeChip<F>,
    ) -> PolyChip<DEG, F>
    where
        [(); DEG + 1]: Sized,
    {
        let mut output = vec![];
        // Enforce that in_assigned[i] % Q = rem_assigned[i]
        // Note that `div_mod` requires the value to be reduced to be at most `num_bits`
        let num_bits = self.max_num_bits;
        for i in 0..=DEG {
            let rem = range
                .div_mod(ctx, self.assigned_coefficients[i], Q, num_bits)
                .1;
            output.push(rem);
        }

        // `max_num_bits` of the output polynomial is equal to the number of bits of `Q` after the reduction
        let max_num_bits = BigInt::from(Q).bits() as usize;

        PolyChip {
            assigned_coefficients: output.try_into().unwrap(),
            max_num_bits,
        }
    }

    /// Enforce that polynomial has coefficients in the range [0, Z] or [Q-Z, Q-1]
    ///
    /// # Generic parameters
    /// * Q is the modulus of the ring R_q (ciphertext space)
    /// * Z is the constant that defines the range
    ///
    /// # Assumptions
    /// * Z < Q
    pub fn check_poly_coefficients_in_range<const Q: u64, const Z: u64>(
        &self,
        ctx: &mut Context<F>,
        range: &RangeChip<F>,
    ) {
        // assert that Z < Q
        assert!(Z < Q);

        // The goal is to check that coeff is in the range [0, Z] OR [Q-Z, Q-1]
        // We split this check into two checks:
        // - Check that coeff is in the range [0, Z] and store the boolean result in in_partial_range_1_vec
        // - Check that coeff is in the range [Q-Z, Q-1] and store the boolean result in in_partial_range_2_vec
        // We then perform (`in_partial_range_1` OR `in_partial_range_2`) to check that coeff is in the range [0, Z] OR [Q-Z, Q-1]
        // The result of this check is stored in the `in_range` vector.
        // All the boolean values in `in_range` are then enforced to be true
        let mut in_range_vec = Vec::with_capacity(DEG + 1);

        // get the number of bits needed to represent the value of Q
        let binary_representation = format!("{:b}", Q);
        let q_bits = binary_representation.len();

        for coeff in self.assigned_coefficients.iter() {
            // First of all, enforce that coefficient is in the [0, 2^q_bits] range
            let bool = range.is_less_than_safe(ctx, *coeff, (1 << q_bits as u64) + 1);
            range.gate().assert_is_const(ctx, &bool, &F::from(1));

            // Check for the range [0, Z]
            // coeff is known are known to have <= `q_bits` bits according to the constraint set above
            // Z + 1 is known to have <= `q_bits` bits according to assumption of the chip
            // Therefore it satisfies the assumption of `is_less_than` chip
            let in_partial_range_1 =
                range.is_less_than(ctx, *coeff, Constant(F::from(Z + 1)), q_bits);

            // Check for the range [Q-Z, Q-1]
            // coeff is known are known to have <= `q_bits` bits according to the constraint set above
            // Q - Z is known to have <= `q_bits` bits according to assumption of the chip
            // Therefore it satisfies the assumption of `is_less_than` chip
            let not_in_range_lower_bound =
                range.is_less_than(ctx, *coeff, Constant(F::from(Q - Z)), q_bits);
            let in_range_lower_bound = range.gate.not(ctx, not_in_range_lower_bound);

            // coeff is known are known to have <= `q_bits` bits according to the constraint set above
            // Q is known to have <= `q_bits` by definition
            // Therefore it satisfies the assumption of `is_less_than` chip
            let in_range_upper_bound =
                range.is_less_than(ctx, *coeff, Constant(F::from(Q)), q_bits);
            let in_partial_range_2 =
                range
                    .gate
                    .and(ctx, in_range_lower_bound, in_range_upper_bound);

            // Combined check for [0, Z] OR [Q-Z, Q-1]
            let in_range = range.gate.or(ctx, in_partial_range_1, in_partial_range_2);
            in_range_vec.push(in_range);
        }

        // Enforce that in_range_vec[i] = true
        for in_range in in_range_vec {
            let bool = range.gate.is_equal(ctx, in_range, Constant(F::from(1)));
            range.gate.assert_is_const(ctx, &bool, &F::from(1));
        }
    }

    /// Enforce that polynomial a of degree DEG is sampled from the distribution chi key. Namely, that the coefficients are in the range [0, 1, Q-1].
    ///
    /// # Generic parameters
    /// * Q is the modulus of the ring R_q (ciphertext space)
    pub fn check_poly_from_distribution_chi_key<const Q: u64>(
        &self,
        ctx: &mut Context<F>,
        gate: &GateChip<F>,
    ) {
        // In order to check that coeff is equal to either 0, 1 or q-1
        // The constraint that we want to enforce is:
        // (coeff - 0) * (coeff - 1) * (coeff - (q-1)) = 0

        // loop over all the coefficients of the polynomial
        for coeff in self.assigned_coefficients.iter() {
            // constrain (a - 0)
            let factor_1 = gate.sub(ctx, *coeff, Constant(F::from(0)));

            // constrain (a - 1)
            let factor_2 = gate.sub(ctx, *coeff, Constant(F::from(1)));

            // constrain (a - (q-1))
            let factor_3 = gate.sub(ctx, *coeff, Constant(F::from(Q - 1)));

            // constrain (a - 0) * (a - 1)
            let factor_1_2 = gate.mul(ctx, factor_1, factor_2);

            // constrain (a - 0) * (a - 1) * (a - (q-1))
            let factor_1_2_3 = gate.mul(ctx, factor_1_2, factor_3);

            // constrain (a - 0) * (a - 1) * (a - (q-1)) = 0
            let bool = gate.is_zero(ctx, factor_1_2_3);
            gate.assert_is_const(ctx, &bool, &F::from(1));
        }
    }
}
