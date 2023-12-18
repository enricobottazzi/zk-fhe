use std::cmp::max;

use crate::poly::Poly;
use axiom_eth::{rlp::rlc::RlcChip, Field};
use halo2_base::{
    gates::GateChip,
    safe_types::{GateInstructions, RangeChip, RangeInstructions},
    utils::{bigint_to_fe, fe_to_bigint},
    AssignedValue, Context,
    QuantumCell::{Constant, Existing},
};
use num_bigint::BigInt;
use num_traits::Num;

/// Chip for polynomial operations in the circuit
/// The polynomial is represented as a vector of AssignedValue
/// `max_num_bits` is the maximum number of bits of the coefficients of the polynomial. We need to keep track of this value to warn for risk of  overflow during the polynomial operations
#[derive(Clone, Debug)]
pub struct PolyChip<F: Field> {
    pub assigned_coefficients: Vec<AssignedValue<F>>,
    pub max_num_bits: usize,
    pub degree: usize,
}

impl<F: Field> PolyChip<F> {
    /// Build `PolyChip` from a `Poly`
    pub fn from_poly(poly: Poly, ctx: &mut Context<F>) -> Self {
        let mut assigned_coefficients = vec![];
        let mut max_num_bits = 0;
        let deg = poly.deg();

        for item in poly.coefficients.iter().take(deg + 1) {
            let num_bits = item.bits();
            let val = bigint_to_fe(item);
            let assigned_coeff = ctx.load_witness(val);
            assigned_coefficients.push(assigned_coeff);
            if num_bits > max_num_bits {
                max_num_bits = num_bits;
            }
        }

        Self {
            assigned_coefficients,
            max_num_bits: max_num_bits as usize,
            degree: deg,
        }
    }

    /// Build `PolyChip` from a vector of `AssignedValue`
    pub fn from_assigned_values(
        assigned_coefficients: Vec<AssignedValue<F>>,
        max_num_bits: usize,
    ) -> Self {
        let deg = assigned_coefficients.len() - 1;
        Self {
            assigned_coefficients,
            max_num_bits,
            degree: deg,
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
    /// * The coefficients of the polynomial `c` have not overflowed the prime field p during the assignment phase. Otherwise, the constraint will fail.
    pub fn constrain_poly_mul(
        &self,
        b: PolyChip<F>,
        c: PolyChip<F>,
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
        let poly_a_trace = rlc.compute_rlc_fixed_len(ctx_rlc, self.assigned_coefficients.clone());
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

    /// Sum polynomial and other
    ///
    /// # Assumptions
    /// * the coefficients are constrained such to avoid overflow during the polynomial addition
    pub fn add(&self, ctx: &mut Context<F>, other: PolyChip<F>, gate: &GateChip<F>) -> PolyChip<F> {
        let mut output = vec![];

        for i in 0..=self.degree {
            let val = gate.add(
                ctx,
                self.assigned_coefficients[i],
                other.assigned_coefficients[i],
            );
            output.push(val);
        }

        let max_num_bits_output = max(self.max_num_bits, other.max_num_bits) + 1;

        PolyChip::from_assigned_values(output, max_num_bits_output)
    }

    /// Multiply polynomial by a scalar
    ///
    /// # Assumptions
    /// * the coefficients are constrained such to avoid overflow during the polynomial scalar multiplication
    pub fn scalar_mul(
        &self,
        ctx: &mut Context<F>,
        scalar: &AssignedValue<F>,
        gate: &GateChip<F>,
    ) -> PolyChip<F> {
        let mut output = vec![];
        let scalar_big_int_bits = fe_to_bigint(scalar.value()).bits();

        for item in self.assigned_coefficients.iter().take(self.degree + 1) {
            let val = gate.mul(ctx, *item, *scalar);
            output.push(val);
        }

        let max_num_bits_output = self.max_num_bits + scalar_big_int_bits as usize;

        // Assert that `max_num_bits` of the coefficients of the polynomial is less than the number of bits of the modulus of the circuit field
        let p = BigInt::from_str_radix(&F::MODULUS[2..], 16).unwrap();
        let p_bits = p.bits();
        assert!(
            max_num_bits_output < p_bits as usize,
            "Risk of overflow detected in scalar_mul"
        );

        PolyChip::from_assigned_values(output, max_num_bits_output)
    }

    /// Reduce the coefficients of the polynomial by `modulus``
    pub fn reduce_by_modulo(
        &self,
        ctx: &mut Context<F>,
        range: &RangeChip<F>,
        modulus: u64,
    ) -> PolyChip<F> {
        let mut output = vec![];
        // Enforce that self.assigned_coefficients[i] % modulus = output[i]
        // Note that `div_mod` requires the value to be reduced to be at most `num_bits`
        let num_bits = self.max_num_bits;
        for i in 0..=self.degree {
            let reduced_coeff = range
                .div_mod(ctx, self.assigned_coefficients[i], modulus, num_bits)
                .1;
            output.push(reduced_coeff);
        }

        // `max_num_bits` of the output polynomial is equal to the number of bits of `modulus` after the reduction
        let max_num_bits = BigInt::from(modulus).bits() as usize;

        PolyChip::from_assigned_values(output, max_num_bits)
    }

    /// Constrain polynomial to be equal to other
    pub fn constrain_poly_equal(
        &self,
        ctx: &mut Context<F>,
        other: PolyChip<F>,
        gate: &GateChip<F>,
    ) {
        for i in 0..=self.degree {
            let bool = gate.is_equal(
                ctx,
                self.assigned_coefficients[i],
                other.assigned_coefficients[i],
            );
            gate.assert_is_const(ctx, &bool, &F::from(1))
        }
    }

    /// Enforce that polynomial has coefficients in the range [0, z] or [y-z, y-1]
    ///
    /// # Assumptions
    /// * z < y
    pub fn check_poly_coefficients_in_range(
        &self,
        ctx: &mut Context<F>,
        range: &RangeChip<F>,
        z: u64,
        y: u64,
    ) {
        // assert that z < y
        assert!(z < y);

        // The goal is to check that coeff is in the range [0, z] OR [y-z, y-1]
        // We split this check into two checks:
        // - Check that coeff is in the range [0, z] and store the boolean result in in_partial_range_1_vec
        // - Check that coeff is in the range [y-z, y-1] and store the boolean result in in_partial_range_2_vec
        // We then perform (`in_partial_range_1` OR `in_partial_range_2`) to check that coeff is in the range [0, z] OR [y-z, y-1]
        // Eventually, enforce that `in_range` = true

        // get the number of bits needed to represent the value of y
        let y_bits = BigInt::from(y).bits();

        let z_plus_one_const = Constant(F::from(z + 1));
        let y_minus_z_const = Constant(F::from(y - z));
        let y_const = Constant(F::from(y));
        let one_const = Constant(F::from(1));

        for coeff in self.assigned_coefficients.iter() {
            // First of all, enforce that coefficient is in the [0, 2^y_bits] range
            let bool = range.is_less_than_safe(ctx, *coeff, (1 << y_bits as u64) + 1);
            range.gate().assert_is_const(ctx, &bool, &F::from(1));

            // Check for the range [0, z]
            // coeff is known are known to have <= `y_bits` bits according to the constraint set above
            // z + 1 is known to have <= `y_bits` bits according to assumption of the function
            // Therefore it satisfies the assumption of `is_less_than` chip
            let in_partial_range_1 =
                range.is_less_than(ctx, *coeff, z_plus_one_const, y_bits as usize);

            // Check for the range [y-z, y-1]
            // coeff is known are known to have <= `y_bits` bits according to the constraint set above
            // y - z is known to have <= `y_bits` bits according to assumption of the function
            // Therefore it satisfies the assumption of `is_less_than` chip
            let not_in_range_lower_bound =
                range.is_less_than(ctx, *coeff, y_minus_z_const, y_bits as usize);
            let in_range_lower_bound = range.gate.not(ctx, not_in_range_lower_bound);

            // coeff is known are known to have <= `y_bits` bits according to the constraint set above
            // y is known to have <= `y_bits` by definition
            // Therefore it satisfies the assumption of `is_less_than` chip
            let in_range_upper_bound = range.is_less_than(ctx, *coeff, y_const, y_bits as usize);
            let in_partial_range_2 =
                range
                    .gate
                    .and(ctx, in_range_lower_bound, in_range_upper_bound);

            // Combined check for [0, z] OR [y-z, y-1]
            let in_range = range.gate.or(ctx, in_partial_range_1, in_partial_range_2);

            // Enforce that in_range = true
            let bool = range.gate.is_equal(ctx, in_range, one_const);
            range.gate.assert_is_const(ctx, &bool, &F::from(1));
        }
    }

    /// Enforce that polynomial is sampled from the distribution chi key. Namely, that the coefficients are in the range [0, 1, Z].
    pub fn check_poly_from_distribution_chi_key(
        &self,
        ctx: &mut Context<F>,
        gate: &GateChip<F>,
        z: u64,
    ) {
        // In order to check that coeff is equal to either 0, 1 or z
        // The constraint that we want to enforce is:
        // (coeff - 0) * (coeff - 1) * (coeff - (z)) = 0

        let zero_const = Constant(F::from(0));
        let one_const = Constant(F::from(1));
        let z_const = Constant(F::from(z));

        // loop over all the coefficients of the polynomial
        for coeff in self.assigned_coefficients.iter() {
            // constrain (a - 0)
            let factor_1 = gate.sub(ctx, *coeff, zero_const);

            // constrain (a - 1)
            let factor_2 = gate.sub(ctx, *coeff, one_const);

            // constrain (a - z)
            let factor_3 = gate.sub(ctx, *coeff, z_const);

            // constrain (a - 0) * (a - 1)
            let factor_1_2 = gate.mul(ctx, factor_1, factor_2);

            // constrain (a - 0) * (a - 1) * (a - z)
            let factor_1_2_3 = gate.mul(ctx, factor_1_2, factor_3);

            // constrain (a - 0) * (a - 1) * (a - z) = 0
            let bool = gate.is_zero(ctx, factor_1_2_3);
            gate.assert_is_const(ctx, &bool, &F::from(1));
        }
    }

    /// Enforce that the coefficients of the polynomial are in the modulus field
    pub fn check_coefficients_in_modulus_field(
        &self,
        ctx: &mut Context<F>,
        range: &RangeChip<F>,
        modulus: u64,
    ) {
        // loop over all the coefficients of the polynomial and check that they are in the modulus field
        for coeff in self.assigned_coefficients.iter() {
            range.check_less_than_safe(ctx, *coeff, modulus);
        }
    }
}
