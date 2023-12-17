use crate::chips::utils::big_uint_to_fp;
use crate::poly::Poly;
use axiom_eth::{rlp::rlc::RlcChip, Field};
use halo2_base::{
    AssignedValue, Context,
    QuantumCell::{Constant, Existing},
};
use num_bigint::BigInt;
use num_traits::Num;

/// Struct to perform polynomial arithmetic operations inside the circuit.
#[derive(Clone, Debug)]
pub struct PolyAssigned<const DEG: usize, F: Field>
where
    [(); DEG + 1]: Sized,
{
    pub assigned_coefficients: [AssignedValue<F>; DEG + 1],
    pub max_num_bits: usize,
}

impl<const DEG: usize, F: Field> PolyAssigned<DEG, F>
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
            let val = big_uint_to_fp(item);
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
        b: PolyAssigned<DEG, F>,
        c: PolyAssigned<{ DEG * 2 }, F>,
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
}
