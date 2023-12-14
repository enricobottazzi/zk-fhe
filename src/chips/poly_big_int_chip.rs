use halo2_base::safe_types::GateInstructions;
use halo2_base::utils::biguint_to_fe;
use halo2_base::{safe_types::RangeChip, utils::ScalarField, Context};
use halo2_ecc::fields::fp::FpChip;
use halo2_ecc::fields::FieldChip;
use halo2_ecc::{bigint::ProperCrtUint, fields::PrimeField};
use num_bigint::BigUint;
use num_traits::Num;

// PolyBigIntChip supports operations on Polynomials where each coefficient is defined in the Wrong Field.
// The polynomial is defined in the cyclotomic ring R_q = Z_q[X]/(X^N + 1) where q is the "wrong" field modulus of the ring R_q (ciphertext space)
#[derive(Clone, Debug)]
pub struct PolyBigIntChip<'range, F: ScalarField, const N: usize, Q: PrimeField> {
    fp_chip: FpChip<'range, F, Q>,
}

impl<'range, F: ScalarField, const N: usize, Q: PrimeField> PolyBigIntChip<'range, F, N, Q> {
    pub fn new(range: &'range RangeChip<F>, limb_bits: usize, num_limbs: usize) -> Self {
        let fp_chip = FpChip::new(range, limb_bits, num_limbs);
        Self { fp_chip }
    }

    pub fn assign_poly_in_ring(
        &self,
        ctx: &mut Context<F>,
        poly: &[BigUint],
    ) -> Vec<ProperCrtUint<F>> {
        let mut output = vec![];

        for coeff in poly.iter().take(N) {
            let coeff = biguint_to_fe::<Q>(coeff);
            let loaded = self.fp_chip.load_private(ctx, coeff);
            output.push(loaded);
        }

        output
    }

    pub fn assing_cyclotomic_poly(
        &self,
        ctx: &mut Context<F>,
        poly: &[BigUint],
    ) -> Vec<ProperCrtUint<F>> {
        let mut output = vec![];

        for coeff in poly.iter().take(N + 1) {
            let coeff = biguint_to_fe::<Q>(coeff);
            let loaded = self.fp_chip.load_private(ctx, coeff);
            output.push(loaded);
        }

        output
    }

    /// Enforce that polynomial a of degree N is sampled from the distribution chi key
    ///
    /// * Namely, that the coefficients are in the range [0, 1, Q-1].
    ///
    /// Assumption: `a` is of degree N
    pub fn check_poly_from_distribution_chi_key(
        &self,
        ctx: &mut Context<F>,
        a: &Vec<ProperCrtUint<F>>,
    ) {
        // assert that the degree of the polynomial a is equal to N - 1
        assert_eq!(a.len() - 1, N - 1);

        // In order to check that coeff is equal to either 0, 1 or q-1
        // The constraint that we want to enforce is:
        // (coeff - 0) * (coeff - 1) * (coeff - (q-1)) = 0

        // load constants
        let zero_loaded_const = self.fp_chip.load_constant_uint(ctx, BigUint::from(0_u64));

        let one_loaded_cons = self.fp_chip.load_constant_uint(ctx, BigUint::from(1_u64));

        let q_minus_one = BigUint::from_str_radix(&Q::MODULUS[2..], 16).unwrap() - 1_u64;
        let q_minus_one_loaded_const = self.fp_chip.load_constant_uint(ctx, q_minus_one);

        // loop over all the coefficients of the polynomial
        for coeff in a {
            // constrain (a - 0)
            let factor_1 = self
                .fp_chip
                .sub_no_carry(ctx, coeff.clone(), zero_loaded_const.clone());

            // constrain (a - 1)
            let factor_2 = self
                .fp_chip
                .sub_no_carry(ctx, coeff.clone(), one_loaded_cons.clone());

            // constrain (a - (q-1))
            let factor_3 =
                self.fp_chip
                    .sub_no_carry(ctx, coeff.clone(), q_minus_one_loaded_const.clone());

            // constrain (a - 0) * (a - 1)
            let factor_1_2 = self.fp_chip.mul_no_carry(ctx, factor_1, factor_2);

            // constrain (a - 0) * (a - 1) * (a - (q-1))
            let factor_1_2_3 = self.fp_chip.mul_no_carry(ctx, factor_1_2, factor_3);

            // constrain (a - 0) * (a - 1) * (a - (q-1)) = 0
            // `is_zero` requires to use a `ProperCrtUint` as input
            let factor_1_2_3_proper = self.fp_chip.carry_mod(ctx, factor_1_2_3);
            let bool = self.fp_chip.is_zero(ctx, factor_1_2_3_proper);

            self.fp_chip.gate().assert_is_const(ctx, &bool, &F::from(1));
        }
    }
}
