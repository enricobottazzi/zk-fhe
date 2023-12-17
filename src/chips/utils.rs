use halo2_base::utils::ScalarField;
use num_bigint::BigInt;

/// Converts a BigInt to a Field Element
pub fn big_uint_to_fp<F: ScalarField>(big_uint: &BigInt) -> F {
    F::from_str_vartime(&big_uint.to_str_radix(10)[..]).unwrap()
}
