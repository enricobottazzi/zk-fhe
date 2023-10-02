use halo2_base::gates::GateChip;
use halo2_base::safe_types::GateInstructions;
use halo2_base::safe_types::RangeChip;
use halo2_base::safe_types::RangeInstructions;
use halo2_base::utils::ScalarField;
use halo2_base::AssignedValue;
use halo2_base::Context;
use halo2_base::QuantumCell::Constant;

/// Enforce that polynomial a of degree N is sampled from the distribution chi error
/// Namely, that the coefficients are in the range [0, B] OR [Q-B, Q-1]
/// N is the degree of the polynomial
pub fn check_poly_from_distribution_chi_error<
    const N: u64,
    const Q: u64,
    const B: u64,
    F: ScalarField,
>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    range: &RangeChip<F>,
) {
    // assert that the degree of the polynomial a is equal to N
    assert_eq!(a.len() - 1, N as usize);

    // The goal is to check that coeff is in the range [0, B] OR [Q-B, Q-1]
    // We split this check into two checks:
    // - Check that coeff is in the range [0, B] and store the boolean result in in_partial_range_1_vec
    // - Check that coeff is in the range [Q-B, Q-1] and store the boolean result in in_partial_range_2_vec
    // We then perform (`in_partial_range_1` OR `in_partial_range_2`) to check that coeff is in the range [0, B] OR [Q-B, Q-1]
    // The result of this check is stored in the `in_range` vector.
    // The bool value of `in_range` is then enforced to be true
    let mut in_range_vec = Vec::with_capacity((N + 1) as usize);

    // get the number of bits needed to represent the value of Q
    let binary_representation = format!("{:b}", Q); // Convert to binary (base-2)
    let q_bits = binary_representation.len();

    for coeff in &a {
        // Check for the range [0, B]
        let in_partial_range_1 = range.is_less_than(ctx, *coeff, Constant(F::from(B + 1)), q_bits);

        // Check for the range [Q-B, Q-1]
        let not_in_range_lower_bound =
            range.is_less_than(ctx, *coeff, Constant(F::from(Q - B)), q_bits);
        let in_range_lower_bound = range.gate.not(ctx, not_in_range_lower_bound);
        let in_range_upper_bound = range.is_less_than(ctx, *coeff, Constant(F::from(Q)), q_bits);
        let in_partial_range_2 = range
            .gate
            .and(ctx, in_range_lower_bound, in_range_upper_bound);

        // Combined check for [0, b] OR [q-b, q-1]
        let in_range = range.gate.or(ctx, in_partial_range_1, in_partial_range_2);
        in_range_vec.push(in_range);
    }

    // Enforce that in_range_vec[i] = true
    for in_range in in_range_vec {
        let bool = range.gate.is_equal(ctx, in_range, Constant(F::from(1)));
        range.gate.assert_is_const(ctx, &bool, &F::from(1));
    }
}

/// Enforce that polynomial a of degree N is sampled from the distribution chi key
/// Namely, that the coefficients are in the range [0, 1, Q-1].
/// N is the degree of the polynomial
pub fn check_poly_from_distribution_chi_key<const N: u64, const Q: u64, F: ScalarField>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    gate: &GateChip<F>,
) {
    // assert that the degree of the polynomial a is equal to N
    assert_eq!(a.len() - 1, N as usize);

    // In order to check that coeff is equal to either 0, 1 or q-1
    // The constraint that we want to enforce is:
    // (coeff - 0) * (coeff - 1) * (coeff - (q-1)) = 0

    // loop over all the coefficients of the polynomial
    for coeff in &a {
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
