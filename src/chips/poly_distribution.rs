use halo2_base::gates::GateChip;
use halo2_base::safe_types::GateInstructions;
use halo2_base::safe_types::RangeChip;
use halo2_base::safe_types::RangeInstructions;
use halo2_base::utils::ScalarField;
use halo2_base::AssignedValue;
use halo2_base::Context;
use halo2_base::QuantumCell::Constant;

/// Enforce that polynomial a of degree DEG is sampled from the distribution chi error
/// Namely, that the coefficients are in the range [0, B] OR [Q-B, Q-1]
/// DEG is the degree of the polynomial
pub fn check_poly_from_distribution_chi_error<
    const DEG: usize,
    const Q: u64,
    const B: u64,
    F: ScalarField,
>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    range: &RangeChip<F>,
) {
    // assert that the degree of the polynomial a is equal to DEG
    assert_eq!(a.len() - 1, DEG);

    // The goal is to check that coeff is in the range [0, B] OR [Q-B, Q-1]
    // We split this check into two checks:
    // - Check that coeff is in the range [0, B] and store the boolean result in in_partial_range_1_vec
    // - Check that coeff is in the range [Q-B, Q-1] and store the boolean result in in_partial_range_2_vec
    // We then perform (`in_partial_range_1` OR `in_partial_range_2`) to check that coeff is in the range [0, B] OR [Q-B, Q-1]
    // The result of this check is stored in the `in_range` vector.
    // The bool value of `in_range` is then enforced to be true
    let mut in_range_vec = Vec::with_capacity((DEG + 1));

    // get the number of bits needed to represent the value of Q
    let binary_representation = format!("{:b}", Q);
    let q_bits = binary_representation.len();

    for coeff in &a {
        // First of all, enforce that coefficients are in the [0, 2^q_bits) range to satisfy the assumption of `is_less_than` chip
        range.is_less_than_safe(ctx, *coeff, 1 << q_bits as u64);

        // Check for the range [0, B]
        // coeff is known are known to have <= `q_bits` bits according to the constraint set above
        // B + 1 is known to have <= `q_bits` bits as public constant
        // Therefore it satisfies the assumption of `is_less_than` chip
        let in_partial_range_1 = range.is_less_than(ctx, *coeff, Constant(F::from(B + 1)), q_bits);

        // Check for the range [Q-B, Q-1]
        // coeff is known are known to have <= `q_bits` bits according to the constraint set above
        // Q - B is known to have <= `q_bits` bits as public constant
        // Therefore it satisfies the assumption of `is_less_than` chip
        let not_in_range_lower_bound =
            range.is_less_than(ctx, *coeff, Constant(F::from(Q - B)), q_bits);
        let in_range_lower_bound = range.gate.not(ctx, not_in_range_lower_bound);

        // coeff is known are known to have <= `q_bits` bits according to the constraint set above
        // Q is known to have <= `q_bits` bits as public constant
        // Therefore it satisfies the assumption of `is_less_than` chip
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

/// Enforce that polynomial a of degree DEG is sampled from the distribution chi key
/// Namely, that the coefficients are in the range [0, 1, Q-1].
/// DEG is the degree of the polynomial
pub fn check_poly_from_distribution_chi_key<const DEG: usize, const Q: u64, F: ScalarField>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    gate: &GateChip<F>,
) {
    // assert that the degree of the polynomial a is equal to DEG
    assert_eq!(a.len() - 1, DEG);

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
