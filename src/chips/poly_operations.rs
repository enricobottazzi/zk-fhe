use crate::chips::utils::big_uint_to_fp;
use crate::chips::PolyWithLength;
use axiom_eth::rlp::rlc::RlcChip;
use axiom_eth::Field;
use halo2_base::{
    gates::{GateChip, GateInstructions},
    safe_types::{RangeChip, RangeInstructions},
    AssignedValue, Context,
    QuantumCell::*,
};
use num_bigint::BigInt;

/// Build the sum of the polynomials a and b as sum of the coefficients
///
/// * DEG is the degree of the input polynomials
/// * Input polynomials are parsed as a vector of assigned coefficients [a_DEG, a_DEG-1, ..., a_1, a_0] where a_0 is the constant term
/// * Assumption: the coefficients are constrained such to avoid overflow during the polynomial addition
pub fn poly_add<const DEG: usize, F: Field>(
    ctx: &mut Context<F>,
    a: &Vec<AssignedValue<F>>,
    b: &Vec<AssignedValue<F>>,
    gate: &GateChip<F>,
) -> Vec<AssignedValue<F>> {
    // assert that the input polynomials have the same degree and this is equal to DEG
    assert_eq!(a.len() - 1, b.len() - 1);
    assert_eq!(a.len() - 1, DEG);

    let mut c = vec![];

    for i in 0..=DEG {
        let val = gate.add(ctx, a[i], b[i]);
        c.push(val);
    }

    // assert that the sum polynomial has degree DEG
    assert_eq!(c.len() - 1, DEG);

    c
}

/// Enforce that a * b = c
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
pub fn constrain_poly_mul<F: Field>(
    a_assigned_with_length: PolyWithLength<F>,
    b_assigned_with_length: PolyWithLength<F>,
    c_assigned_with_length: PolyWithLength<F>,
    ctx_gate: &mut Context<F>,
    ctx_rlc: &mut Context<F>,
    rlc: &RlcChip<F>,
) {
    // `compute_rlc` evaluates the polynomial at gamma and returns the evaluation
    let poly_a_trace = rlc.compute_rlc_fixed_len(ctx_rlc, a_assigned_with_length.assigned_poly);
    let poly_a_eval_assigned = poly_a_trace.rlc_val;

    let poly_b_trace = rlc.compute_rlc_fixed_len(ctx_rlc, b_assigned_with_length.assigned_poly);
    let poly_b_eval_assigned = poly_b_trace.rlc_val;

    let poly_c_trace = rlc.compute_rlc_fixed_len(ctx_rlc, c_assigned_with_length.assigned_poly);
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

/// Build the scalar multiplication of the polynomials a and the scalar k as scalar multiplication of the coefficients of a and k
///
/// * DEG is the degree of the polynomial
/// * Input polynomial is parsed as a vector of assigned coefficients [a_DEG, a_DEG-1, ..., a_1, a_0] where a_0 is the constant term
/// * Assumption: the coefficients are constrained such to avoid overflow during the polynomial scalar multiplication
pub fn poly_scalar_mul<const DEG: usize, F: Field>(
    ctx: &mut Context<F>,
    a: &Vec<AssignedValue<F>>,
    b: &AssignedValue<F>,
    gate: &GateChip<F>,
) -> Vec<AssignedValue<F>> {
    // assert that the degree of the polynomial a is equal to DEG
    assert_eq!(a.len() - 1, DEG);

    let mut c = vec![];

    for item in a.iter().take(DEG + 1) {
        let val = gate.mul(ctx, *item, *b);
        c.push(val);
    }

    // assert that the product polynomial has degree DEG
    assert_eq!(c.len() - 1, DEG);

    c
}

/// Takes a polynomial represented by its coefficients in a vector and output a new polynomial reduced by applying modulo Q to each coefficient
///
/// * DEG is the degree of the polynomial
/// * Input polynomial is parsed as a vector of assigned coefficients [a_DEG, a_DEG-1, ..., a_1, a_0] where a_0 is the constant term
/// * Assumption: the coefficients of the input polynomial can be expressed in at most num_bits bits
pub fn poly_reduce_by_modulo_q<const DEG: usize, const Q: u64, F: Field>(
    ctx: &mut Context<F>,
    input: &Vec<AssignedValue<F>>,
    range: &RangeChip<F>,
    num_bits: usize,
) -> Vec<AssignedValue<F>> {
    // Assert that degree of input polynomial is equal to the constant DEG
    assert_eq!(input.len() - 1, DEG);

    let mut rem_assigned = vec![];

    // Enforce that in_assigned[i] % Q = rem_assigned[i]
    for i in 0..=DEG {
        let rem = range.div_mod(ctx, input[i], Q, num_bits).1;
        rem_assigned.push(rem);
    }

    // assert that the reduced polynomial has degree DEG
    assert_eq!(rem_assigned.len() - 1, DEG);

    rem_assigned
}

/// Takes a polynomial represented by its coefficients in a vector of BigInt and output the same polynomial represented by its coefficients in a vector of assigned values of length DEG + 1
pub fn poly_big_int_assign<const DEG: usize, F: Field>(
    ctx: &mut Context<F>,
    poly: &Vec<BigInt>,
) -> Vec<AssignedValue<F>> {
    // assert that the degree of the input polynomial is equal to DEG
    assert_eq!(poly.len() - 1, DEG);

    let mut output = vec![];

    for item in poly.iter().take(DEG + 1) {
        let val = big_uint_to_fp(item);
        let assigned_val = ctx.load_witness(val);
        output.push(assigned_val);
    }

    // assert that the degree of the output polynomial is equal to DEG
    assert_eq!(output.len() - 1, DEG);

    output
}

/// Takes as inputs:
/// - a polynomial `poly` to be reduced (dividend)
/// - a `cyclo` polynomial (divisor)
/// - the `quotient` polynomial (which is the result of dividing `poly` by `cyclo`)
/// - the `quotient_times_cyclo` polynomial (which is the result of multiplying `quotient` by `cyclo`)
/// - the `remainder` polynomial of the division
/// Note that the `remainder` polynomial is the reduced version of `poly` by `cyclo`.
/// Enforce that `poly` = `quotient` * `cyclo` + `remainder`
///
/// * DEG_DVD is the degree of the dividend polynomial (poly)
/// * DEG_DVS is the degree of the divisor polynomial (cyclo)
/// * Q is the modulus of the ring R_q (cipher text space)
///
/// * Assumption: the coefficients of quotient have to be in the range [0, Q - 1]
/// * Assumption: the coefficients of remainder have to be in the range [0, Q - 1]
/// * Assumption: the coefficients are constrained such to avoid overflow during the polynomial addition between `quotient_times_cyclo` and `remainder`
pub fn constrain_poly_reduction_by_cyclo<
    const DEG_DVD: usize,
    const DEG_DVS: usize,
    const Q: u64,
    F: Field,
>(
    poly: &Vec<AssignedValue<F>>,
    cyclo: PolyWithLength<F>,
    quotient: PolyWithLength<F>,
    quotient_times_cyclo: PolyWithLength<F>,
    remainder: PolyWithLength<F>,
    range: &RangeChip<F>,
    ctx_gate: &mut Context<F>,
    ctx_rlc: &mut Context<F>,
    rlc: &RlcChip<F>,
) {
    // Note: these are not constraints but just sanity check
    assert_eq!(poly.len() - 1, DEG_DVD);
    assert_eq!(cyclo.assigned_poly.len() - 1, DEG_DVS);
    assert_eq!(poly.len() - 1, (2 * DEG_DVS) - 2);
    assert_eq!(quotient.assigned_poly.len() - 1, DEG_DVD - DEG_DVS);
    assert_eq!(quotient_times_cyclo.assigned_poly.len() - 1, DEG_DVD);
    assert_eq!(remainder.assigned_poly.len() - 1, DEG_DVD);

    // DEG_DVS must be strictly less than DEG_DVD
    assert!(DEG_DVS < DEG_DVD);

    // Check that poly = quotient * cyclo + remainder

    // First, let's constrain that quotient * cyclo = quotient_times_cyclo
    constrain_poly_mul(
        quotient,
        cyclo,
        quotient_times_cyclo.clone(),
        ctx_gate,
        ctx_rlc,
        rlc,
    );

    // Perform the addition between quotient_times_cyclo and remainder
    let sum = poly_add::<DEG_DVD, F>(
        ctx_gate,
        quotient_times_cyclo.get_poly(),
        remainder.get_poly(),
        range.gate(),
    );

    // We can reduce the coefficients of sum modulo Q to make them in the range [0, Q - 1]
    // get the number of bits needed to represent the value of (Q-1) * (DEG_DVD - DEG_DVS + 1)] + Q-1
    let binary_representation = format!(
        "{:b}",
        (Q - 1) * (DEG_DVD as u64 - DEG_DVS as u64 + 1) + (Q - 1)
    ); // Convert to binary (base-2)
    let num_bits = binary_representation.len();

    // The coefficients of sum are in the range [0, (Q-1) * (DEG_DVD - DEG_DVS + 1)] + Q-1] according to the polynomial addition constraint set above. Proof left as an exercise to the reader.
    // Therefore the coefficients of sum are known to have <= `num_bits` bits, therefore they satisfy the assumption of the `poly_reduce` chip
    let sum_mod = poly_reduce_by_modulo_q::<DEG_DVD, Q, F>(ctx_gate, &sum, range, num_bits);

    // Enforce that sum_mod = poly
    for i in 0..=DEG_DVD {
        let bool = range.gate().is_equal(ctx_gate, sum_mod[i], poly[i]);
        range.gate().assert_is_const(ctx_gate, &bool, &F::from(1))
    }
}
