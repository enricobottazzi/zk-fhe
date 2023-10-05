use crate::chips::utils::{div_euclid, vec_assigned_to_vec_u64};
use halo2_base::gates::GateChip;
use halo2_base::gates::GateInstructions;
use halo2_base::safe_types::RangeChip;
use halo2_base::safe_types::RangeInstructions;
use halo2_base::utils::ScalarField;
use halo2_base::AssignedValue;
use halo2_base::Context;
use halo2_base::QuantumCell;

/// Build the sum of the polynomials a and b as sum of the coefficients
/// DEG is the degree of the input polynomials
/// Input polynomials are parsed as a vector of assigned coefficients [a_DEG, a_DEG-1, ..., a_1, a_0] where a_0 is the constant term
pub fn poly_add<const DEG: usize, F: ScalarField>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    b: Vec<AssignedValue<F>>,
    gate: &GateChip<F>,
) -> Vec<AssignedValue<F>> {
    // assert that the input polynomials have the same degree and this is equal to DEG
    assert_eq!(a.len() - 1, b.len() - 1);
    assert_eq!(a.len() - 1, DEG);

    let c: Vec<AssignedValue<F>> = a
        .iter()
        .zip(b.iter())
        .take(2 * (DEG) - 1)
        .map(|(&a, &b)| gate.add(ctx, a, b))
        .collect();

    // assert that the sum polynomial has degree DEG
    assert_eq!(c.len() - 1, DEG);

    c
}

/// Build the product of the polynomials a and b as dot product of the coefficients of a and b
/// Compared to `poly_mul_diff_deg`, this function assumes that the polynomials have the same degree and therefore optimizes the computation
/// DEG is the degree of the input polynomials
/// Input polynomials are parsed as a vector of assigned coefficients [a_DEG, a_DEG-1, ..., a_1, a_0] where a_0 is the constant term
pub fn poly_mul_equal_deg<const DEG: usize, F: ScalarField>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    b: Vec<AssignedValue<F>>,
    gate: &GateChip<F>,
) -> Vec<AssignedValue<F>> {
    // assert that the input polynomials have the same degree and this is equal to DEG
    assert_eq!(a.len() - 1, b.len() - 1);
    assert_eq!(a.len() - 1, DEG);

    let mut c = vec![];

    for i in 0..(2 * DEG + 1) {
        let mut coefficient_accumaltor = vec![];

        if i < (DEG + 1) {
            for a_idx in 0..=i {
                let a_coef = a[a_idx];
                let b_coef = b[(i - a_idx)];
                coefficient_accumaltor.push(gate.mul(ctx, a_coef, b_coef));
            }
        } else {
            for a_idx in (i - DEG)..=DEG {
                let a_coef = a[a_idx];
                let b_coef = b[(i - a_idx)];
                coefficient_accumaltor.push(gate.mul(ctx, a_coef, b_coef));
            }
        }

        let c_val = coefficient_accumaltor
            .iter()
            .fold(ctx.load_witness(F::zero()), |acc, x| gate.add(ctx, acc, *x));

        c.push(c_val);
    }

    // assert that the product polynomial has degree 2*DEG
    assert_eq!(c.len() - 1, 2 * DEG);

    c
}

/// Build the product of the polynomials a and b as dot product of the coefficients of a and b
/// Compared to `poly_mul_equal_deg`, this function doesn't assume that the polynomials have the same degree. Therefore the computation is less efficient.
/// Input polynomials are parsed as a vector of assigned coefficients [a_n, a_n-1, ..., a_1, a_0] where a_0 is the constant term and n is the degree of the polynomial
pub fn poly_mul_diff_deg<F: ScalarField>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    b: Vec<AssignedValue<F>>,
    gate: &GateChip<F>,
) -> Vec<AssignedValue<F>> {
    let a_deg = a.len() - 1;
    let b_deg = b.len() - 1;
    let c_deg = a_deg + b_deg;

    let mut c = vec![];

    for i in 0..=c_deg {
        let mut coefficient_accumaltor = vec![];

        for j in 0..=i {
            if j <= a_deg && (i - j) <= b_deg {
                let a_coef = a[j];
                let b_coef = b[i - j];

                // Update the accumulator
                coefficient_accumaltor.push(gate.mul(ctx, a_coef, b_coef));
            }
        }

        let c_val = coefficient_accumaltor
            .iter()
            .fold(ctx.load_witness(F::zero()), |acc, x| gate.add(ctx, acc, *x));

        c.push(c_val);
    }

    // assert that the product polynomial has degree c_deg
    assert_eq!(c.len() - 1, c_deg);

    c
}

/// Build the scalar multiplication of the polynomials a and the scalar k as scalar multiplication of the coefficients of a and k
/// DEG is the degree of the polynomial
/// Input polynomial is parsed as a vector of assigned coefficients [a_DEG, a_DEG-1, ..., a_1, a_0] where a_0 is the constant term
pub fn poly_scalar_mul<const DEG: usize, F: ScalarField>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    b: QuantumCell<F>,
    gate: &GateChip<F>,
) -> Vec<AssignedValue<F>> {
    // assert that the degree of the polynomial a is equal to DEG
    assert_eq!(a.len() - 1, DEG);

    let c: Vec<AssignedValue<F>> = a.iter().map(|&a| gate.mul(ctx, a, b)).collect();

    // assert that the product polynomial has degree DEG
    assert_eq!(c.len() - 1, DEG);

    c
}

/// Takes a polynomial represented by its coefficients in a vector and output a new polynomial reduced by applying modulo Q to each coefficient
/// DEG is the degree of the polynomial
/// Input polynomial is parsed as a vector of assigned coefficients [a_DEG, a_DEG-1, ..., a_1, a_0] where a_0 is the constant term
/// It assumes that the coefficients of the input polynomial can be expressed in at most num_bits bits
pub fn poly_reduce<const DEG: usize, const Q: u64, F: ScalarField>(
    ctx: &mut Context<F>,
    input: Vec<AssignedValue<F>>,
    range: &RangeChip<F>,
    num_bits: usize,
) -> Vec<AssignedValue<F>> {
    // Assert that degree of input polynomial is equal to the constant DEG
    assert_eq!(input.len() - 1, DEG);

    // Enforce that in_assigned[i] % Q = rem_assigned[i]
    let rem_assigned: Vec<AssignedValue<F>> = input
        .iter()
        .take(2 * DEG - 1)
        .map(|&x| range.div_mod(ctx, x, Q, num_bits).1)
        .collect();

    // assert that the reduced polynomial has degree DEG
    assert_eq!(rem_assigned.len() - 1, DEG);

    rem_assigned
}

/// Takes a polynomial `divisor` represented by its coefficients in a vector
/// Takes a cyclotomic polynomial `dividend` f(x)=x^m+1 (m is a power of 2) of the form represented by its coefficients in a vector
/// Output the remainder of the division of `dividend` by `dividend` as a vector of coefficients
/// DEG_DVD is the degree of the `dividend` polynomial
/// DEG_DVS is the degree of the `divisor` polynomial
/// Q is the modulus of the Ring
/// Input polynomials is parsed as a vector of assigned coefficients [a_DEG, a_DEG-1, ..., a_1, a_0] where a_0 is the constant term
/// Assumes that the coefficients of `dividend` are in the range [0, Q - 1]
/// Assumes that divisor is a cyclotomic polynomial with coefficients either 0 or 1
pub fn poly_divide_by_cyclo<
    const DEG_DVD: usize,
    const DEG_DVS: usize,
    const Q: u64,
    F: ScalarField,
>(
    ctx: &mut Context<F>,
    dividend: Vec<AssignedValue<F>>,
    divisor: Vec<AssignedValue<F>>,
    range: &RangeChip<F>,
) -> Vec<AssignedValue<F>> {
    // Assert that degree of dividend polynomial is equal to the constant DEG_DVD
    assert_eq!(dividend.len() - 1, DEG_DVD);
    // Assert that degree of divisor poly is equal to the constant DEG_DVS
    assert_eq!(divisor.len() - 1, DEG_DVS);

    // DEG_DVS must be strictly less than DEG_DVD
    assert!(DEG_DVS < DEG_DVD);

    // long division operation performed outside the circuit
    // Need to convert the dividend and divisor into a vector of u64
    let dividend_to_u64 = vec_assigned_to_vec_u64(&dividend);
    let divisor_to_u64 = vec_assigned_to_vec_u64(&divisor);

    let (quotient_to_u64, remainder_to_u64) =
        div_euclid::<DEG_DVD, DEG_DVS, Q>(&dividend_to_u64, &divisor_to_u64);

    // After the division, the degree of the quotient should be equal to DEG_DVD - DEG_DVS
    assert_eq!(quotient_to_u64.len() - 1, DEG_DVD - DEG_DVS);

    // Furthermore, the degree of the remainder must be strictly less than the degree of the divisor
    assert!(remainder_to_u64.len() - 1 < DEG_DVS);

    // Pad the remainder with 0s to make its degree equal to DEG_DVS - 1
    let mut remainder_to_u64 = remainder_to_u64;
    while remainder_to_u64.len() - 1 < DEG_DVS - 1 {
        remainder_to_u64.push(0);
    }

    // Now remainder must be of degree DEG_DVS - 1
    assert_eq!(remainder_to_u64.len() - 1, DEG_DVS - 1);

    // Later we need to perform the operation remainder + prod where prod is of degree DEG_DVD
    // In order to perform the operation inside the circuit we need to pad the remainder with 0s at the beginning to make its degree equal to DEG_DVD
    let mut remainder_to_u64 = remainder_to_u64;
    while remainder_to_u64.len() - 1 < DEG_DVD {
        remainder_to_u64.insert(0, 0);
    }

    // Now remainder must be of degree DEG_DVD
    assert_eq!(remainder_to_u64.len() - 1, DEG_DVD);

    // Assign the quotient and remainder to the circuit
    let mut quotient = vec![];
    let mut remainder = vec![];

    for i in 0..DEG_DVD - DEG_DVS + 1 {
        let val = F::from(quotient_to_u64[i]);
        let assigned_val = ctx.load_witness(val);
        quotient.push(assigned_val);
    }

    for i in 0..DEG_DVD + 1 {
        let val = F::from(remainder_to_u64[i]);
        let assigned_val = ctx.load_witness(val);
        remainder.push(assigned_val);
    }

    // Quotient is obtained by dividing the coefficients of the dividend by the highest degree coefficient of divisor
    // The coefficients of dividend are in the range [0, Q - 1] by assumption.
    // The leading coefficient of divisor is 1 by assumption.
    // Therefore, the coefficients of quotient have to be in the range [0, Q - 1]
    // Since the quotient is computed outside the circuit, we need to enforce this constraint
    for i in 0..(DEG_DVD - DEG_DVS + 1) {
        range.check_less_than_safe(ctx, quotient[i], Q);
    }

    // Remainder is equal to dividend - (quotient * divisor).
    // The coefficients of dividend are in the range [0, Q - 1] by assumption.
    // The coefficients of quotient are in the range [0, Q - 1] by constraint set above.
    // The coefficients of divisior are either 0, 1 by assumption of the cyclotomic polynomial.
    // It follows that the coefficients of quotient * divisor are in the range [0, Q - 1]
    // The remainder might have coefficients that are negative. In that case we add Q to them to make them positive.
    // Therefore, the coefficients of remainder are in the range [0, Q - 1]
    // Since the remainder is computed outside the circuit, we need to enforce this constraint
    for i in 0..DEG_DVS {
        range.check_less_than_safe(ctx, remainder[i], Q);
    }

    // check that quotient * divisor + remainder = dividend

    // DEGREE ANALYSIS
    // Quotient is of degree DEG_DVD - DEG_DVS
    // Divisor is of degree DEG_DVS
    // Quotient * divisor is of degree DEG_DVD
    // Remainder is of degree DEG_DVS - 1
    // Quotient * divisor + rem is of degree DEG_DVD
    // Dividend is of degree DEG_DVD

    // Perform the multiplication between quotient and divisor
    // No risk of overflowing the circuit prime here when performing multiplication across coefficients
    // The coefficients of quotient are in the range [0, Q - 1] by constraint set above.
    // The coefficients of divisor are either 0, 1 by assumption.

    // We use a polynomial multiplication algorithm that does not require the input polynomials to be of the same degree

    let prod = poly_mul_diff_deg(ctx, quotient, divisor, range.gate());

    // The coefficients of Divisor are either 0 or 1 given that this is a cyclotomic polynomial
    // The coefficients of Quotient are in the range [0, Q - 1] as per constraint set above.
    // Therefore, the coefficients of prod are in the range [0, Q - 1]

    // The degree of prod is DEG_DVD
    assert_eq!(prod.len() - 1, DEG_DVD);

    // Perform the addition between prod and remainder
    // The degree of prod is DEG_DVD
    // The degree of remainder is DEG_DVS - 1

    // prod + remainder

    // No risk of overflowing the circuit prime here when performing addition across coefficients
    // The coefficients of prod are in the range [0, Q - 1] by the operation above.
    // The coefficients of remainder are in the range [0, Q - 1] by constraint set above.

    let sum = poly_add::<DEG_DVD, F>(ctx, prod, remainder.clone(), range.gate());

    // assert that the degree of sum is DEG_DVD
    assert_eq!(sum.len() - 1, DEG_DVD);

    // The coefficients of prod are in the range [0, Q - 1]
    // The coefficients of remainder are in the range [0, Q - 1] by constraint set above.
    // Therefore, the coefficients of sum are in the range [0, 2 * (Q - 1)].
    // We can reduce the coefficients of sum modulo Q to make them in the range [0, Q - 1]

    // get the number of bits needed to represent the value of 2 * (Q - 1)
    let binary_representation = format!("{:b}", (2 * (Q - 1))); // Convert to binary (base-2)
    let num_bits = binary_representation.len();

    let sum_mod = poly_reduce::<DEG_DVD, Q, F>(ctx, sum, range, num_bits);

    // assert that the degree of sum_mod is DEG_DVD
    assert_eq!(sum_mod.len() - 1, DEG_DVD);

    // Enforce that sum_mod = dividend
    for i in 0..=DEG_DVD {
        let bool = range.gate().is_equal(ctx, sum_mod[i], dividend[i]);
        range.gate().assert_is_const(ctx, &bool, &F::from(1))
    }

    remainder
}
