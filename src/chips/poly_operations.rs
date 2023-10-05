use halo2_base::gates::GateChip;
use halo2_base::gates::GateInstructions;
use halo2_base::safe_types::RangeChip;
use halo2_base::safe_types::RangeInstructions;
use halo2_base::utils::ScalarField;
use halo2_base::AssignedValue;
use halo2_base::Context;
use halo2_base::QuantumCell;
use zk_fhe::chips::utils::div_euclid;

/// Build the sum of the polynomials a and b as sum of the coefficients
/// DEG is the degree of the input polynomials
/// Input polynomials are parsed as a vector of assigned coefficients [a_DEG, a_DEG-1, ..., a_1, a_0] where a_0 is the constant term
pub fn poly_add<const DEG: u64, F: ScalarField>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    b: Vec<AssignedValue<F>>,
    gate: &GateChip<F>,
) -> Vec<AssignedValue<F>> {
    // assert that the input polynomials have the same degree and this is equal to DEG
    assert_eq!(a.len() - 1, b.len() - 1);
    assert_eq!(a.len() - 1, DEG as usize);

    let c: Vec<AssignedValue<F>> = a
        .iter()
        .zip(b.iter())
        .take(2 * (DEG as usize) - 1)
        .map(|(&a, &b)| gate.add(ctx, a, b))
        .collect();

    // assert that the sum polynomial has degree DEG
    assert_eq!(c.len() - 1, DEG as usize);

    c
}

/// Build the product of the polynomials a and b as dot product of the coefficients of a and b
/// Compared to `poly_mul_diff_deg`, this function assumes that the polynomials have the same degree and therefore optimizes the computation
/// DEG is the degree of the input polynomials
/// Input polynomials are parsed as a vector of assigned coefficients [a_DEG, a_DEG-1, ..., a_1, a_0] where a_0 is the constant term
pub fn poly_mul_equal_deg<const DEG: u64, F: ScalarField>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    b: Vec<AssignedValue<F>>,
    gate: &GateChip<F>,
) -> Vec<AssignedValue<F>> {
    // assert that the input polynomials have the same degree and this is equal to DEG
    assert_eq!(a.len() - 1, b.len() - 1);
    assert_eq!(a.len() - 1, DEG as usize);

    let mut c = vec![];

    for i in 0..(2 * DEG + 1) {
        let mut coefficient_accumaltor = vec![];

        if i < (DEG + 1) {
            for a_idx in 0..=i {
                let a_coef = a[a_idx as usize];
                let b_coef = b[(i - a_idx) as usize];
                coefficient_accumaltor.push(gate.mul(ctx, a_coef, b_coef));
            }
        } else {
            for a_idx in (i - DEG)..=DEG {
                let a_coef = a[a_idx as usize];
                let b_coef = b[(i - a_idx) as usize];
                coefficient_accumaltor.push(gate.mul(ctx, a_coef, b_coef));
            }
        }

        let c_val = coefficient_accumaltor
            .iter()
            .fold(ctx.load_witness(F::zero()), |acc, x| gate.add(ctx, acc, *x));

        c.push(c_val);
    }

    // assert that the product polynomial has degree 2*DEG
    assert_eq!(c.len() - 1, 2 * (DEG as usize));

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
/// Input polynomials is parsed as a vector of assigned coefficients [a_DEG, a_DEG-1, ..., a_1, a_0] where a_0 is the constant term
pub fn poly_scalar_mul<const DEG: u64, F: ScalarField>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    b: QuantumCell<F>,
    gate: &GateChip<F>,
) -> Vec<AssignedValue<F>> {
    // assert that the degree of the polynomial a is equal to DEG
    assert_eq!(a.len() - 1, DEG as usize);

    let c: Vec<AssignedValue<F>> = a.iter().map(|&a| gate.mul(ctx, a, b)).collect();

    // assert that the product polynomial has degree DEG
    assert_eq!(c.len() - 1, DEG as usize);

    c
}

// Takes a polynomial represented by its coefficients in a vector and output a new polynomial reduced by applying modulo Q
// Q is the modulus
// N is the degree of the polynomials
pub fn poly_reduce<const Q: u64, const N: u64, F: ScalarField>(
    ctx: &mut Context<F>,
    input: Vec<AssignedValue<F>>,
    range: &RangeChip<F>,
) -> Vec<AssignedValue<F>> {

	// Assert that degree is equal to the constant N
    assert_eq!(input.len() - 1, N as usize);

    // Assign the input polynomials to the circuit
    let in_assigned: Vec<AssignedValue<F>> = input
        .iter()
        .map(|x| {
            let result = F::from(*x as u64);
            ctx.load_witness(result)})
        .collect();

    // Enforce that in_assigned[i] % Q = rem_assigned[i]
    // coefficients of input polynomials are guaranteed to be at most 16 bits by assumption
    let rem_assigned: Vec<AssignedValue<F>> =
        in_assigned.iter().take(2 * N - 1).map(|&x| range.div_mod(ctx, x, Q, 16).1).collect();

    // assert that the reduced polynomial has degree N
    assert_eq!(rem_assigned.len() - 1, N as usize);

    rem_assigned
}

// Takes a polynomial represented by its coefficients in a vector
// and output a remainder polynomial divided by divisor polynomial f(x)=x^m+1 where m is a power of 2 (public output)
// N is the degree of the dividend polynomial
// M is the degree of the divisor polynomial
// Q is the modulus
pub fn poly_divide_by_cyclo<const N: u64, const M: u64, const Q: u64, F: ScalarField>(
    ctx: &mut Context<F>,
    dividend: Vec<AssignedValue<F>>,
    divisor: Vec<AssignedValue<F>>,
    range: &RangeChip<F>,
) -> Vec<AssignedValue<F>> {
    // Assert that degree of dividend poly is equal to the constant N
    assert_eq!(dividend.len() - 1, N as usize);
    // Assert that degree of divisor poly is equal to the constant M
    assert_eq!(divisor.len() - 1, M as usize);

    // M must be less than N
    assert!(M < N);

    let mut dividend_assigned = vec![];
    let mut divisor_assigned = vec![];

    for i in 0..N + 1 {
        let val = F::from(dividend[i]);

        let assigned_val = ctx.load_witness(val);

        dividend_assigned.push(assigned_val);
    }

    for i in 0..M + 1 {
        let val = F::from(divisor[i]);

        let assigned_val = ctx.load_witness(val);

        divisor_assigned.push(assigned_val);
    }

    // assert the correct length of the assigned polynomails
    assert_eq!(dividend_assigned.len() - 1, N as usize);
    assert_eq!(divisor_assigned.len() - 1, M as usize);

    // long division operation performed outside the circuit
    let (quotient, remainder) = div_euclid(&dividend, &divisor);

    // After the division, the degree of the quotient is N - M
    assert_eq!(quotient.len() - 1, N - M);

    // Furthermore, the degree of the remainder must be strictly less than the degree of the divisor which is M
    assert!(remainder.len() - 1 < M as usize);

    // Pad the remainder with 0s to make its degree equal to M - 1
    let mut remainder = remainder;
    while remainder.len() - 1 < M - 1 {
        remainder.push(0);
    }

    // Now remainder must be of degree M - 1
    assert_eq!(remainder.len() - 1, M - 1);

    // Assign the quotient and remainder to the circuit
    let mut quot_assigned = vec![];
    let mut remainder_assigned = vec![];

    for i in 0..N - M + 1 {
        let val = F::from(quotient[i]);
        let assigned_val = ctx.load_witness(val);
        quot_assigned.push(assigned_val);
    }

    for i in 0..M {
        let val = F::from(remainder[i]);
        let assigned_val = ctx.load_witness(val);
        remainder_assigned.push(assigned_val);
    }

    // Quotient is obtained by dividing the coefficients of the dividend by the highest degree coefficient of divisor
    // The coefficients of dividend are in the range [0, Q - 1] by assumption.
    // The leading coefficient of divisor is 1 by assumption.
    // Therefore, the coefficients of quotient have to be in the range [0, Q - 1]
    // Since the quotient is computed outside the circuit, we need to enforce this constraint
    for i in 0..(N - M + 1) as usize {
        range.check_less_than_safe(ctx, quot_assigned[i], Q);
    }

    // Remainder is equal to dividend - (quotient * divisor).
    // The coefficients of dividend are in the range [0, Q - 1] by assumption.
    // The coefficients of quotient are in the range [0, Q - 1] by constraint set above.
    // The coefficients of divisior are either 0, 1 by assumption.
    // It follows that the coefficients of quotient * divisor are in the range [0, Q - 1]
    // The remainder might have coefficients that are negative. In that case we add Q to them to make them positive.
    // Therefore, the coefficients of remainder are in the range [0, Q - 1]
    // Since the remainder is computed outside the circuit, we need to enforce this constraint
    for i in 0..M as usize {
        range.check_less_than_safe(ctx, remainder_assigned[i], Q);
    }

    // ---- constraint check -----
    // check that quotient * divisor + rem = dividend

    // DEGREE ANALYSIS
    // Quotient is of degree N - M
    // Divisor is of degree M
    // Quotient * divisor is of degree N
    // Rem is of degree M - 1
    // Quotient * divisor + rem is of degree N
    // Dividend is of degree N

    // Perform the multiplication between quotient and divisor
    // We use a polynomial multiplication algorithm that does not require the input polynomials to be of the same degree

    // Initialize the resulting polynomial prod
    // The degree of prod is N
    let mut prod = vec![];

    // degree of quotient
    let a = (N - M) as usize;
    // degree of divisor
    let b = M as usize;

    // No risk of overflowing the prime here when performing multiplication across coefficients
    // The coefficients of quotient are in the range [0, Q - 1] by constraint set above.
    // The coefficients of divisor are either 0, 1 by assumption.

    // Calculate the coefficients of the product polynomial
    for i in 0..=a + b {
        // Init the accumulator for the i-th coefficient of c
        let mut coefficient_accumulator = vec![];

        for j in 0..=i {
            if j <= a && (i - j) <= b {
                let quotient_coef = quot_assigned[j];
                let divisor_coef = divisor_assigned[i - j];

                // Update the accumulator
                coefficient_accumulator.push(range.gate().mul(ctx, quotient_coef, divisor_coef));
            }
        }

        let val = coefficient_accumulator
            .iter()
            .fold(ctx.load_witness(F::zero()), |acc, x| range.gate().add(ctx, acc, *x));

        // Assign the accumulated value to the i-th coefficient of c
        prod.push(val);
    }

    // assert that the degree of prod is N
    assert_eq!(prod.len() - 1, N as usize);

    // The coefficients of Divisor are either 0 or 1 given that this is a cyclotomic polynomial
    // The coefficients of Quotient are in the range [0, Q - 1] as per constraint set above.
    // Therefore, the coefficients of prod are in the range [0, Q - 1]

    // Perform the addition between prod and remainder_assigned
    // The degree of prod is N
    // The degree of remainder_assigned is M - 1

    // We can pad the remainder_assigned with 0s to make its degree equal to N
    let mut remainder_assigned = remainder_assigned;
    while remainder_assigned.len() - 1 < N as usize {
        remainder_assigned.push(ctx.load_witness(F::zero()));
    }

    // Now remainder_assigned is of degree N
    assert_eq!(remainder_assigned.len() - 1, N as usize);

    // prod + rem_assigned

    // No risk of overflowing the prime here when performing addition across coefficients
    // The coefficients of prod are in the range [0, Q - 1] by the operation above.
    // The coefficients of rem_assigned are in the range [0, Q - 1] by constraint set above.

    let sum: Vec<AssignedValue<F>> = prod
        .iter()
        .zip(remainder_assigned.iter())
        .take(2 * N as usize - 1)
        .map(|(&a, &b)| range.gate().add(ctx, a, b))
        .collect();

    // assert that the degree of sum is N
    assert_eq!(sum.len() - 1, N as usize);

    // The coefficients of prod are in the range [0, Q - 1]
    // The coefficients of rem_assigned are in the range [0, Q - 1] by constraint set above.
    // Therefore, the coefficients of sum are in the range [0, 2 * (Q - 1)].
    // We can reduce the coefficients of sum modulo Q to make them in the range [0, Q - 1]

    // Reduce the values of sum modulo Q
    let mut sum_mod = vec![];

    // get the number of bits needed to represent the value of 2 * (Q - 1)
    let binary_representation = format!("{:b}", (2 * (Q - 1))); // Convert to binary (base-2)
    let num_bits = binary_representation.len();

    // `div_mod` assumes that sum[i] lives in q_bits which is true given the above analysis
    for i in 0..=N as usize {
        let assigned_val = range.div_mod(ctx, sum[i], Q, num_bits).1;
        sum_mod.push(assigned_val);
    }

    // assert that the degree of sum_mod is N
    assert_eq!(sum_mod.len() - 1, N as usize);

    // Enforce that sum_mod = dividend_assigned
    for i in 0..=N as usize{
        range.gate().is_equal(ctx, sum_mod[i], dividend_assigned[i]);
    }

    // ---- constraint check -----

    remainder_assigned
}

