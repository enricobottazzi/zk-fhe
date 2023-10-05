use halo2_base::gates::GateChip;
use halo2_base::gates::GateInstructions;
use halo2_base::utils::ScalarField;
use halo2_base::AssignedValue;
use halo2_base::Context;
use halo2_base::QuantumCell;

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
