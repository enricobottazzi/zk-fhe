use halo2_base::gates::GateChip;
use halo2_base::gates::GateInstructions;
use halo2_base::utils::ScalarField;
use halo2_base::AssignedValue;
use halo2_base::Context;
use halo2_base::QuantumCell;

/// Build the sum of the polynomials a and b as sum of the coefficients
/// N is the degree of the polynomials
pub fn poly_add<const N: u64, F: ScalarField>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    b: Vec<AssignedValue<F>>,
    gate: &GateChip<F>,
) -> Vec<AssignedValue<F>> {
    // assert that the input polynomials have the same degree and this is equal to N
    assert_eq!(a.len() - 1, b.len() - 1);
    assert_eq!(a.len() - 1, N as usize);

    let c = a
        .iter()
        .zip(b.iter())
        .take(2 * (N as usize) - 1)
        .map(|(&a, &b)| gate.add(ctx, a, b))
        .collect();

    c
}

/// Build the product of the polynomials a and b as dot product of the coefficients of a and b
/// N is the degree of the polynomials
pub fn poly_mul<const N: u64, F: ScalarField>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    b: Vec<AssignedValue<F>>,
    gate: &GateChip<F>,
) -> Vec<AssignedValue<F>> {
    // assert that the input polynomials have the same degree and this is equal to N
    assert_eq!(a.len() - 1, b.len() - 1);
    assert_eq!(a.len() - 1, N as usize);

    let mut c = vec![];

    for i in 0..(2 * N + 1) {
        let mut coefficient_accumaltor = vec![];

        if i < (N + 1) {
            for a_idx in 0..=i {
                let a_coef = a[a_idx as usize];
                let b_coef = b[(i - a_idx) as usize];
                coefficient_accumaltor.push(gate.mul(ctx, a_coef, b_coef));
            }
        } else {
            for a_idx in (i - N)..=N {
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

    // assert that the product polynomial has degree 2N
    assert_eq!(c.len() - 1, 2 * (N as usize));

    c
}

/// Build the scalar multiplication of the polynomials a and the scalar k as scalar multiplication of the coefficients of a and k
/// N is the degree of the polynomial
pub fn poly_scalar_mul<const N: u64, F: ScalarField>(
    ctx: &mut Context<F>,
    a: Vec<AssignedValue<F>>,
    b: QuantumCell<F>,
    gate: &GateChip<F>,
) -> Vec<AssignedValue<F>> {
    // assert that the degree of the polynomial a is equal to N
    assert_eq!(a.len() - 1, N as usize);

    let c = a.iter().map(|&a| gate.mul(ctx, a, b)).collect();

    c
}
