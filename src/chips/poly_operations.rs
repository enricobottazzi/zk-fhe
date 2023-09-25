use halo2_base::gates::GateChip;
use halo2_base::gates::GateInstructions;
use halo2_base::halo2_proofs::halo2curves::bn256::Fr;
use halo2_base::AssignedValue;
use halo2_base::Context;

// TO DO:
// - Add assumptions
// - Add tests
// - Add documentation
// - Assert that the polynomials are of the same degree
// - Work with polynomials at the example level
// - Specify what is Q in poly lib

// Build the sum of the polynomials a and b as sum of the coefficients
fn poly_add(
    ctx: &mut Context<Fr>,
    a: Vec<AssignedValue<Fr>>,
    b: Vec<AssignedValue<Fr>>,
    gate: &GateChip<Fr>,
) -> Vec<AssignedValue<Fr>> {
    let degree = a.len() - 1;

    // assert that the polynomials are of the same degree
    assert_eq!(degree, b.len() - 1);

    let c: Vec<AssignedValue<Fr>> = a
        .iter()
        .zip(b.iter())
        .take(2 * degree - 1)
        .map(|(&a, &b)| gate.add(ctx, a, b))
        .collect();

    c
}

// Build the product of the polynomials a and b as dot product of the coefficients of a and b
fn poly_mul(
    ctx: &mut Context<Fr>,
    a: Vec<AssignedValue<Fr>>,
    b: Vec<AssignedValue<Fr>>,
    gate: &GateChip<Fr>,
) -> Vec<AssignedValue<Fr>> {
    let degree = a.len() - 1;

    // assert that the polynomials are of the same degree
    assert_eq!(degree, b.len() - 1);

    let c: Vec<AssignedValue<Fr>> = vec![gate.mul(ctx, a[0], b[0])];

    let mut prod_val: Vec<AssignedValue<Fr>> = vec![];

    for i in 0..(2 * degree + 1) {
        let mut coefficient_accumaltor: Vec<AssignedValue<Fr>> = vec![];

        if i < degree + 1 {
            for a_idx in 0..=i {
                let a_coef = a[a_idx];
                let b_coef = b[i - a_idx];
                coefficient_accumaltor.push(gate.mul(ctx, a_coef, b_coef));
            }
        } else {
            for a_idx in (i - degree)..=degree {
                let a_coef = a[a_idx];
                let b_coef = b[i - a_idx];
                coefficient_accumaltor.push(gate.mul(ctx, a_coef, b_coef));
            }
        }

        let prod_value = coefficient_accumaltor
            .iter()
            .fold(ctx.load_witness(Fr::zero()), |acc, x| {
                gate.add(ctx, acc, *x)
            });

        prod_val.push(prod_value);
    }
    c
}

// Build the scalar multiplication of the polynomials a and the scalar k as scalar multiplication of the coefficients of a and k
fn poly_scalar_mul(
    ctx: &mut Context<Fr>,
    a: Vec<AssignedValue<Fr>>,
    b: AssignedValue<Fr>,
    gate: &GateChip<Fr>,
) -> Vec<AssignedValue<Fr>> {
    let c: Vec<AssignedValue<Fr>> = a.iter().map(|&a| gate.mul(ctx, a, b)).collect();

    c
}
