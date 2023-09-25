use halo2_base::gates::GateChip;
use halo2_base::gates::GateInstructions;
use halo2_base::halo2_proofs::halo2curves::bn256::Fr;
use halo2_base::AssignedValue;
use halo2_base::Context;

// Build the sum of the polynomials a and b as sum of the coefficients
// N is the degree of the polynomial
fn poly_add<const N: u64>(
    ctx: &mut Context<Fr>,
    a: Vec<AssignedValue<Fr>>,
    b: Vec<AssignedValue<Fr>>,
    gate: &GateChip<Fr>,
) -> Vec<AssignedValue<Fr>> {
    // assert that the input polynomials have the same degree and this is equal to N
    assert_eq!(a.len() - 1, b.len() - 1);
    assert_eq!(a.len() - 1, N as usize);

    let c: Vec<AssignedValue<Fr>> = a
        .iter()
        .zip(b.iter())
        .take(2 * (N as usize) - 1)
        .map(|(&a, &b)| gate.add(ctx, a, b))
        .collect();

    c
}

// Build the product of the polynomials a and b as dot product of the coefficients of a and b
fn poly_mul<const N: u64>(
    ctx: &mut Context<Fr>,
    a: Vec<AssignedValue<Fr>>,
    b: Vec<AssignedValue<Fr>>,
    gate: &GateChip<Fr>,
) -> Vec<AssignedValue<Fr>> {
    // assert that the input polynomials have the same degree and this is equal to N
    assert_eq!(a.len() - 1, b.len() - 1);
    assert_eq!(a.len() - 1, N as usize);

    let c: Vec<AssignedValue<Fr>> = vec![gate.mul(ctx, a[0], b[0])];

    let mut prod_val: Vec<AssignedValue<Fr>> = vec![];

    for i in 0..(2 * (N + 1)) {
        let mut coefficient_accumaltor: Vec<AssignedValue<Fr>> = vec![];

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
fn poly_scalar_mul<const N: u64>(
    ctx: &mut Context<Fr>,
    a: Vec<AssignedValue<Fr>>,
    b: AssignedValue<Fr>,
    gate: &GateChip<Fr>,
) -> Vec<AssignedValue<Fr>> {
    // assert that the degree of the polynomial a is equal to N
    assert_eq!(a.len() - 1, N as usize);

    let c: Vec<AssignedValue<Fr>> = a.iter().map(|&a| gate.mul(ctx, a, b)).collect();

    c
}
