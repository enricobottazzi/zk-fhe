use ark_bn254::Fr;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use clap::Parser;
use halo2_base::gates::GateChip;
use halo2_base::utils::ScalarField;
use halo2_base::AssignedValue;
use halo2_base::Context;
use halo2_scaffold::scaffold::cmd::Cli;
use halo2_scaffold::scaffold::run;
use serde::{Deserialize, Serialize};
use zk_fhe::chips::poly_operations::{poly_add, poly_mul, poly_scalar_mul};
use zk_fhe::poly_utils::poly_to_scalar_field_coeffs;

const N: u64 = 3; // degree of the polynomial

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CircuitInput<const N: u64> {
    pub a: Vec<u8>, // polynomial coefficients [a_0, a_1, ..., a_N] where a_i is the coefficient of x^i
    pub b: Vec<u8>, // polynomial coefficients [a_0, a_1, ..., a_N] where a_i is the coefficient of x^i
}

fn bfv_encryption_circuit<F: ScalarField>(
    ctx: &mut Context<F>,
    input: CircuitInput<N>,
    make_public: &mut Vec<AssignedValue<F>>,
) {
    // Assign the input polynomials to the circuit
    let a: Vec<AssignedValue<F>> = input
        .a
        .iter()
        .map(|x| {
            let result = F::from(*x as u64);
            ctx.load_witness(result)
        })
        .collect();

    let b: Vec<AssignedValue<F>> = input
        .b
        .iter()
        .map(|x| {
            let result = F::from(*x as u64);
            ctx.load_witness(result)
        })
        .collect();

    let gate = GateChip::<F>::default();

    // Perorm the addition of the polynomials a and b
    let c = poly_add::<N, F>(ctx, a, b, &gate);

    // Perform the scalar multiplication of the polynomial d by 2
    let d = poly_scalar_mul::<N, F>(ctx, c.clone(), 2, &gate);

    // Perform the multiplication of the polynomials d and c
    let e = poly_mul::<N, F>(ctx, d, c, &gate);

    // iterate over the coefficients of the polynomial e and make them public
    for i in 0..e.len() {
        make_public.push(e[i]);
    }

    // TEST:
    // Perform the same poly operations outside the circuit (using arkworks) and see if the results match
    let a_ark = DensePolynomial::<Fr>::from_coefficients_vec(
        input
            .a
            .iter()
            .map(|x| Fr::from(*x as u64))
            .collect::<Vec<Fr>>(),
    );

    let b_ark = DensePolynomial::<Fr>::from_coefficients_vec(
        input
            .b
            .iter()
            .map(|x| Fr::from(*x as u64))
            .collect::<Vec<Fr>>(),
    );

    let c_ark: DensePolynomial<Fr> = a_ark + b_ark;

    let d_ark: DensePolynomial<Fr> = &c_ark * Fr::from(2);

    let e_ark: DensePolynomial<Fr> = &d_ark * &c_ark;

    let e_coeffs = poly_to_scalar_field_coeffs::<F>(e_ark);

    // Compare the result of the circuit with the result of the multiplication
    for i in 0..e.len() {
        assert_eq!(e[i].value(), &e_coeffs[i]);
    }
}

fn main() {
    env_logger::init();

    let args = Cli::parse();

    // run different zk commands based on the command line arguments
    run(bfv_encryption_circuit, args);
}
