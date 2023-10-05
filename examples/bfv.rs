use std::env::var;
use std::vec;

use clap::Parser;
use halo2_base::safe_types::RangeChip;
use halo2_base::safe_types::RangeInstructions;
use halo2_base::utils::ScalarField;
use halo2_base::AssignedValue;
use halo2_base::Context;
use halo2_base::QuantumCell::Constant;
use halo2_scaffold::scaffold::cmd::Cli;
use halo2_scaffold::scaffold::run;
use serde::{Deserialize, Serialize};
use zk_fhe::chips::poly_distribution::{
    check_poly_from_distribution_chi_error, check_poly_from_distribution_chi_key,
};
use zk_fhe::chips::poly_operations::{
    poly_add, poly_divide_by_cyclo, poly_mul_equal_deg, poly_reduce, poly_scalar_mul,
};

/// Circuit inputs for BFV encryption operations
///
/// # Type Parameters
///
/// * `DEG`: Degree of the cyclotomic polynomial `cyclo` of the polynomial ring R_q.
/// * `Q`: Modulus of the cipher text field
/// * `T`: Modulus of the plaintext field
/// * `DELTA` : Q/T rounded to the lower integer
/// * `B`: Upper bound of the Gaussian distribution Chi Error. It is defined as 6 * sigma
///
/// # Fields
///
/// * `pk0`: Public key 0 polynomial coefficients of degree DEG-1 [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term
/// * `pk1`: Public key 1 polynomial coefficients of degree DEG-1 [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term
/// * `m`: Plaintext polynomial of degree DEG-1 [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term. Represents the message to be encrypted
/// * `u`: Ephemeral key polynomial coefficients from the distribution ChiKey [a_DEG-1, a_DEG-2, ..., a_1, a_0]
/// * `e0`: Error polynomial coefficients from the distribution ChiError [a_DEG-1, a_DEG-2, ..., a_1, a_0]
/// * `e1`: Error polynomial coefficients from the distribution ChiError [a_DEG-1, a_DEG-2, ..., a_1, a_0]

/// # Assumes that the following checks have been performed outside the circuit
/// - `DEG` must be a power of 2
/// - `Q` must be a prime number
/// - `Q` must be greater than 1.
/// -  If n is the number of bits of Q, and m is the number of bits of the prime field of the circuit. n must be set such that (n * 2) + 2 < m to avoid overflow of the coefficients of the polynomials
/// - `T` must be a prime number and must be greater than 1 and less than `Q`
/// - `B` must be a positive integer
/// - `pk0` and `pk1` must be polynomials in the R_q ring. The ring R_q is defined as R_q = Z_q[x]/(x^DEG + 1)
/// - `cyclo` must be the cyclotomic polynomial of degree `DEG` => x^DEG + 1

// DEG and Q Parameters of the BFV encryption scheme chosen according to TABLES of RECOMMENDED PARAMETERS for 128-bits security level
// https://homomorphicencryption.org/wp-content/uploads/2018/11/HomomorphicEncryptionStandardv1.1.pdf
// B is the upper bound of the distribution Chi Error. We pick standard deviation 𝜎 ≈ 3.2 according to the HomomorphicEncryptionStandardv1 paper.
// T has been picked according to Lattigo (https://github.com/tuneinsight/lattigo/blob/master/bfv/params.go) implementation
// As suggest by https://eprint.iacr.org/2021/204.pdf (paragraph 2) we take B = 6σerr
const DEG: usize = 4;
const Q: u64 = 4637;
const T: u64 = 7;
const B: u64 = 18;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CircuitInput<const DEG: usize, const Q: u64, const T: u64, const B: u64> {
    pub pk0: Vec<u64>, // polynomial coefficients [a_N-1, a_N-2, ..., a_1, a_0] where a_0 is the constant term. Should live in R_q (to be checked outside the circuit)
    pub pk1: Vec<u64>, // polynomial coefficients [a_N-1, a_N-2, ..., a_1, a_0] where a_0 is the constant term. Should live in R_q (to be checked outside the circuit)
    pub m: Vec<u64>, // polynomial coefficients [a_N-1, a_N-2, ..., a_1, a_0] where a_0 is the constant term. Should in R_t (checked inside the circuit)
    pub u: Vec<u64>, // polynomial coefficients [a_N-1, a_N-2, ..., a_1, a_0] where a_0 is the constant term. Lives in R_q (checked inside the circuit)
    pub e0: Vec<u64>, // polynomial coefficients [a_N-1, a_N-2, ..., a_1, a_0] where a_0 is the constant term. Lives in R_q (checked inside the circuit)
    pub e1: Vec<u64>, // polynomial coefficients [a_N-1, a_N-2, ..., a_1, a_0] where a_0 is the constant term. Lives in R_q (checked inside the circuit)
    pub c0: Vec<u64>, // polynomial coefficients [a_N-1, a_N-2, ..., a_1, a_0] where a_0 is the constant term. It is compared to the ciphertext c0 generated as output by the circuit
    pub c1: Vec<u64>, // polynomial coefficients [a_N-1, a_N-2, ..., a_1, a_0] where a_0 is the constant term. It is compared to the ciphertext c1 generated as output by the circuit
}

fn bfv_encryption_circuit<F: ScalarField>(
    ctx: &mut Context<F>,
    input: CircuitInput<DEG, Q, T, B>,
    make_public: &mut Vec<AssignedValue<F>>,
) {
    // assert that the input polynomials have the same degree and this is equal to DEG - 1
    assert_eq!(input.pk0.len() - 1, DEG - 1);
    assert_eq!(input.pk1.len() - 1, DEG - 1);
    assert_eq!(input.m.len() - 1, DEG - 1);
    assert_eq!(input.u.len() - 1, DEG - 1);
    assert_eq!(input.e0.len() - 1, DEG - 1);
    assert_eq!(input.e1.len() - 1, DEG - 1);
    assert_eq!(input.c0.len() - 1, DEG - 1);
    assert_eq!(input.c1.len() - 1, DEG - 1);

    let mut pk0 = vec![];
    let mut pk1 = vec![];
    let mut u = vec![];
    let mut m = vec![];
    let mut e0 = vec![];
    let mut e1 = vec![];

    // Assign the input polynomials to the circuit
    // Using a for loop from 0 to DEG - 1 enforces that the assigned input polynomials have the same degree and this is equal to DEG - 1
    for i in 0..DEG {
        let pk0_val = F::from(input.pk0[i]);
        let pk1_val = F::from(input.pk1[i]);
        let u_val = F::from(input.u[i]);
        let m_val = F::from(input.m[i]);
        let e0_val = F::from(input.e0[i]);
        let e1_val = F::from(input.e1[i]);

        let pk0_assigned_value = ctx.load_witness(pk0_val);
        let pk1_assigned_value = ctx.load_witness(pk1_val);
        let u_assigned_value = ctx.load_witness(u_val);
        let m_assigned_value = ctx.load_witness(m_val);
        let e0_assigned_value = ctx.load_witness(e0_val);
        let e1_assigned_value = ctx.load_witness(e1_val);

        pk0.push(pk0_assigned_value);
        pk1.push(pk1_assigned_value);
        u.push(u_assigned_value);
        m.push(m_assigned_value);
        e0.push(e0_assigned_value);
        e1.push(e1_assigned_value);
    }

    assert!(pk0.len() - 1 == DEG - 1);
    assert!(pk1.len() - 1 == DEG - 1);
    assert!(u.len() - 1 == DEG - 1);
    assert!(m.len() - 1 == DEG - 1);
    assert!(e0.len() - 1 == DEG - 1);
    assert!(e1.len() - 1 == DEG - 1);

    const DELTA: u64 = Q / T; // Q/T rounded to the lower integer

    // lookup bits must agree with the size of the lookup table, which is specified by an environmental variable
    let lookup_bits = var("LOOKUP_BITS")
        .unwrap_or_else(|_| panic!("LOOKUP_BITS not set"))
        .parse()
        .unwrap();

    let range = RangeChip::default(lookup_bits);

    // Assign the cyclotomic polynomial to the circuit -> x^DEG + 1
    // Performing the assignemnt for the index 0, using a for loop from 1 to DEG - 1, and performing the assignment for the index DEG enforces that the degree of the polynomial is DEG
    let mut cyclo = vec![];

    let leading_coefficient = F::from(1);
    let leading_coefficient_assigned_value = ctx.load_witness(leading_coefficient);
    cyclo.push(leading_coefficient_assigned_value);

    for _i in 1..DEG {
        let cyclo_val = F::from(0);
        let cyclo_assigned_value = ctx.load_witness(cyclo_val);
        cyclo.push(cyclo_assigned_value);
    }

    let constant_term = F::from(1);
    let constant_term_assigned_value = ctx.load_witness(constant_term);
    cyclo.push(constant_term_assigned_value);

    assert!(cyclo.len() - 1 == DEG);

    /* Constraints on e0
        - e0 must be a polynomial in the R_q ring => Coefficients must be in the [0, Q) range and the degree of e0 must be DEG - 1
        - e0 must be sampled from the distribution ChiError

        Approach:
        - `check_poly_from_distribution_chi_error` chip guarantees that the coefficients of e0 are in the range [0, b] OR [q-b, q-1]
        - As this range is a subset of the [0, Q) range, the coefficients of e0 are in the [0, Q) range
        - The assignment for loop above guarantees that the degree of e0 is DEG - 1
    */

    /* Constraints on e1
        Same as e0
    */

    check_poly_from_distribution_chi_error::<{ DEG - 1 }, Q, B, F>(ctx, e0.clone(), &range);
    check_poly_from_distribution_chi_error::<{ DEG - 1 }, Q, B, F>(ctx, e1.clone(), &range);

    /* Constraints on u
        - u must be a polynomial in the R_q ring => Coefficients must be in the [0, Q) range and the degree of u must be DEG - 1
        - u must be sampled from the distribution ChiKey

        Approach:
        - `check_poly_from_distribution_chi_key` chip guarantees that the coefficients of u are in the range [0, 1, Q-1]
        - As this range is a subset of the [0, Q) range, the coefficients of u are in the [0, Q) range
        - The assignment for loop above guarantees that the degree of u is DEG - 1
    */

    check_poly_from_distribution_chi_key::<{ DEG - 1 }, Q, F>(ctx, u.clone(), range.gate());

    /* Constraints on m
        - m must be a polynomial in the R_t ring => Coefficients must be in the [0, T/2] OR [Q - T/2, Q - 1] range and the degree of m must be DEG - 1

        Approach:
        - Perform a range check on the coefficients of m to be in the [0, T/2] OR [Q - T/2, Q - 1] range
        - The assignment for loop above guarantees that the degree of m is DEG - 1
    */

    // TO DO: apply constraint on the coefficients of m!

    // 1. COMPUTE C0

    // pk0 * u

    // Perform the polynomial multiplication between pk0 and u.

    // OVERFLOW ANALYSIS
    // The coefficients of pk0 are in the range [0, Q) according to the check to be performed outside the circuit.
    // The coefficients of u are either [0, 1, Q-1] according to the constraints set above.
    // The maximum value of the coffiecient of pk0_u is (Q-1) * (Q-1) = Q^2 - 2Q + 1.
    // Q needs to be chosen such that Q^2 - 2Q + 1 < p where p is the prime field of the circuit in order to avoid overflow during the multiplication.

    let pk0_u = poly_mul_equal_deg::<{ DEG - 1 }, F>(ctx, pk0.clone(), u.clone(), &range.gate());

    // pk0_u is a polynomial of degree (DEG - 1) * 2 = 2*DEG - 2
    // pk0_u has coefficients in the [0, Q^2 - 2Q + 1] range
    // Reduce the coefficients by modulo `Q`

    // get the number of bits needed to represent the value of Q^2 - 2Q + 1
    let binary_representation = format!("{:b}", (Q.pow(2) - (2 * Q) + 1));
    let num_bits_1 = binary_representation.len();

    let pk0_u = poly_reduce::<{ 2 * DEG - 2 }, Q, F>(ctx, pk0_u, &range, num_bits_1);

    // pk0_u is a polynomial of degree (DEG - 1) * 2 = 2*DEG - 2
    // pk0_u has coefficients in the [0, Q) range
    // cyclo is a polynomial of degree DEG
    // Reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1 to get a polynomial of degree DEG - 1
    let pk0_u =
        poly_divide_by_cyclo::<{ 2 * DEG - 2 }, DEG, Q, F>(ctx, pk0_u, cyclo.clone(), &range);

    // assert that the degree of pk0_u is 2*DEG - 2
    assert_eq!(pk0_u.len() - 1, 2 * DEG - 2);

    // But actually, the degree of pk0_u is DEG - 1, the first DEG - 1 coefficients are just zeroes
    // Therefore, we need to trim the first DEG - 1 coefficients

    let mut pk0_u_trimmed = vec![];
    for i in DEG - 1..pk0_u.len() {
        pk0_u_trimmed.push(pk0_u[i]);
    }

    // assert that the degree of pk0_u_trimmed is DEG - 1
    assert_eq!(pk0_u_trimmed.len() - 1, DEG - 1);

    // pk0_u_trimmed is a polynomial in the R_q ring!

    // m * delta

    // Perform the polynomial scalar multiplication between m and delta.

    // OVERFLOW ANALYSIS
    // The coefficients of m are in the range [0, T/2] OR [Q - T/2, Q - 1] according to the constaints set above.
    // Delta is a constant in the range [0, Q) as it is defined as Q/T rounded to the lower integer and T < Q and T > 1.
    // The maximum value of the coffiecient of m_delta is (Q-1) * (Q-1) = Q^2 - 2Q + 1.
    // T has to be less than Q (check performed outside the circuit).
    // If the previous condition (Q^2 - 2Q + 1 < p) is satisfied there is no risk of overflow during the scalar multiplication.

    let m_delta =
        poly_scalar_mul::<{ DEG - 1 }, F>(ctx, m.clone(), Constant(F::from(DELTA)), range.gate());

    // Reduce the coefficients of `m_delta` by modulo `Q`
    // Note: Scalar multiplication does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1
    let m_delta = poly_reduce::<{ DEG - 1 }, Q, F>(ctx, m_delta, &range, num_bits_1);
    // m_delta is a polynomial in the R_q ring

    // pk0_u_trimmed + m_delta

    // Perform the polynomial addition between pk0_u_trimmed and m_delta.
    // Note: Addition does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1

    // OVERFLOW ANALYSIS
    // The coefficients of pk0_u_trimmed and m_delta are in the [0, Q) range according to the constraints set above.
    // The maximum value of the coffiecient of pk0_u_trimmed_plus_m_delta is (Q-1) + (Q-1) = 2Q - 2.
    // If the previous condition (Q^2 - 2Q + 1 < p) is satisfied there is no risk of overflow during the addition.

    let pk0_u_trimmed_plus_m_delta =
        poly_add::<{ DEG - 1 }, F>(ctx, pk0_u_trimmed, m_delta, range.gate());

    // Reduce the coefficients of `pk0_u_trimmed_plus_m_delta` by modulo `Q`
    let pk0_u_trimmed_plus_m_delta =
        poly_reduce::<{ DEG - 1 }, Q, F>(ctx, pk0_u_trimmed_plus_m_delta, &range, num_bits_1);

    // pk0_u_trimmed_plus_m_delta is a polynomial in the R_q ring

    // c0 = pk0_u_trimmed_plus_m_delta + e0

    // Perform the polynomial addition between pk0_u_trimmed_plus_m_delta and e0.

    // OVERFLOW ANALYSIS
    // The coefficients of pk0_u_trimmed_plus_m_delta are in the [0, Q) range according to the constraints set above.
    // The coefficients of e0 are in the range [0, b] OR [q-b, q-1] according to the constraints set above.
    // The maximum value of the coffiecient of c0 is (Q-1) + (Q-1) = 2Q - 2.
    // If the previous condition (Q^2 - 2Q + 1 < p) is satisfied there is no risk of overflow during the addition.

    let c0 = poly_add::<{ DEG - 1 }, F>(ctx, pk0_u_trimmed_plus_m_delta, e0, range.gate());

    // get the number of bits needed to represent the value of 2Q - 2.
    let binary_representation = format!("{:b}", 2 * Q - 2);
    let num_bits_2 = binary_representation.len();

    // Reduce the coefficients by modulo `Q`
    // Note: Addition does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1
    let c0 = poly_reduce::<{ DEG - 1 }, Q, F>(ctx, c0, &range, num_bits_2);

    // c0 is a polynomial in the R_q ring

    // 1. COMPUTE C1

    // pk1 * u

    // Perform the polynomial multiplication between pk1 and u.

    // OVERFLOW ANALYSIS
    // The coefficients of pk1 are in the range [0, Q) according to the check to be performed outside the circuit.
    // The coefficients of u are either [0, 1, Q-1] according to the constraints set above.
    // The maximum value of the coffiecient of pk1_u is (Q-1) * (Q-1) = Q^2 - 2Q + 1.
    // If the previous condition (Q^2 - 2Q + 1 < p) is satisfied there is no risk of overflow during the multiplication.

    let pk1_u = poly_mul_equal_deg::<{ DEG - 1 }, F>(ctx, pk1.clone(), u, range.gate());

    // // pk1_u is a polynomial of degree (DEG - 1) * 2 = 2*DEG - 2
    // // pk1_u has coefficients in the [0, Q^2 - 2Q + 1] range
    // // Reduce the coefficients by modulo `Q`
    let pk1_u = poly_reduce::<{ 2 * DEG - 2 }, Q, F>(ctx, pk1_u, &range, num_bits_1);

    // pk1_u is a polynomial of degree (DEG - 1) * 2 = 2*DEG - 2
    // pk1_u has coefficients in the [0, Q) range
    // cyclo is a polynomial of degree DEG
    // Reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1 to get a polynomial of degree DEG - 1
    let pk1_u =
        poly_divide_by_cyclo::<{ 2 * DEG - 2 }, DEG, Q, F>(ctx, pk1_u, cyclo.clone(), &range);

    // assert that the degree of pk1_u is 2*DEG - 2
    assert_eq!(pk1_u.len() - 1, 2 * DEG - 2);

    // But actually, the degree of pk1_u is DEG - 1, the first DEG - 1 coefficients are just zeroes
    // Therefore, we need to trim the first DEG - 1 coefficients
    let mut pk1_u_trimmed = vec![];
    for i in DEG - 1..pk1_u.len() {
        pk1_u_trimmed.push(pk1_u[i]);
    }

    // assert that the degree of pk1_u_trimmed is DEG - 1
    assert_eq!(pk1_u_trimmed.len() - 1, DEG - 1);

    // pk1_u_trimmed is a polynomial in the R_q ring

    // c1 = pk1_u_trimmed + e0

    // OVERFLOW ANALYSIS
    // The coefficients of pk1_u are in the [0, Q) range according to the constraints set above.
    // The coefficients of e1 are in the range [0, b] OR [q-b, q-1] according to the constraints set above.
    // The maximum value of the coffiecient of c1 is (Q-1) + (Q-1) = 2Q - 2.
    // If the previous condition (Q^2 - 2Q + 1 < p) is satisfied there is no risk of overflow during the addition.

    // Perform the polynomial addition between pk1_u and e1.
    let c1 = poly_add::<{ DEG - 1 }, F>(ctx, pk1_u_trimmed, e1, range.gate());

    // Reduce the coefficients by modulo `Q`
    let c1 = poly_reduce::<{ DEG - 1 }, Q, F>(ctx, c1, &range, num_bits_2);

    // Note: Addition does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1

    // c1 is a polynomial in the R_q ring

    // That that c0 and c1 computed inside the circuit are equal to the ciphertexts provided as input in the test vector json file
    // Check outside the circuit that the remainder matches the expected one
    for i in 0..DEG {
        assert_eq!(*c0[i].value(), F::from(input.c0[i]));
        assert_eq!(*c1[i].value(), F::from(input.c1[i]));
    }

    // Expose to the public the coefficients of c0 and c1
    for i in 0..DEG {
        make_public.push(c0[i]);
    }

    for i in 0..DEG {
        make_public.push(c1[i]);
    }

    // Expose to the public pk0 and pk1
    for i in 0..DEG {
        make_public.push(pk0[i]);
    }

    for i in 0..DEG {
        make_public.push(pk1[i]);
    }

    // Expose to the public `cyclo`
    for i in 0..DEG + 1 {
        make_public.push(cyclo[i]);
    }
}

fn main() {
    env_logger::init();

    let args = Cli::parse();

    run(bfv_encryption_circuit, args);
}
