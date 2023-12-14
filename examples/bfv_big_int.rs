#![feature(adt_const_params)]
use axiom_eth::Field;
use clap::Parser;
use halo2_base::halo2_proofs::halo2curves::secp256k1::Fq as Q;
use halo2_base::safe_types::RangeChip;
use halo2_base::{AssignedValue, Context};
use halo2_scaffold::scaffold::cmd::Cli;
use halo2_scaffold::scaffold::run;
use serde::{Deserialize, Serialize};
use std::env::var;
use zk_fhe::chips::poly_big_int_chip::PolyBigIntChip;
use zk_fhe::chips::utils::vec_string_to_vec_biguint;

// TO DO:
// - [] Fix the prime field to not be something that implements PrimeField
// - [] Update assumptions at the top of the file
// - [] input generated from the python script `python3 cli.py -n 4 -q 115792089237316195423570985008687907852837564279074904382605163141518161494337 -t 257 --output input.json`
// - [] run the circuit `LOOKUP_BITS=8 cargo run --example bfv_big_int -- --name bfv_3 -k 9  mock`

/// Circuit inputs for BFV encryption operations
///
/// # Type Parameters
///
/// * `DEG`: Degree of the cyclotomic polynomial `cyclo` of the polynomial ring R_q.
/// * `T`: Modulus of the plaintext field
/// * `B`: Upper bound of the Gaussian distribution Chi Error. It is defined as 6 * 𝜎
///
/// # Fields
///
/// * `pk0`: Public key 0 - polynomial of degree DEG-1 living in ciphertext space R_q
/// * `pk1`: Public key 1 - polynomial of degree DEG-1 living in ciphertext space R_q
/// * `m`: Plaintext message to be encrypted - polynomial of degree DEG-1 living in plaintext space R_t
/// * `u`: Ephemeral key - polynomial of degree DEG-1 living in ciphertext space R_q - its coefficients are sampled from the distribution ChiKey
/// * `e0`: Error - polynomial of degree DEG-1 living in ciphertext space R_q - its coefficients are sampled from the distribution ChiError
/// * `e1`: Error - polynomial of degree DEG-1 living in ciphertext space R_q - its coefficients are sampled from the distribution ChiError
/// * `c0`: First ciphertext component - polynomial of degree DEG-1 living in ciphertext space R_q
/// * `c1`: Second ciphertext component - polynomial of degree DEG-1 living in ciphertext space R_q
/// * `cyclo`: Cyclotomic polynomial of degree DEG in the form x^DEG + 1
///
/// Note: all the polynomials are expressed by their coefficients in the form [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term
///
/// # Assumptions (to be checked on the public inputs outside the circuit)
///
/// * `DEG` must be a power of 2
/// * `Q` must be a prime number and be greater than 1.
/// * `T` must be a prime number and must be greater than 1 and less than `Q`
/// * `B` must be a positive integer and must be less than `Q`
/// * `cyclo` must be the cyclotomic polynomial of degree `DEG` in the form x^DEG + 1
/// * `pk0` and `pk1` must be polynomials in the R_q ring. The ring R_q is defined as R_q = Z_q[x]/(x^DEG + 1)
/// *  Q and DEG must be chosen such that (Q-1) * (Q-1) * (DEG+1) + (Q-1) < p, where p is the modulus of the circuit field to avoid overflow during polynomial addition inside the circuit
/// *  Q, T must be chosen such that (Q-1) * (Q-T) + (Q-1) + (Q-1) < p, where p is the modulus of the circuit field.. This is required to avoid overflow during polynomial scalar multiplication inside the circuit
/// *  Q must be chosen such that 2Q - 2 < p, where p is the modulus of the circuit field. This is required to avoid overflow during polynomial addition inside the circuit

// For real world applications, the parameters should be chosen according to the security level required.
// DEG and Q Parameters of the BFV encryption scheme should be chosen according to TABLES of RECOMMENDED PARAMETERS for 128-bits security level
// https://homomorphicencryption.org/wp-content/uploads/2018/11/HomomorphicEncryptionStandardv1.1.pdf
// B is the upper bound of the distribution Chi Error. Pick standard deviation 𝜎 ≈ 3.2 according to the HomomorphicEncryptionStandardv1 paper.
// T is picked according to Lattigo (https://github.com/tuneinsight/lattigo/blob/master/schemes/bfv/example_parameters.go) implementation
// As suggest by https://eprint.iacr.org/2021/204.pdf (paragraph 2) B = 6σerr
// These are just parameters used for fast testing purpose - to match with input file `data/bfv.in`
const DEG: usize = 4;
const T: u64 = 7;
const B: u64 = 19;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CircuitInput<const DEG: usize, const T: u64, const B: u64> {
    pub pk0: Vec<String>, // PUBLIC INPUT. Should live in R_q according to assumption
    pub pk1: Vec<String>, // PUBLIC INPUT. Should live in R_q according to assumption
    pub m: Vec<String>,   // PRIVATE INPUT. Should in R_t (enforced inside the circuit)
    pub u: Vec<String>, // PRIVATE INPUT. Should live in R_q and be sampled from the distribution ChiKey (enforced inside the circuit)
    pub e0: Vec<String>, // PRIVATE INPUT. Should live in R_q and be sampled from the distribution ChiError (enforced inside the circuit)
    pub e1: Vec<String>, // PRIVATE INPUT. Should live in R_q and be sampled from the distribution ChiError (enforced inside the circuit)
    pub c0: Vec<String>, // PUBLIC INPUT. Should live in R_q. We constraint equality between c0 and computed_c0 namely the ciphertext computed inside the circuit
    pub c1: Vec<String>, // PUBLIC INPUT. Should live in R_q. We constraint equality between c1 and computed_c1 namely the ciphertext computed inside the circuit
    pub cyclo: Vec<String>, // PUBLIC INPUT. Should be the cyclotomic polynomial of degree DEG in the form x^DEG + 1 according to assumption
}

fn bfv_encryption_circuit<F: Field>(
    ctx: &mut Context<F>,
    input: CircuitInput<DEG, T, B>,
    make_public: &mut Vec<AssignedValue<F>>,
) {
    // Note: this is not a constraint enforced inside the circuit, just a sanity check
    assert_eq!(input.pk0.len() - 1, DEG - 1);
    assert_eq!(input.pk1.len() - 1, DEG - 1);
    assert_eq!(input.m.len() - 1, DEG - 1);
    assert_eq!(input.u.len() - 1, DEG - 1);
    assert_eq!(input.e0.len() - 1, DEG - 1);
    assert_eq!(input.e1.len() - 1, DEG - 1);
    assert_eq!(input.c0.len() - 1, DEG - 1);
    assert_eq!(input.c1.len() - 1, DEG - 1);
    assert_eq!(input.cyclo.len() - 1, DEG);

    // Transform the input polynomials from strings to BigUints
    let pk0_big_int = vec_string_to_vec_biguint(&input.pk0);
    let pk1_big_int = vec_string_to_vec_biguint(&input.pk1);
    let m_big_int = vec_string_to_vec_biguint(&input.m);
    let u_big_int = vec_string_to_vec_biguint(&input.u);
    let e0_big_int = vec_string_to_vec_biguint(&input.e0);
    let e1_big_int = vec_string_to_vec_biguint(&input.e1);
    let c0_big_int = vec_string_to_vec_biguint(&input.c0);
    let c1_big_int = vec_string_to_vec_biguint(&input.c1);
    let cyclo_big_int = vec_string_to_vec_biguint(&input.cyclo);

    // lookup bits must agree with the size of the lookup table, which is specified by an environmental variable
    let lookup_bits = var("LOOKUP_BITS")
        .unwrap_or_else(|_| panic!("LOOKUP_BITS not set"))
        .parse()
        .unwrap();

    // create a Range chip that contains methods for basic arithmetic operations
    let range = RangeChip::default(lookup_bits);

    // The prime, in this case, is 256 bits. Therefore, we use 4 limbs of 64 bits each to represent the prime in the circuit
    let limb_bits = 64;
    let num_limbs = 4;

    // Create Field Chip
    let poly_big_int_chip = PolyBigIntChip::<F, DEG, Q>::new(&range, limb_bits, num_limbs);

    // Assign the polynomials to the circuit
    let pk0 = poly_big_int_chip.assign_poly_in_ring(ctx, &pk0_big_int);
    let pk1 = poly_big_int_chip.assign_poly_in_ring(ctx, &pk1_big_int);
    let m = poly_big_int_chip.assign_poly_in_ring(ctx, &m_big_int);
    let u = poly_big_int_chip.assign_poly_in_ring(ctx, &u_big_int);
    let e0 = poly_big_int_chip.assign_poly_in_ring(ctx, &e0_big_int);
    let e1 = poly_big_int_chip.assign_poly_in_ring(ctx, &e1_big_int);
    let c0 = poly_big_int_chip.assign_poly_in_ring(ctx, &c0_big_int);
    let c1 = poly_big_int_chip.assign_poly_in_ring(ctx, &c1_big_int);
    let cyclo = poly_big_int_chip.assing_cyclotomic_poly(ctx, &cyclo_big_int);

    /* constraint on u
        - u must be a polynomial in the R_q ring => Coefficients must be in the [0, Q-1] range and the degree of u must be DEG - 1
        - u must be sampled from the distribution ChiKey, namely the coefficients of u must be either 0, 1 or Q-1

        Approach:
        - `check_poly_from_distribution_chi_key` chip guarantees that the coefficients of u are either 0, 1 or Q-1
        - As this range is a subset of the [0, Q-1] range, the coefficients of u are guaranteed to be in the [0, Q-1] range
        - The assignment for loop above guarantees that the degree of u is DEG - 1
    */
    poly_big_int_chip.check_poly_from_distribution_chi_key(ctx, &u);
}

fn main() {
    env_logger::init();

    let args = Cli::parse();

    run(bfv_encryption_circuit, args);
}
