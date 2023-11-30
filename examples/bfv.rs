use axiom_eth::keccak::KeccakChip;
use axiom_eth::{EthChip, Field};
use clap::Parser;
use halo2_base::safe_types::{GateInstructions, RangeInstructions};
use halo2_base::{AssignedValue, Context, QuantumCell::Constant};
use halo2_scaffold::scaffold::{cmd::Cli, run_eth};
use num_bigint::BigInt;
use num_traits::Num;
use serde::{Deserialize, Serialize};
use zk_fhe::chips::poly_distribution::{
    check_poly_coefficients_in_range, check_poly_from_distribution_chi_key,
};
use zk_fhe::chips::poly_operations::{
    constrain_poly_mul, constrain_poly_reduction_by_cyclo, poly_add, poly_big_int_assign,
    poly_reduce_by_modulo_q, poly_scalar_mul,
};
use zk_fhe::chips::utils::{
    div_euclid, poly_mul, reduce_poly_by_modulo_q, vec_string_to_vec_bigint,
};
use zk_fhe::chips::PolyWithLength;

/// Circuit inputs for BFV encryption operations
///
/// # Type Parameters
///
/// * `DEG`: Degree of the cyclotomic polynomial `cyclo` of the polynomial ring R_q.
/// * `Q`: Modulus of the cipher text field
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
/// *  Q and T must be chosen such that (Q-1) * (Q/T) < p, where p is the modulus of the circuit field. This is required to avoid overflow during polynomial scalar multiplication inside the circuit
/// *  Q must be chosen such that 2Q - 2 < p, where p is the modulus of the circuit field. This is required to avoid overflow during polynomial addition inside the circuit

// For real world applications, the parameters should be chosen according to the security level required.
// DEG and Q Parameters of the BFV encryption scheme should be chosen according to TABLES of RECOMMENDED PARAMETERS for 128-bits security level
// https://homomorphicencryption.org/wp-content/uploads/2018/11/HomomorphicEncryptionStandardv1.1.pdf
// B is the upper bound of the distribution Chi Error. Pick standard deviation 𝜎 ≈ 3.2 according to the HomomorphicEncryptionStandardv1 paper.
// T is picked according to Lattigo (https://github.com/tuneinsight/lattigo/blob/master/schemes/bfv/example_parameters.go) implementation
// As suggest by https://eprint.iacr.org/2021/204.pdf (paragraph 2) B = 6σerr
// These are just parameters used for fast testing purpose - to match with input file `data/bfv.in`
const DEG: usize = 4;
const Q: u64 = 4637;
const T: u64 = 7;
const B: u64 = 18;

// These are the parameters used for the real world application - to match with input file `data/bfv_2.in`
// const DEG: usize = 1024;
// const Q: u64 = 536870909;
// const T: u64 = 7;
// const B: u64 = 18;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CircuitInput<const DEG: usize, const Q: u64, const T: u64, const B: u64> {
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
    _eth_chip: &EthChip<F>,
    _keccak: &mut KeccakChip<F>,
    input: CircuitInput<DEG, Q, T, B>,
    make_public: &mut Vec<AssignedValue<F>>,
) -> impl FnOnce(&mut Context<F>, &mut Context<F>, &EthChip<F>) + Clone {
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

    // Transform the input polynomials from strings to BigInts
    let pk0_big_int = vec_string_to_vec_bigint(&input.pk0);
    let pk1_big_int = vec_string_to_vec_bigint(&input.pk1);
    let m_big_int = vec_string_to_vec_bigint(&input.m);
    let u_big_int = vec_string_to_vec_bigint(&input.u);
    let e0_big_int = vec_string_to_vec_bigint(&input.e0);
    let e1_big_int = vec_string_to_vec_bigint(&input.e1);
    let c0_big_int = vec_string_to_vec_bigint(&input.c0);
    let c1_big_int = vec_string_to_vec_bigint(&input.c1);
    let cyclo_big_int = vec_string_to_vec_bigint(&input.cyclo);

    // The circuit logic requires to access some random value
    // In order to draw randomness within the circuit we use Axiom's Challenge API (https://hackmd.io/@axiom/SJw3p-qX3)
    // Challenge API requires a Phase 0 of witness generation. During this phase, all the input polynomials are assigned to the circuit witness table.
    // A commitment from the witness generated during Phase 0 is extracted and then hashed to generate the random value according to Fiat-Shamir heuristic.
    // This random challenge can be then used as part of witness generation during Phase 1. We will need this to perform efficient polynomial multiplication.
    // Note that if you wanna verify something with the challenge API (eg enforcing polynomial multiplcation), the stuffs you verify (namely the input polynomials)
    // must be assigned in phase 0 so their values can be part of the Phase 0 commtiment and contribute to the challenge (gamma) used in Phase 1.

    // Phase 0: Assign the input polynomials to the circuit witness table
    let pk0 = poly_big_int_assign::<{ DEG - 1 }, F>(ctx, &pk0_big_int);
    let pk1 = poly_big_int_assign::<{ DEG - 1 }, F>(ctx, &pk1_big_int);
    let m = poly_big_int_assign::<{ DEG - 1 }, F>(ctx, &m_big_int);
    let u = poly_big_int_assign::<{ DEG - 1 }, F>(ctx, &u_big_int);
    let e0 = poly_big_int_assign::<{ DEG - 1 }, F>(ctx, &e0_big_int);
    let e1 = poly_big_int_assign::<{ DEG - 1 }, F>(ctx, &e1_big_int);
    let c0 = poly_big_int_assign::<{ DEG - 1 }, F>(ctx, &c0_big_int);
    let c1 = poly_big_int_assign::<{ DEG - 1 }, F>(ctx, &c1_big_int);
    let cyclo = poly_big_int_assign::<{ DEG }, F>(ctx, &cyclo_big_int);

    // DELTA is equal to Q/T rounded to the lower integer from BFV paper
    const DELTA: u64 = Q / T;

    // assign constant DELTA to the circuit
    let delta = ctx.load_witness(F::from(DELTA));

    // Expose to the public pk0 and pk1, c0, c1, cyclo and delta
    for &assigned_coefficient_pk0 in pk0.iter().take(DEG) {
        make_public.push(assigned_coefficient_pk0);
    }

    for &assigned_coefficient_pk1 in pk1.iter().take(DEG) {
        make_public.push(assigned_coefficient_pk1);
    }

    for &assigned_coefficient_c0 in c0.iter().take(DEG) {
        make_public.push(assigned_coefficient_c0);
    }

    for &assigned_coefficient_c1 in c1.iter().take(DEG) {
        make_public.push(assigned_coefficient_c1);
    }

    for &assigned_coefficient_cyclo in cyclo.iter().take(DEG + 1) {
        make_public.push(assigned_coefficient_cyclo);
    }

    make_public.push(delta);

    // Assign the length of the input polynomials pk0, pk1 and u to the circuit.
    let poly_len = ctx.load_witness(F::from(DEG as u64));
    let pk0_with_length = PolyWithLength::new(pk0, poly_len);
    let pk1_with_length = PolyWithLength::new(pk1, poly_len);
    let u_with_length = PolyWithLength::new(u.clone(), poly_len);

    // Assign the length of the cyclotomic polynomial to the circuit. This is equal to DEG + 1
    let cyclo_len = ctx.load_witness(F::from((DEG + 1) as u64));
    let cyclo_with_length = PolyWithLength::new(cyclo, cyclo_len);

    // PRECOMPUTATION

    // In this section we perform some precomputations outside the circuit
    // The resulting polynomials are then assigned to the circuit witness table
    // We are gonna need these polynomials to enforce some constraints inside the circuit in phase 1

    // Compute the polynomial pk0 * u outside the circuit
    let pk0_u_unassigned = poly_mul(&pk0_big_int, &u_big_int);

    // Compute the polynomial pk1 * u outside the circuit
    let pk1_u_unassigned = poly_mul(&pk1_big_int, &u_big_int);

    // Assert that all the coefficients of pk0_u_unassigned and pk1_u_unassigned are in the [0, p - 1] range, where p is the modulus of the circuit field.
    // If that is not the case, the value assigned to the circuit witness will be equal to the value modulo p, which is not what we want. This will eventually fail the `constrain_poly_mul` constraint.
    let p = BigInt::from_str_radix(&F::MODULUS[2..], 16).unwrap();
    for coeff in pk0_u_unassigned.iter() {
        assert!(coeff < &p);
    }

    for coeff in pk1_u_unassigned.iter() {
        assert!(coeff < &p);
    }

    // Reduce pk0_u_unassigned by modulo Q
    let pk0_u_reduced = reduce_poly_by_modulo_q::<Q>(&pk0_u_unassigned);

    // Reduce pk1_u_unassigned by modulo Q
    let pk1_u_reduced = reduce_poly_by_modulo_q::<Q>(&pk1_u_unassigned);

    // Compute the division between pk0_u_reduced and cyclo outside the circuit
    let (quotient_0, mut remainder_0) =
        div_euclid::<{ 2 * DEG - 2 }, DEG>(&pk0_u_reduced, &cyclo_big_int);

    // Compute the division between pk1_u_reduced and cyclo outside the circuit
    let (quotient_1, mut remainder_1) =
        div_euclid::<{ 2 * DEG - 2 }, DEG>(&pk1_u_reduced, &cyclo_big_int);

    // Note: that this is not a constraint enforced inside the circuit
    assert_eq!(quotient_0.len() - 1, DEG - 2);
    assert_eq!(quotient_1.len() - 1, DEG - 2);
    assert_eq!(remainder_0.len() - 1, DEG - 1);
    assert_eq!(remainder_1.len() - 1, DEG - 1);

    // Pad remainder_0 and remainder_1 with zeroes at the beginning to make their degree equal to 2 * DEG - 2 (required by the `constrain_poly_reduction_by_cyclo` chip used later in phase 1)
    while remainder_0.len() - 1 < 2 * DEG - 2 {
        remainder_0.insert(0, BigInt::from(0u32));
    }

    while remainder_1.len() - 1 < 2 * DEG - 2 {
        remainder_1.insert(0, BigInt::from(0u32));
    }

    assert_eq!(remainder_0.len() - 1, 2 * DEG - 2);
    assert_eq!(remainder_1.len() - 1, 2 * DEG - 2);

    // quotient_0 and quotient_1 have coefficients in the [0, Q-1] range
    // remainder_0 and remainder_1 have coefficients in the [-(Q-1), Q-1] range
    // Why is that? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Division section

    // Reduce the coefficients of remainder_0 and remainder_1 by modulo Q to make them in the [0, Q-1] range
    let remainder_0 = reduce_poly_by_modulo_q::<Q>(&remainder_0);
    let remainder_1 = reduce_poly_by_modulo_q::<Q>(&remainder_1);

    // Compute the polynomial multiplication between quotient_0 * cyclo outside the circuit
    let quotient_0_times_cyclo = poly_mul(&quotient_0, &cyclo_big_int);

    // Compute the polynomial multiplication between quotient_1 * cyclo outside the circuit
    let quotient_1_times_cyclo = poly_mul(&quotient_1, &cyclo_big_int);

    // Assert that all the coefficients of quotient_0_times_cyclo and quotient_1_times_cyclo are in the [0, p - 1] range.
    // If that is not the case, the value assigned to the circuit witness will be equal to the value modulo p, which is not what we want. This will eventually fail the `constrain_poly_mul` constraint.
    let p = BigInt::from_str_radix(&F::MODULUS[2..], 16).unwrap();
    for coeff in quotient_0_times_cyclo.iter() {
        assert!(coeff < &p);
    }

    for coeff in quotient_1_times_cyclo.iter() {
        assert!(coeff < &p);
    }

    // Precomputation is over, now we can assign the resulting polynomials to the circuit witness table
    // Note: we are still in Phase 0 of the witness generation

    // assign pk0_u_unassigned and pk1_u_unassigned to the circuit
    let pk0_u = poly_big_int_assign::<{ 2 * DEG - 2 }, F>(ctx, &pk0_u_unassigned);
    let pk1_u = poly_big_int_assign::<{ 2 * DEG - 2 }, F>(ctx, &pk1_u_unassigned);

    // Note: this is not a constraint enforced inside the circuit. The constraint on the degree of pk0_u and pk1_u is enforced inside `poly_big_int_assign`.
    assert!(pk0_u.len() - 1 == 2 * DEG - 2);
    assert!(pk1_u.len() - 1 == 2 * DEG - 2);

    // Assign the length of the polynomial pk0_u (and pk1_u) to the circuit -> this is equal to 2 * DEG - 1
    let len = ctx.load_witness(F::from((2 * DEG - 1) as u64));

    let pk0_u_with_length = PolyWithLength::new(pk0_u.clone(), len);
    let pk1_u_with_length = PolyWithLength::new(pk1_u.clone(), len);

    // Assign quotient_0 and quotient_1 to the circuit
    let quotient_0 = poly_big_int_assign::<{ DEG - 2 }, F>(ctx, &quotient_0);
    let quotient_1 = poly_big_int_assign::<{ DEG - 2 }, F>(ctx, &quotient_1);

    // assert that the degree of quotient_0 and quotient_1 is DEG - 2. The constraint on the degree of quotient_0 and quotient_1 is enforced inside `poly_big_int_assign`
    assert!(quotient_0.len() - 1 == DEG - 2);
    assert!(quotient_1.len() - 1 == DEG - 2);

    // assign quotient_0 (and quotient_1) length to the circuit -> this is equal to DEG - 1
    let quotient_len = ctx.load_witness(F::from((DEG - 1) as u64));

    let quotient_0_with_length = PolyWithLength::new(quotient_0, quotient_len);
    let quotient_1_with_length = PolyWithLength::new(quotient_1, quotient_len);

    // Assign quotient_0_times_cyclo and quotient_1_times_cyclo to the circuit
    let quotient_0_times_cyclo =
        poly_big_int_assign::<{ 2 * DEG - 2 }, F>(ctx, &quotient_0_times_cyclo);

    let quotient_1_times_cyclo =
        poly_big_int_assign::<{ 2 * DEG - 2 }, F>(ctx, &quotient_1_times_cyclo);

    // assert that the degree of quotient_0_times_cyclo and quotient_1_times_cyclo is 2 * DEG - 2. The constraint on the degree of quotient_0_times_cyclo and quotient_1_times_cyclo is enforced inside `poly_big_int_assign`
    assert!(quotient_0_times_cyclo.len() - 1 == 2 * DEG - 2);
    assert!(quotient_1_times_cyclo.len() - 1 == 2 * DEG - 2);

    // Assign the length of the polynomial quotient_0_times_cyclo (and quotient_1_times_cyclo) to the circuit -> this is equal to 2 * DEG - 1
    let len = ctx.load_witness(F::from((2 * DEG - 1) as u64));

    let quotient_0_times_cyclo_with_length = PolyWithLength::new(quotient_0_times_cyclo, len);
    let quotient_1_times_cyclo_with_length = PolyWithLength::new(quotient_1_times_cyclo, len);

    // Assign the remainder_0 and remainder_1 to the circuit
    let remainder_0 = poly_big_int_assign::<{ 2 * DEG - 2 }, F>(ctx, &remainder_0);
    let remainder_1 = poly_big_int_assign::<{ 2 * DEG - 2 }, F>(ctx, &remainder_1);

    // assert that the degree of quotient_0_times_cyclo and quotient_1_times_cyclo is 2 * DEG - 1. The constraint on the degree of quotient_0_times_cyclo and quotient_1_times_cyclo is enforced inside `poly_big_int_assign`
    assert!(remainder_0.len() - 1 == 2 * DEG - 2);
    assert!(remainder_1.len() - 1 == 2 * DEG - 2);

    let len = ctx.load_witness(F::from((2 * DEG - 1) as u64));

    let remainder_0_with_length = PolyWithLength::new(remainder_0, len);
    let remainder_1_with_length = PolyWithLength::new(remainder_1, len);

    // Phase 0 is over, we can now move to Phase 1, in which we will leverage the random challenge generated during Phase 0.
    // According to the design of this API, all the constraints must be written inside a callback function.

    #[allow(clippy::let_and_return)]
    let callback = move |ctx_gate: &mut Context<F>,
                         ctx_rlc: &mut Context<F>,
                         eth_chip: &EthChip<F>| {
        let range = eth_chip.range();
        let rlc = eth_chip.rlc();

        /* constraint on e0
            - e0 must be a polynomial in the R_q ring => Coefficients must be in the [0, Q-1] range and the degree of e0 must be DEG - 1
            - e0 must be sampled from the distribution ChiError, namely the coefficients of e0 must be in the [0, b] OR [q-b, q-1] range

            Approach:
            - `check_poly_from_distribution_chi_error` chip guarantees that the coefficients of e0 are in the range [0, b] OR [q-b, q-1]
            - As this range is a subset of the [0, Q-1] range, the coefficients of e0 are guaranteed to be in the [0, Q-1] range
            - The assignment for loop above enforces that the degree of e0 is DEG - 1
        */

        /* constraint on e1
            Same as e0
        */

        check_poly_coefficients_in_range::<{ DEG - 1 }, Q, B, F>(ctx_gate, &e0, range);
        check_poly_coefficients_in_range::<{ DEG - 1 }, Q, B, F>(ctx_gate, &e1, range);

        /* constraint on u
            - u must be a polynomial in the R_q ring => Coefficients must be in the [0, Q-1] range and the degree of u must be DEG - 1
            - u must be sampled from the distribution ChiKey, namely the coefficients of u must be either 0, 1 or Q-1

            Approach:
            - `check_poly_from_distribution_chi_key` chip guarantees that the coefficients of u are either 0, 1 or Q-1
            - As this range is a subset of the [0, Q-1] range, the coefficients of u are guaranteed to be in the [0, Q-1] range
            - The assignment for loop above guarantees that the degree of u is DEG - 1
        */

        check_poly_from_distribution_chi_key::<{ DEG - 1 }, Q, F>(ctx_gate, &u, range.gate());

        /* constraint on m
            - m must be a polynomial in the R_t ring => Coefficients must be in the [0, T/2] OR [Q - T/2, Q - 1] range and the degree of m must be DEG - 1

            Approach:
            - Perform a range check on the coefficients of m to be in the [0, T/2] OR [Q - T/2, Q - 1] range
            - The assignment for loop above guarantees that the degree of m is DEG - 1
        */

        check_poly_coefficients_in_range::<{ DEG - 1 }, Q, { T / 2 }, F>(ctx_gate, &m, range);

        // 1. COMPUTE C0 (c0 is the first ciphertext component)

        // 1.1 pk0 * u
        // Constrain the polynomial multiplication between pk0 and u to be equal to pk0_u using the `constrain_poly_mul` chip

        constrain_poly_mul(
            pk0_with_length,
            u_with_length.clone(),
            pk0_u_with_length.clone(),
            ctx_gate,
            ctx_rlc,
            rlc,
            range.gate(),
        );

        // pk0_u is a polynomial of degree (DEG - 1) * 2 = 2*DEG - 2
        // pk0_u has coefficients in the [0, (Q-1) * (Q-1) * (DEG+1)] range. This is the maximum value that a coefficient of pk0_u can take. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Multiplication section

        // 1.2 Reduce the coefficients of pk0_u by modulo `Q`

        // get the number of bits needed to represent the value of (Q-1) * (Q-1) * (DEG+1)
        let q_minus_1 = BigInt::from(Q) - BigInt::from(1u32);
        let deg_plus_one = BigInt::from(DEG as u64) + BigInt::from(1u32);
        let binary_representation = format!("{:b}", (q_minus_1.clone() * q_minus_1 * deg_plus_one));
        let num_bits_1 = binary_representation.len();

        // The coefficients of pk0_u are in the range [0, (Q-1) * (Q-1) * (DEG+1)] according to the above analysis.
        // Therefore the coefficients of pk0_u are known to have <= `num_bits_1` bits, therefore satisfying the assumption of the `poly_reduce_by_modulo_q` chip

        let pk0_u =
            poly_reduce_by_modulo_q::<{ 2 * DEG - 2 }, Q, F>(ctx_gate, &pk0_u, range, num_bits_1);

        // pk0_u is a polynomial of degree (DEG - 1) * 2 = 2 * DEG - 2
        // pk0_u now has coefficients in the [0, Q-1] range given the constraint above
        // cyclo is a polynomial of degree DEG

        // 1.3 Reduce pk0_u by the cyclo polynomial

        for i in 0..DEG - 1 {
            range.check_less_than_safe(ctx_gate, quotient_0_with_length.get_poly()[i], Q);
        }

        for i in 0..2 * DEG - 1 {
            range.check_less_than_safe(ctx_gate, remainder_0_with_length.get_poly()[i], Q);
        }

        // The assumption of the `constrain_poly_reduction_by_cyclo` chip are that:
        // * Assumption: the coefficients of quotient are constrained to be in the range [0, Q - 1] -> met according to the constraint just set above
        // * Assumption: the coefficients of remainder are constrained to be in the range [0, Q - 1] -> met according to the constraint just set above
        // * Assumption: the coefficients are constrained such to avoid overflow during the polynomial addition between `quotient_times_cyclo` and `remainder`
        //      -> remainder has coefficients in the [0, Q - 1] range according to the constraint just set above
        //      -> quotient_times_cyclo has coefficients in the [0, (Q-1) * (Q-1) * (DEG+1)] range according to the above analysis
        //      -> (Q-1) * (Q-1) * (DEG+1) + (Q-1) is the maximum value that a quotient of remainder + quotient_times_cyclo can take. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Addition section
        //         Therefore Q and DEG must be chosen such that (Q-1) * (Q-1) * (DEG+1) + (Q-1) < p, where p is the modulus of the circuit field.

        constrain_poly_reduction_by_cyclo::<{ 2 * DEG - 2 }, DEG, Q, F>(
            &pk0_u,
            cyclo_with_length.clone(),
            quotient_0_with_length,
            quotient_0_times_cyclo_with_length,
            remainder_0_with_length.clone(),
            range,
            ctx_gate,
            ctx_rlc,
            rlc,
        );

        let pk0_u = remainder_0_with_length.get_poly().clone();

        assert_eq!(pk0_u.len() - 1, 2 * DEG - 2);

        // pk0_u is a polynomial of degree 2*DEG - 2
        // pk0_u now has coefficients in the [0, Q-1] range

        // But actually, the degree of pk0_u should be DEG - 1 after reduction by the cyclo polynomial, the first DEG - 1 coefficients are just zeroes

        // 1.4 Enforce that the first DEG - 1 coefficients of pk0_u are zeroes
        for pk0_u_element in pk0_u.iter().take(DEG - 1) {
            let bool = range
                .gate()
                .is_equal(ctx_gate, *pk0_u_element, Constant(F::from(0)));
            range.gate().assert_is_const(ctx_gate, &bool, &F::from(1));
        }

        // Therefore, we can safely trim the first DEG - 1 coefficients from pk0_u

        let pk0_u_trimmed: Vec<_> = pk0_u.iter().skip(DEG - 1).cloned().collect();

        // assert that the degree of pk0_u_trimmed is DEG - 1
        assert_eq!(pk0_u_trimmed.len() - 1, DEG - 1);

        // pk0_u_trimmed is a polynomial in the R_q ring!

        // 1.5 m * delta

        // Perform the polynomial scalar multiplication between m and delta.
        // The assumption of the `poly_scalar_mul` chip is that the coefficients of the input polynomials are constrained such to avoid overflow during the polynomial scalar multiplication
        // m has coefficients in the [0, T/2] OR [Q - T/2, Q - 1] range according to the constraint just set above
        // delta is a constant equal to Q/T rounded to the lower integer from BFV paper
        // The coefficients of m_delta are in the [0, (Q-1) * (Q/T)] range. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Scalar Multiplication section
        // Q and T must be chosen such that (Q-1) * (Q/T) < p, where p is the modulus of the circuit field.

        let m_delta = poly_scalar_mul::<{ DEG - 1 }, F>(ctx_gate, &m, &delta, range.gate());

        // m_delta is a polynomial of degree DEG - 1

        // 1.6 Reduce the coefficients of `m_delta` by modulo `Q`

        // get the number of bits needed to represent the value of (Q-1) * (Q/T)

        let binary_representation = format!("{:b}", ((Q - 1) * (Q / T)));
        let num_bits_2 = binary_representation.len();

        // The coefficients of m_delta are in the range [0,  (Q-1) * (Q/T)] according to the polynomial scalar multiplication constraint set above.
        // Therefore the coefficients of m_delta are known to have <= `num_bits_2` bits, therefore they satisfy the assumption of the `poly_reduce_by_modulo_q` chip

        let m_delta =
            poly_reduce_by_modulo_q::<{ DEG - 1 }, Q, F>(ctx_gate, &m_delta, range, num_bits_2);

        // Note: Scalar multiplication does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1
        // m_delta is a polynomial in the R_q ring

        // 1.7 pk0_u_trimmed + m_delta

        // Perform the polynomial addition between pk0_u_trimmed and m_delta.
        // The assumption of the `poly_add` chip is that the coefficients of the input polynomials are constrained such to avoid overflow during the polynomial addition
        // `m_delta` has coefficients in the [0, Q-1] range according to the constraint just set above
        // `pk0_u_trimmed` has coefficients in the [0, Q-1] range according to the constraint just set above
        // The coefficients of pk0_u_trimmed_plus_m_delta are in the [0, 2Q - 2] range. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Addition section
        // Q must be chosen such that 2Q - 2 < p, where p is the modulus of the circuit field.
        let pk0_u_trimmed_plus_m_delta =
            poly_add::<{ DEG - 1 }, F>(ctx_gate, &pk0_u_trimmed, &m_delta, range.gate());

        // 1.8 Reduce the coefficients of `pk0_u_trimmed_plus_m_delta` by modulo `Q`

        // get the number of bits needed to represent the value of 2Q - 2

        let binary_representation = format!("{:b}", (2 * Q - 2));
        let num_bits_3 = binary_representation.len();

        // The coefficients of pk0_u_trimmed_plus_m_delta are in the range [0, 2Q - 2] according to the polynomial addition constraint set above.
        // Therefore the coefficients of m_delta are known to have <= `num_bits_3` bits, therefore they satisfy the assumption of the `poly_reduce_by_modulo_q` chip

        let pk0_u_trimmed_plus_m_delta = poly_reduce_by_modulo_q::<{ DEG - 1 }, Q, F>(
            ctx_gate,
            &pk0_u_trimmed_plus_m_delta,
            range,
            num_bits_3,
        );

        // Note: Addition does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1
        // pk0_u_trimmed_plus_m_delta is a polynomial in the R_q ring

        // 1.9 c0 = pk0_u_trimmed_plus_m_delta + e0

        // Perform the polynomial addition between pk0_u_trimmed_plus_m_delta and e0.
        // The assumption of the `poly_add` chip is that the coefficients of the input polynomials are constrained such to avoid overflow during the polynomial addition
        // `pk0_u_trimmed_plus_m_delta` has coefficients in the [0, Q-1] range according to the constraint just set above
        // `e0` has coefficients in the [0, b] OR [q-b, q-1] range according to the constraint just set above
        // The coefficients of computed_c0 are in the [0, 2Q - 2] range. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Addition section
        // Q must be chosen such that 2Q - 2 < p, where p is the modulus of the circuit field.
        let computed_c0 =
            poly_add::<{ DEG - 1 }, F>(ctx_gate, &pk0_u_trimmed_plus_m_delta, &e0, range.gate());

        // Coefficients of computed_c0 are in the [0, 2Q - 2] range according to the polynomial addition constraint set above.
        // Therefore the coefficients of computed_c0 are known to have <= `num_bits_3` bits, therefore they satisfy the assumption of the `poly_reduce_by_modulo_q` chip

        // 1.10 Reduce the coefficients of `pk0_u_trimmed_plus_m_delta` by modulo `Q`

        let computed_c0 =
            poly_reduce_by_modulo_q::<{ DEG - 1 }, Q, F>(ctx_gate, &computed_c0, range, num_bits_3);

        // Note: Addition does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1
        // computed_c0 is a polynomial in the R_q ring!

        // 2. COMPUTE C1 (c1 is the second ciphertext component)

        // 2.1 pk1 * u
        // Constrain the polynomial multiplication between pk1 and u to be equal to pk§_u using the `constrain_poly_mul` chip

        constrain_poly_mul(
            pk1_with_length,
            u_with_length,
            pk1_u_with_length,
            ctx_gate,
            ctx_rlc,
            rlc,
            range.gate(),
        );

        // pk1_u is a polynomial of degree (DEG - 1) * 2 = 2*DEG - 2
        // pk1_u has coefficients in the [0, (Q-1) * (Q-1) * (DEG+1)] range. This is the maximum value that a coefficient of pk0_u can take. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Multiplication section

        // 2.2 Reduce the coefficients by modulo `Q`

        // The coefficients of pk1_u are in the range [0, (Q-1) * (Q-1) * (DEG+1)] according to the above analysis.
        // Therefore the coefficients of pk1_u are known to have <= `num_bits_1` bits, therefore satisfying the assumption of the `poly_reduce_by_modulo_q` chip

        let pk1_u =
            poly_reduce_by_modulo_q::<{ 2 * DEG - 2 }, Q, F>(ctx_gate, &pk1_u, range, num_bits_1);

        // pk0_u is a polynomial of degree (DEG - 1) * 2 = 2 * DEG - 2
        // pk0_u now has coefficients in the [0, Q-1] range given the constraint above
        // cyclo is a polynomial of degree DEG

        // 2.3 Reduce pk1_u by the cyclo polynomial

        for i in 0..DEG - 1 {
            range.check_less_than_safe(ctx_gate, quotient_1_with_length.get_poly()[i], Q);
        }

        for i in 0..2 * DEG - 1 {
            range.check_less_than_safe(ctx_gate, remainder_1_with_length.get_poly()[i], Q);
        }

        // For the analysis of the assumptions of the `constrain_poly_reduction_by_cyclo` chip, see the comments above in the `1.3 Reduce pk0_u by the cyclo polynomial` section
        constrain_poly_reduction_by_cyclo::<{ 2 * DEG - 2 }, DEG, Q, F>(
            &pk1_u,
            cyclo_with_length.clone(),
            quotient_1_with_length,
            quotient_1_times_cyclo_with_length,
            remainder_1_with_length.clone(),
            range,
            ctx_gate,
            ctx_rlc,
            rlc,
        );

        let pk1_u = remainder_1_with_length.get_poly().clone();

        // assert that the degree of pk1_u is 2*DEG - 2

        // pk1_u is a polynomial of degree 2*DEG - 2
        // pk1_u now has coefficients in the [0, Q-1] range

        // But actually, the degree of pk1_u should be DEG - 1 after reduction by the cyclo polynomial, the first DEG - 1 coefficients are just zeroes

        // 2.4 Enforce that the first DEG - 1 coefficients of pk0_u are zeroes
        for pk1_u_element in pk1_u.iter().take(DEG - 1) {
            let bool = range
                .gate()
                .is_equal(ctx_gate, *pk1_u_element, Constant(F::from(0)));
            range.gate().assert_is_const(ctx_gate, &bool, &F::from(1));
        }

        // Therefore, we can safely trim the first DEG - 1 coefficients from pk1_u

        let pk1_u_trimmed: Vec<_> = pk1_u.iter().skip(DEG - 1).cloned().collect();

        // assert that the degree of pk1_u_trimmed is DEG - 1
        assert_eq!(pk1_u_trimmed.len() - 1, DEG - 1);

        // pk1_u_trimmed is a polynomial in the R_q ring!

        // 2.5 c1 = pk1_u_trimmed + e1

        // Perform the polynomial addition between pk1_u_trimmed and e1.
        // The assumption of the `poly_add` chip is that the coefficients of the input polynomials are constrained such to avoid overflow during the polynomial addition
        // `pk1_u_trimmed` has coefficients in the [0, Q-1] range according to the constraint just set above
        // `e1` has coefficients in the [0, b] OR [q-b, q-1] range according to the constraint just set above
        // The coefficients of computed_c1 are in the [0, 2Q - 2] range. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Addition section
        // Q must be chosen such that 2Q - 2 < p, where p is the modulus of the circuit field.
        let computed_c1 = poly_add::<{ DEG - 1 }, F>(ctx_gate, &pk1_u_trimmed, &e1, range.gate());

        // The coefficients of computed_c1 are in the range [0, 2Q - 2] according to the polynomial addition performed above.
        // Therefore the coefficients of computed_c1 are known to have <= `num_bits_3` bits, therefore they satisfy the assumption of the `poly_reduce_by_modulo_q` chip

        // 2.6 Reduce the coefficients of `computed_c1` by modulo `Q`

        let computed_c1 =
            poly_reduce_by_modulo_q::<{ DEG - 1 }, Q, F>(ctx_gate, &computed_c1, range, num_bits_3);

        // Note: Addition does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1
        // computed_c1 is a polynomial in the R_q ring

        // Enforce equality between `c0` and `computed_c0` using equality check
        // Enfroce equality between `c1` and `computed_c1` using equality check
        for i in 0..DEG {
            let bool_0 = range.gate().is_equal(ctx_gate, c0[i], computed_c0[i]);
            range.gate().assert_is_const(ctx_gate, &bool_0, &F::from(1));

            let bool_1 = range.gate().is_equal(ctx_gate, c1[i], computed_c1[i]);
            range.gate().assert_is_const(ctx_gate, &bool_1, &F::from(1));
        }
    };

    callback
}

fn main() {
    env_logger::init();

    let args = Cli::parse();

    run_eth(bfv_encryption_circuit, args);
}
