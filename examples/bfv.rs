use axiom_eth::keccak::KeccakChip;
use axiom_eth::{EthChip, Field};
use clap::Parser;
use halo2_base::safe_types::RangeInstructions;
use halo2_base::{AssignedValue, Context};
use halo2_scaffold::scaffold::{cmd::Cli, run_eth};
use serde::{Deserialize, Serialize};
use zk_fhe::chips::poly_assigned::PolyAssigned;
use zk_fhe::chips::poly_distribution::{
    check_poly_coefficients_in_range, check_poly_from_distribution_chi_key, reduce_by_modulo_q,
};
use zk_fhe::poly::Poly;

/// Circuit inputs for BFV encryption operations
///
/// # Type Parameters
///
/// * `N`: Degree of the cyclotomic polynomial `cyclo` of the polynomial ring.
/// * `Q`: Modulus of the ciphertext space
/// * `T`: Modulus of the plaintext space
/// * `B`: Upper bound of the Gaussian distribution Chi Error. It is defined as 6 * 𝜎
///
/// # Fields
///
/// * `pk0`: Public key 0 - polynomial of degree N-1 living in ciphertext space R_q
/// * `pk1`: Public key 1 - polynomial of degree N-1 living in ciphertext space R_q
/// * `m`: Plaintext message to be encrypted - polynomial of degree N-1 living in plaintext space R_t
/// * `u`: Ephemeral key - polynomial of degree N-1 living in ciphertext space R_q - its coefficients are sampled from the distribution ChiKey
/// * `e0`: Error - polynomial of degree N-1 living in ciphertext space R_q - its coefficients are sampled from the distribution ChiError
/// * `e1`: Error - polynomial of degree N-1 living in ciphertext space R_q - its coefficients are sampled from the distribution ChiError
/// * `c0`: First ciphertext component - polynomial of degree N-1 living in ciphertext space R_q
/// * `c1`: Second ciphertext component - polynomial of degree N-1 living in ciphertext space R_q
/// * `cyclo`: Cyclotomic polynomial of degree N in the form x^N + 1
///
/// Note: all the polynomials are expressed by their coefficients in the form [a_N-1, a_N-2, ..., a_1, a_0] where a_0 is the constant term
///
/// # Assumptions (to be checked on the public inputs outside the circuit)
///
/// * `N` must be a power of 2
/// * `Q` must be a prime number and be greater than 1.
/// * `T` must be a prime number and must be greater than 1 and less than `Q`
/// * `B` must be a positive integer and must be less than `Q`
/// * `cyclo` must be the cyclotomic polynomial of degree `N` in the form x^N + 1
/// * `pk0` and `pk1` must be polynomials in the R_q ring. The ring R_q is defined as R_q = Z_q[x]/(x^N + 1)
/// *  Q and N must be chosen such that (Q-1) * (Q-1) * (N+1) + (Q-1) < p, where p is the modulus of the circuit field to avoid overflow during polynomial addition inside the circuit
/// *  Q, T must be chosen such that (Q-1) * (Q-T) + (Q-1) + (Q-1) < p, where p is the modulus of the circuit field.. This is required to avoid overflow during polynomial scalar multiplication inside the circuit
/// *  Q must be chosen such that 2Q - 2 < p, where p is the modulus of the circuit field. This is required to avoid overflow during polynomial addition inside the circuit

// For real world applications, the parameters should be chosen according to the security level required.
// N and Q Parameters of the BFV encryption scheme should be chosen according to TABLES of RECOMMENDED PARAMETERS for 128-bits security level
// https://homomorphicencryption.org/wp-content/uploads/2018/11/HomomorphicEncryptionStandardv1.1.pdf
// B is the upper bound of the distribution Chi Error. Pick standard deviation 𝜎 ≈ 3.2 according to the HomomorphicEncryptionStandardv1 paper.
// T is picked according to Lattigo (https://github.com/tuneinsight/lattigo/blob/master/schemes/bfv/example_parameters.go) implementation
// As suggest by https://eprint.iacr.org/2021/204.pdf (paragraph 2) B = 6σerr
// These are just parameters used for fast testing purpose - to match with input file `data/bfv.in`
const N: usize = 4;
const Q: u64 = 4637;
const T: u64 = 7;
const B: u64 = 19;

// These are the parameters used for the real world application - to match with input file `data/bfv_2.in`
// const N: usize = 1024;
// const Q: u64 = 536870909;
// const T: u64 = 7;
// const B: u64 = 19;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CircuitInput<const N: usize, const Q: u64, const T: u64, const B: u64> {
    pub pk0: Vec<String>, // PUBLIC INPUT. Should live in R_q according to assumption
    pub pk1: Vec<String>, // PUBLIC INPUT. Should live in R_q according to assumption
    pub m: Vec<String>,   // PRIVATE INPUT. Should in R_t (enforced inside the circuit)
    pub u: Vec<String>, // PRIVATE INPUT. Should live in R_q and be sampled from the distribution ChiKey (enforced inside the circuit)
    pub e0: Vec<String>, // PRIVATE INPUT. Should live in R_q and be sampled from the distribution ChiError (enforced inside the circuit)
    pub e1: Vec<String>, // PRIVATE INPUT. Should live in R_q and be sampled from the distribution ChiError (enforced inside the circuit)
    pub c0: Vec<String>, // PUBLIC INPUT. Should live in R_q. We constraint equality between c0 and computed_c0 namely the ciphertext computed inside the circuit
    pub c1: Vec<String>, // PUBLIC INPUT. Should live in R_q. We constraint equality between c1 and computed_c1 namely the ciphertext computed inside the circuit
    pub cyclo: Vec<String>, // PUBLIC INPUT. Should be the cyclotomic polynomial of degree N in the form x^N + 1 according to assumption
}

fn bfv_encryption_circuit<F: Field>(
    ctx: &mut Context<F>,
    _eth_chip: &EthChip<F>,
    _keccak: &mut KeccakChip<F>,
    input: CircuitInput<N, Q, T, B>,
    make_public: &mut Vec<AssignedValue<F>>,
) -> impl FnOnce(&mut Context<F>, &mut Context<F>, &EthChip<F>) + Clone {
    // Transform the input polynomials strings into Polys
    let pk0_unassigned = Poly::from_string(input.pk0);
    let pk1_unassigned = Poly::from_string(input.pk1);
    let m_unassigned = Poly::from_string(input.m);
    let u_unassigned = Poly::from_string(input.u);
    let e0_unassigned = Poly::from_string(input.e0);
    let e1_unassigned = Poly::from_string(input.e1);
    let c0_unassigned = Poly::from_string(input.c0);
    let c1_unassigned = Poly::from_string(input.c1);
    let cyclo_unassigned = Poly::from_string(input.cyclo);

    // The circuit logic requires to access some random value
    // In order to draw randomness within the circuit we use Axiom's Challenge API (https://hackmd.io/@axiom/SJw3p-qX3)
    // Challenge API requires a Phase 0 of witness generation. During this phase, all the input polynomials are assigned to the circuit witness table.
    // A commitment from the witness generated during Phase 0 is extracted and then hashed to generate the random value according to Fiat-Shamir heuristic.
    // This random challenge can be then used as part of witness generation during Phase 1. We will need this to perform efficient polynomial multiplication.
    // Note that if you wanna verify something with the challenge API (eg enforcing polynomial multiplcation), the stuffs you verify (namely the input polynomials)
    // must be assigned in phase 0 so their values can be part of the Phase 0 commtiment and contribute to the challenge (gamma) used in Phase 1.

    // // Phase 0: Assign the input polynomials to the circuit witness table
    let pk0 = PolyAssigned::<{ N - 1 }, F>::new(pk0_unassigned.clone(), ctx);
    let pk1 = PolyAssigned::<{ N - 1 }, F>::new(pk1_unassigned.clone(), ctx);
    let m = PolyAssigned::<{ N - 1 }, F>::new(m_unassigned, ctx);
    let u = PolyAssigned::<{ N - 1 }, F>::new(u_unassigned.clone(), ctx);
    let e0 = PolyAssigned::<{ N - 1 }, F>::new(e0_unassigned, ctx);
    let e1 = PolyAssigned::<{ N - 1 }, F>::new(e1_unassigned, ctx);
    let c0 = PolyAssigned::<{ N - 1 }, F>::new(c0_unassigned, ctx);
    let c1 = PolyAssigned::<{ N - 1 }, F>::new(c1_unassigned, ctx);
    let cyclo = PolyAssigned::<{ N }, F>::new(cyclo_unassigned.clone(), ctx);

    // DELTA is equal to Q/T rounded to the lower integer from BFV paper
    const DELTA: u64 = Q / T;

    // assign constant DELTA to the circuit
    let delta = ctx.load_constant(F::from(DELTA));

    // Make pk0, pk1, c0, c1, cyclo public
    pk0.to_public(make_public);
    pk1.to_public(make_public);
    c0.to_public(make_public);
    c1.to_public(make_public);
    cyclo.to_public(make_public);

    // PRECOMPUTATION

    // In this section we perform some precomputations outside the circuit
    // The resulting polynomials are then assigned to the circuit witness table
    // We are gonna need these polynomials to enforce some constraints inside the circuit in phase 1

    // Compute pk0 * u outside the circuit
    let mut pk0_u_unassigned = pk0_unassigned.mul(&u_unassigned);

    // Compute pk1 * u outside the circuit
    let mut pk1_u_unassigned = pk1_unassigned.mul(&u_unassigned);

    // assign pk0_u_unassigned and pk1_u_unassigned to the circuit
    let pk0_u = PolyAssigned::<{ 2 * N - 2 }, F>::new(pk0_u_unassigned.clone(), ctx);
    let pk1_u = PolyAssigned::<{ 2 * N - 2 }, F>::new(pk1_u_unassigned.clone(), ctx);

    // Reduce pk0_u_unassigned and pk1_u_unassigned by modulo Q
    pk0_u_unassigned.reduce_by_modulus(Q);
    pk1_u_unassigned.reduce_by_modulus(Q);

    // Compute the division between pk0_u_unassigned and cyclo outside the circuit
    let (quotient_0_unassigned, mut remainder_0_unassigned) =
        pk0_u_unassigned.div_euclid(&cyclo_unassigned);

    // Compute the division between pk1_u_unassigned and cyclo outside the circuit
    let (quotient_1_unassigned, mut remainder_1_unassigned) =
        pk1_u_unassigned.div_euclid(&cyclo_unassigned);

    // quotient_0 and quotient_1 have coefficients in the [0, Q-1] range
    // remainder_0 and remainder_1 have coefficients in the [-(Q-1), Q-1] range
    // Why is that? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Division section

    // Reduce the coefficients of remainder_0 and remainder_1 by modulo Q to make them in the [0, Q-1] range
    remainder_0_unassigned.reduce_by_modulus(Q);
    remainder_1_unassigned.reduce_by_modulus(Q);

    // Compute the polynomial multiplication between quotient_0_unassigned * cyclo outside the circuit
    let quotient_0_times_cyclo_unassigned = quotient_0_unassigned.mul(&cyclo_unassigned);

    // Compute the polynomial multiplication between quotient_1_unassigned * cyclo outside the circuit
    let quotient_1_times_cyclo_unassigned = quotient_1_unassigned.mul(&cyclo_unassigned);

    // Precomputation is over, now we can assign the resulting polynomials to the circuit witness table
    // Note: we are still in Phase 0 of the witness generation

    // Assign quotient_0 and quotient_1 to the circuit
    let quotient_0 = PolyAssigned::<{ N - 2 }, F>::new(quotient_0_unassigned, ctx);
    let quotient_1 = PolyAssigned::<{ N - 2 }, F>::new(quotient_1_unassigned, ctx);

    // // Assign quotient_0_times_cyclo_unassigned and quotient_1_times_cyclo_unassigned to the circuit
    let quotient_0_times_cyclo =
        PolyAssigned::<{ 2 * N - 2 }, F>::new(quotient_0_times_cyclo_unassigned, ctx);
    let quotient_1_times_cyclo =
        PolyAssigned::<{ 2 * N - 2 }, F>::new(quotient_1_times_cyclo_unassigned, ctx);

    // Assign the remainder_0 and remainder_1 to the circuit
    let remainder_0 = PolyAssigned::<{ 2 * N - 2 }, F>::new(remainder_0_unassigned, ctx);
    let remainder_1 = PolyAssigned::<{ 2 * N - 2 }, F>::new(remainder_1_unassigned, ctx);

    // // Phase 0 is over, we can now move to Phase 1, in which we will leverage the random challenge generated during Phase 0.
    // // According to the design of this API, all the constraints must be written inside a callback function.

    #[allow(clippy::let_and_return)]
    let callback =
        move |ctx_gate: &mut Context<F>, ctx_rlc: &mut Context<F>, eth_chip: &EthChip<F>| {
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
            check_poly_coefficients_in_range::<{ N - 1 }, Q, B, F>(ctx_gate, &e0, range);
            check_poly_coefficients_in_range::<{ N - 1 }, Q, B, F>(ctx_gate, &e1, range);

            /* constraint on u
                - u must be a polynomial in the R_q ring => Coefficients must be in the [0, Q-1] range and the degree of u must be DEG - 1
                - u must be sampled from the distribution ChiKey, namely the coefficients of u must be either 0, 1 or Q-1

                Approach:
                - `check_poly_from_distribution_chi_key` chip guarantees that the coefficients of u are either 0, 1 or Q-1
                - As this range is a subset of the [0, Q-1] range, the coefficients of u are guaranteed to be in the [0, Q-1] range
                - The assignment for loop above guarantees that the degree of u is DEG - 1
            */
            check_poly_from_distribution_chi_key::<{ N - 1 }, Q, F>(ctx_gate, &u, range.gate());

            /* constraint on m
                - m must be a polynomial in the R_t ring => Coefficients must be in the [0, T/2] OR [Q - T/2, Q - 1] range and the degree of m must be DEG - 1

                Approach:
                - Perform a range check on the coefficients of m to be in the [0, T/2] OR [Q - T/2, Q - 1] range
                - The assignment for loop above guarantees that the degree of m is DEG - 1
            */
            check_poly_coefficients_in_range::<{ N - 1 }, Q, { T / 2 }, F>(ctx_gate, &m, range);

            // // 1. COMPUTE C0 (c0 is the first ciphertext component)

            // 1.1 pk0 * u
            // Constrain the polynomial multiplication between pk0 and u to be equal to pk0_u using the `constrain_poly_mul` chip

            pk0.constrain_poly_mul(u, pk0_u.clone(), ctx_gate, ctx_rlc, rlc);

            // pk0_u is a polynomial of degree (DEG - 1) * 2 = 2*DEG - 2
            // pk0_u has coefficients in the [0, (Q-1) * (Q-1) * (DEG+1)] range. This is the maximum value that a coefficient of pk0_u can take. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Multiplication section

            // 1.2 Reduce the coefficients of pk0_u by modulo `Q`
            let pk0_u = reduce_by_modulo_q::<{ 2 * N - 2 }, Q, F>(ctx_gate, range, &pk0_u);

            // pk0_u is a polynomial of degree (DEG - 1) * 2 = 2 * DEG - 2
            // pk0_u now has coefficients in the [0, Q-1] after reduction by modulo Q
            // cyclo is a polynomial of degree DEG

            // // 1.3 Reduce pk0_u by the cyclo polynomial

            // for i in 0..DEG - 1 {
            //     range.check_less_than_safe(ctx_gate, quotient_0_with_length.get_poly()[i], Q);
            // }

            // for i in 0..2 * DEG - 1 {
            //     range.check_less_than_safe(ctx_gate, remainder_0_with_length.get_poly()[i], Q);
            // }

            // // Assert that all the coefficients of quotient_0_times_cyclo and quotient_1_times_cyclo are in the [0, p - 1] range.
            // // If that is not the case, the value assigned to the circuit witness will be equal to the value modulo p, which is not what we want. This will eventually fail the `constrain_poly_mul` constraint.
            // let p = BigInt::from_str_radix(&F::MODULUS[2..], 16).unwrap();
            // for coeff in quotient_0_times_cyclo.iter() {
            //     assert!(coeff < &p);
            // }

            // for coeff in quotient_1_times_cyclo.iter() {
            //     assert!(coeff < &p);
            // }

            // // The assumption of the `constrain_poly_reduction_by_cyclo` chip are that:
            // // * Assumption: the coefficients of quotient are constrained to be in the range [0, Q - 1] -> met according to the constraint just set above
            // // * Assumption: the coefficients of remainder are constrained to be in the range [0, Q - 1] -> met according to the constraint just set above
            // // * Assumption: the coefficients are constrained such to avoid overflow during the polynomial addition between `quotient_times_cyclo` and `remainder`
            // //      -> remainder has coefficients in the [0, Q - 1] range according to the constraint just set above
            // //      -> quotient_times_cyclo has coefficients in the [0, (Q-1) * (Q-1) * (DEG+1)] range according to the above analysis
            // //      -> (Q-1) * (Q-1) * (DEG+1) + (Q-1) is the maximum value that a quotient of remainder + quotient_times_cyclo can take. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Addition section
            // //         Therefore Q and DEG must be chosen such that (Q-1) * (Q-1) * (DEG+1) + (Q-1) < p, where p is the modulus of the circuit field.

            // constrain_poly_reduction_by_cyclo::<{ 2 * DEG - 2 }, DEG, Q, F>(
            //     &pk0_u,
            //     cyclo_with_length.clone(),
            //     quotient_0_with_length,
            //     quotient_0_times_cyclo_with_length,
            //     remainder_0_with_length.clone(),
            //     range,
            //     ctx_gate,
            //     ctx_rlc,
            //     rlc,
            // );

            // let pk0_u = remainder_0_with_length.get_poly().clone();

            // assert_eq!(pk0_u.len() - 1, 2 * DEG - 2);

            // // pk0_u is a polynomial of degree 2*DEG - 2
            // // pk0_u now has coefficients in the [0, Q-1] range

            // // But actually, the degree of pk0_u should be DEG - 1 after reduction by the cyclo polynomial, the first DEG - 1 coefficients are just zeroes

            // // 1.4 Enforce that the first DEG - 1 coefficients of pk0_u are zeroes
            // for pk0_u_element in pk0_u.iter().take(DEG - 1) {
            //     let bool = range
            //         .gate()
            //         .is_equal(ctx_gate, *pk0_u_element, Constant(F::from(0)));
            //     range.gate().assert_is_const(ctx_gate, &bool, &F::from(1));
            // }

            // // Therefore, we can safely trim the first DEG - 1 coefficients from pk0_u

            // let pk0_u_trimmed: Vec<_> = pk0_u.iter().skip(DEG - 1).cloned().collect();

            // // assert that the degree of pk0_u_trimmed is DEG - 1
            // assert_eq!(pk0_u_trimmed.len() - 1, DEG - 1);

            // // pk0_u_trimmed is a polynomial in the R_q ring!

            // // 1.5 m * delta

            // // Perform the polynomial scalar multiplication between m and delta.
            // // The assumption of the `poly_scalar_mul` chip is that the coefficients of the input polynomials are constrained such to avoid overflow during the polynomial scalar multiplication
            // // m has coefficients in the [0, T/2] OR [Q - T/2, Q - 1] range according to the constraint set above
            // // delta is a constant equal to Q/T rounded to the lower integer from BFV paper
            // // The coefficients of m_delta are in the [0, (Q-1) * (Q/T)] range. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Scalar Multiplication section
            // // Q and T must be chosen such that (Q-1) * (Q/T) < p, where p is the modulus of the circuit field.

            // let m_delta = poly_scalar_mul::<{ DEG - 1 }, F>(ctx_gate, &m, &delta, range.gate());

            // // m_delta is a polynomial of degree DEG - 1

            // // 1.6 pk0_u_trimmed + m_delta

            // // Perform the polynomial addition between pk0_u_trimmed and m_delta.
            // // The assumption of the `poly_add` chip is that the coefficients of the input polynomials are constrained such to avoid overflow during the polynomial addition
            // // `m_delta` has coefficients in the [0, (Q-1) * (Q/T)] range according to the constraint set above
            // // `pk0_u_trimmed` has coefficients in the [0, Q-1] range according to the constraint set above
            // // The coefficients of pk0_u_trimmed_plus_m_delta are in the [0, (Q-1) * (Q-T) + (Q-1)] range. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Addition section
            // // Q and T must be chosen such that (Q-1) * (Q-T) + (Q-1) < p, where p is the modulus of the circuit field.
            // let pk0_u_trimmed_plus_m_delta =
            //     poly_add::<{ DEG - 1 }, F>(ctx_gate, &pk0_u_trimmed, &m_delta, range.gate());

            // // Note: Addition does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1
            // // pk0_u_trimmed_plus_m_delta is a polynomial in the R_q ring

            // // 1.7 c0 = pk0_u_trimmed_plus_m_delta + e0

            // // Perform the polynomial addition between pk0_u_trimmed_plus_m_delta and e0.
            // // The assumption of the `poly_add` chip is that the coefficients of the input polynomials are constrained such to avoid overflow during the polynomial addition
            // // `pk0_u_trimmed_plus_m_delta` has coefficients in the [0, (Q-1) * (Q-T) + (Q-1)] range according to the constraint set above
            // // `e0` has coefficients in the [0, B] OR [Q-B, Q-1] range according to the constraint set above
            // // The coefficients of computed_c0 are in the [0, (Q-1) * (Q-T) + (Q-1) + (Q-1)] range. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Addition section
            // // Q, T must be chosen such that (Q-1) * (Q-T) + (Q-1) + (Q-1) < p, where p is the modulus of the circuit field.
            // let computed_c0 =
            //     poly_add::<{ DEG - 1 }, F>(ctx_gate, &pk0_u_trimmed_plus_m_delta, &e0, range.gate());

            // // get the number of bits needed to represent the value of (Q-1) * (Q-T) + (Q-1) + (Q-1)
            // let q_minus_1 = BigInt::from(Q) - BigInt::from(1u32);
            // let q_minus_t = BigInt::from(Q) - BigInt::from(T as u32);
            // let binary_representation = format!(
            //     "{:b}",
            //     (q_minus_1.clone() * q_minus_t + q_minus_1.clone() + q_minus_1)
            // );
            // let num_bits_3 = binary_representation.len();

            // // Coefficients of computed_c0 are in the [0, (Q-1) * (Q-T) + (Q-1) + (Q-1)] range according to the polynomial addition constraint set above.
            // // Therefore the coefficients of computed_c0 are known to have <= `num_bits_4` bits, therefore they satisfy the assumption of the `poly_reduce_by_modulo_q` chip

            // // 1.8 Reduce the coefficients of `pk0_u_trimmed_plus_m_delta` by modulo `Q`

            // let computed_c0 =
            //     poly_reduce_by_modulo_q::<{ DEG - 1 }, Q, F>(ctx_gate, &computed_c0, range, num_bits_3);

            // // Note: Addition does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1
            // // computed_c0 is a polynomial in the R_q ring!

            // // 2. COMPUTE C1 (c1 is the second ciphertext component)

            // // 2.1 pk1 * u
            // // Constrain the polynomial multiplication between pk1 and u to be equal to pk§_u using the `constrain_poly_mul` chip

            // constrain_poly_mul(
            //     pk1_with_length,
            //     u_with_length,
            //     pk1_u_with_length,
            //     ctx_gate,
            //     ctx_rlc,
            //     rlc,
            // );

            // // pk1_u is a polynomial of degree (DEG - 1) * 2 = 2*DEG - 2
            // // pk1_u has coefficients in the [0, (Q-1) * (Q-1) * (DEG+1)] range. This is the maximum value that a coefficient of pk0_u can take. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Multiplication section

            // // 2.2 Reduce the coefficients by modulo `Q`

            // // The coefficients of pk1_u are in the range [0, (Q-1) * (Q-1) * (DEG+1)] according to the above analysis.
            // // Therefore the coefficients of pk1_u are known to have <= `num_bits_1` bits, therefore satisfying the assumption of the `poly_reduce_by_modulo_q` chip

            // let pk1_u =
            //     poly_reduce_by_modulo_q::<{ 2 * DEG - 2 }, Q, F>(ctx_gate, &pk1_u, range, num_bits_1);

            // // pk0_u is a polynomial of degree (DEG - 1) * 2 = 2 * DEG - 2
            // // pk0_u now has coefficients in the [0, Q-1] range after reduction by modulo Q
            // // cyclo is a polynomial of degree DEG

            // // 2.3 Reduce pk1_u by the cyclo polynomial

            // for i in 0..DEG - 1 {
            //     range.check_less_than_safe(ctx_gate, quotient_1_with_length.get_poly()[i], Q);
            // }

            // for i in 0..2 * DEG - 1 {
            //     range.check_less_than_safe(ctx_gate, remainder_1_with_length.get_poly()[i], Q);
            // }

            // // For the analysis of the assumptions of the `constrain_poly_reduction_by_cyclo` chip, see the comments above in the `1.3 Reduce pk0_u by the cyclo polynomial` section
            // constrain_poly_reduction_by_cyclo::<{ 2 * DEG - 2 }, DEG, Q, F>(
            //     &pk1_u,
            //     cyclo_with_length.clone(),
            //     quotient_1_with_length,
            //     quotient_1_times_cyclo_with_length,
            //     remainder_1_with_length.clone(),
            //     range,
            //     ctx_gate,
            //     ctx_rlc,
            //     rlc,
            // );

            // let pk1_u = remainder_1_with_length.get_poly().clone();

            // // assert that the degree of pk1_u is 2*DEG - 2

            // // pk1_u is a polynomial of degree 2*DEG - 2
            // // pk1_u now has coefficients in the [0, Q-1] range

            // // But actually, the degree of pk1_u should be DEG - 1 after reduction by the cyclo polynomial, the first DEG - 1 coefficients are just zeroes

            // // 2.4 Enforce that the first DEG - 1 coefficients of pk0_u are zeroes
            // for pk1_u_element in pk1_u.iter().take(DEG - 1) {
            //     let bool = range
            //         .gate()
            //         .is_equal(ctx_gate, *pk1_u_element, Constant(F::from(0)));
            //     range.gate().assert_is_const(ctx_gate, &bool, &F::from(1));
            // }

            // // Therefore, we can safely trim the first DEG - 1 coefficients from pk1_u

            // let pk1_u_trimmed: Vec<_> = pk1_u.iter().skip(DEG - 1).cloned().collect();

            // // assert that the degree of pk1_u_trimmed is DEG - 1
            // assert_eq!(pk1_u_trimmed.len() - 1, DEG - 1);

            // // pk1_u_trimmed is a polynomial in the R_q ring!

            // // 2.5 c1 = pk1_u_trimmed + e1

            // // Perform the polynomial addition between pk1_u_trimmed and e1.
            // // The assumption of the `poly_add` chip is that the coefficients of the input polynomials are constrained such to avoid overflow during the polynomial addition
            // // `pk1_u_trimmed` has coefficients in the [0, Q-1] range according to the constraint set above
            // // `e1` has coefficients in the [0, B] OR [Q-B, Q-1] range according to the constraint set above
            // // The coefficients of computed_c1 are in the [0, 2Q - 2] range. Why? Answer is here -> https://hackmd.io/@letargicus/Bk4KtYkSp - Polynomial Addition section
            // // Q must be chosen such that 2Q - 2 < p, where p is the modulus of the circuit field.
            // let computed_c1 = poly_add::<{ DEG - 1 }, F>(ctx_gate, &pk1_u_trimmed, &e1, range.gate());

            // // get the number of bits needed to represent the value of 2Q - 2
            // // get the number of bits needed to represent the value of 2Q - 2
            // let binary_representation = format!("{:b}", (2 * Q - 2));
            // let num_bits_4 = binary_representation.len();

            // // The coefficients of computed_c1 are in the range [0, 2Q - 2] according to the polynomial addition performed above.
            // // Therefore the coefficients of computed_c1 are known to have <= `num_bits_4` bits, therefore they satisfy the assumption of the `poly_reduce_by_modulo_q` chip

            // // 2.6 Reduce the coefficients of `computed_c1` by modulo `Q`

            // let computed_c1 =
            //     poly_reduce_by_modulo_q::<{ DEG - 1 }, Q, F>(ctx_gate, &computed_c1, range, num_bits_4);

            // // Note: Addition does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1
            // // computed_c1 is a polynomial in the R_q ring

            // // Enforce equality between `c0` and `computed_c0` using equality check
            // // Enfroce equality between `c1` and `computed_c1` using equality check
            // for i in 0..DEG {
            //     let bool_0 = range.gate().is_equal(ctx_gate, c0[i], computed_c0[i]);
            //     range.gate().assert_is_const(ctx_gate, &bool_0, &F::from(1));

            //     let bool_1 = range.gate().is_equal(ctx_gate, c1[i], computed_c1[i]);
            //     range.gate().assert_is_const(ctx_gate, &bool_1, &F::from(1));
            // }
        };

    callback
}

fn main() {
    env_logger::init();

    let args = Cli::parse();

    run_eth(bfv_encryption_circuit, args);
}
