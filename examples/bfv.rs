use axiom_eth::keccak::KeccakChip;
use axiom_eth::{EthChip, Field};
use clap::Parser;
use halo2_base::safe_types::{GateInstructions, RangeInstructions};
use halo2_base::{AssignedValue, Context, QuantumCell::Constant};
use halo2_scaffold::scaffold::{cmd::Cli, run_eth};
use serde::{Deserialize, Serialize};
use zk_fhe::chips::poly_distribution::{
    check_poly_coefficients_in_range, check_poly_from_distribution_chi_key,
};
use zk_fhe::chips::poly_operations::{
    constrain_poly_mul, poly_add, poly_divide_by_cyclo, poly_reduce, poly_scalar_mul,
};
use zk_fhe::chips::utils::poly_mul;

/// Circuit inputs for BFV encryption operations
///
/// # Type Parameters
///
/// * `DEG`: Degree of the cyclotomic polynomial `cyclo` of the polynomial ring R_q.
/// * `Q`: Modulus of the cipher text field
/// * `T`: Modulus of the plaintext field
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
/// * `c0`: First ciphertext component polynomial coefficients of degree DEG-1 [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term
/// * `c1`: Second ciphertext component polynomial coefficients of degree DEG-1 [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term

/// # Assumptions (to be checked outside the circuit)
///
/// * `DEG` must be a power of 2
/// * `Q` must be a prime number and be greater than 1.
/// * `Q` must be less than 2^64
/// * `Q` is less than (Q-1) * (Q-1) * DEG < p where p is the prime field of the circuit
/// * `T` must be a prime number and must be greater than 1 and less than `Q`
/// * `B` must be a positive integer and must be less than `Q`
/// * `cyclo` must be the cyclotomic polynomial of degree `DEG` => x^DEG + 1 (this is a public output of the circuit)
/// * `pk0` and `pk1` must be polynomials in the R_q ring. The ring R_q is defined as R_q = Z_q[x]/(x^DEG + 1)

// For real world applications, the parameters should be chosen according to the security level required.
// DEG and Q Parameters of the BFV encryption scheme should be chosen according to TABLES of RECOMMENDED PARAMETERS for 128-bits security level
// https://homomorphicencryption.org/wp-content/uploads/2018/11/HomomorphicEncryptionStandardv1.1.pdf
// B is the upper bound of the distribution Chi Error. Pick standard deviation 𝜎 ≈ 3.2 according to the HomomorphicEncryptionStandardv1 paper.
// T should be be picked according to Lattigo (https://github.com/tuneinsight/lattigo/blob/master/bfv/params.go) implementation
// As suggest by https://eprint.iacr.org/2021/204.pdf (paragraph 2) B = 6σerr
// These are just parameters used for fast testing
const DEG: usize = 4;
const Q: u64 = 4637;
const T: u64 = 7;
const B: u64 = 18;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CircuitInput<const DEG: usize, const Q: u64, const T: u64, const B: u64> {
    pub pk0: Vec<u64>, // polynomial coefficients [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term. PUBLIC. Should live in R_q according to assumption
    pub pk1: Vec<u64>, // polynomial coefficients [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term. PUBLIC. Should live in R_q according to assumption
    pub m: Vec<u64>, // polynomial coefficients [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term. PRIVATE. Should in R_t (enforced inside the circuit)
    pub u: Vec<u64>, // polynomial coefficients [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term. PRIVATE. Should live in R_q and be sampled from the distribution ChiKey (checked inside the circuit)
    pub e0: Vec<u64>, // polynomial coefficients [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term. PRIVATE. Should live in R_q and be sampled from the distribution ChiError (checked inside the circuit)
    pub e1: Vec<u64>, // polynomial coefficients [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term. PRIVATE. Should live in R_q and be sampled from the distribution ChiError (checked inside the circuit)
    pub c0: Vec<u64>, // polynomial coefficients [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term. Should live in R_q. This is just a test value compared to the ciphertext c0 generated as (public) output by the circuit
    pub c1: Vec<u64>, // polynomial coefficients [a_DEG-1, a_DEG-2, ..., a_1, a_0] where a_0 is the constant term. Should live in R_q. This is just a test value compared to the ciphertext c1 generated as (public) output by the circuit
}

fn bfv_encryption_circuit<F: Field>(
    ctx: &mut Context<F>,
    _eth_chip: &EthChip<F>,
    _keccak: &mut KeccakChip<F>,
    input: CircuitInput<DEG, Q, T, B>,
    make_public: &mut Vec<AssignedValue<F>>,
) -> impl FnOnce(&mut Context<F>, &mut Context<F>, &EthChip<F>) + Clone {
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

    // The circuit logic requires to access some random value
    // In order to draw randomness within the circuit we use Axiom's Challenge API (https://hackmd.io/@axiom/SJw3p-qX3)
    // Challenge API requires a Phase 0 of witness generation. A commitment from the witness generated during Phase 0 is then hashed to generate the random value according to Fiat-Shamir heuristic.
    // This random challenge can be then used as part of witness generation during Phase 1. We will need this to perform efficient polynomial multiplication.

    // Phase 0: Assign the input polynomials to the circuit witness table
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

    // Assign the cyclotomic polynomial to the circuit -> x^DEG + 1
    // Performing the assignemnt for the index 0, using a for loop from 1 to DEG - 1, and performing the assignment for the index DEG enforces that:
    // - the degree of the polynomial is DEG
    // - the leading coefficient is 1
    // - the constant term is 1
    // - all the other coefficients are 0
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

    // Assign the length of the input polynomials pk0 and u to the circuit
    let pk0_len = ctx.load_witness(F::from(input.pk0.len() as u64));
    let u_len = ctx.load_witness(F::from(input.u.len() as u64));

    // Compute the polynomial pk0 * u outside the circuit
    let pk0_u_u64 = poly_mul(&input.pk0, &input.u);

    // Assign the polynomial pk0_u_u64 to the circuit
    let mut pk0_u = vec![];

    // This should be a polynomial of degree 2*(DEG - 1)
    for &item in pk0_u_u64.iter().take(2 * (DEG - 1) + 1) {
        let pk0_u_val = F::from(item);
        let pk0_u_assigned_value = ctx.load_witness(pk0_u_val);
        pk0_u.push(pk0_u_assigned_value);
    }

    assert!(pk0_u.len() - 1 == 2 * (DEG - 1));

    // Assign the length of the polynomial pk0_u to the circuit
    let pk0_u_len = ctx.load_witness(F::from(pk0_u.len() as u64));

    // Assign the length of the input polynomials pk1 and u to the circuit
    let pk1_len = ctx.load_witness(F::from(input.pk1.len() as u64));

    // Compute the polynomial pk1 * u outside the circuit
    let pk1_u_u64 = poly_mul(&input.pk1, &input.u);

    // Assign the polynomial pk1_u_u64 to the circuit
    let mut pk1_u = vec![];

    // This should be a polynomial of degree 2*(DEG - 1)
    for &item in pk1_u_u64.iter().take(2 * (DEG - 1) + 1) {
        let pk1_u_val = F::from(item);
        let pk1_u_assigned_value = ctx.load_witness(pk1_u_val);
        pk1_u.push(pk1_u_assigned_value);
    }

    assert!(pk1_u.len() - 1 == 2 * (DEG - 1));

    // Assign the length of the polynomial pk1_u to the circuit
    let pk1_u_len = ctx.load_witness(F::from(pk1_u.len() as u64));

    // // Expose to the public the coefficients of c0 and c1
    // for i in 0..DEG {
    //     make_public.push(c0[i]);
    // }

    // for i in 0..DEG {
    //     make_public.push(c1[i]);
    // }

    // // Expose to the public pk0 and pk1
    // for i in 0..DEG {
    //     make_public.push(pk0[i]);
    // }

    // for i in 0..DEG {
    //     make_public.push(pk1[i]);
    // }

    // // Expose to the public `cyclo`
    // for i in 0..DEG + 1 {
    //     make_public.push(cyclo[i]);
    // }

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
            - e0 must be sampled from the distribution ChiError

            Approach:
            - `check_poly_from_distribution_chi_error` chip guarantees that the coefficients of e0 are in the range [0, b] OR [q-b, q-1]
            - As this range is a subset of the [0, Q-1] range, the coefficients of e0 are guaranteed to be in the [0, Q-1] range
            - The assignment for loop above enforces that the degree of e0 is DEG - 1
        */

        /* constraint on e1
            Same as e0
        */

        // Assumption for the chip is that B < Q which is satisfied by circuit assumption
        check_poly_coefficients_in_range::<{ DEG - 1 }, Q, B, F>(ctx_gate, &e0, range);
        check_poly_coefficients_in_range::<{ DEG - 1 }, Q, B, F>(ctx_gate, &e1, range);

        /* constraint on u
            - u must be a polynomial in the R_q ring => Coefficients must be in the [0, Q-1] range and the degree of u must be DEG - 1
            - u must be sampled from the distribution ChiKey

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

        // // pk0 * u

        // // Perform the polynomial multiplication between pk0 and u.

        // // DEGREE ANALYSIS
        // // The degree of pk0 is DEG - 1 according to the constraint set above
        // // The degree of u is DEG - 1 according to the constraint set above
        // // The degree of pk0_u is constrained to be DEG - 1 + DEG - 1 = 2*DEG - 2 according to the logic of the `poly_mul_equal_deg` chip

        // // COEFFICIENTS OVERFLOW ANALYSIS
        // // The coefficients of pk0 are in the range [0, Q-1] according to the assumption of the circuit
        // // The coefficients of u are either 0, 1 or Q-1 according to the constraint set above.
        // // The coefficients of pk0_u are calculated as $c_{k} = \sum_{i=0}^{k} pk0[i] * u[k - i]$. Where k is the index of the coefficient c of pk0_u.
        // // For two polynomials of the same degree n, the maximum number of multiplications in the sum is for k = n. Namely for the coefficient c_n.
        // // The number of multiplications in the sum for the coefficient c_n is n + 1.
        // // Given that the input polynomials are of degree DEG - 1, the maximum number of multiplications in the sum is for k = DEG - 1.
        // // In that case there are max DEG multiplications in the sum.
        // // It follows that the maximum value that a coefficient of pk0_u can have is (Q-1) * (Q-1) * DEG.
        // // Q needs to be chosen such that (Q-1) * (Q-1) * DEG < p where p is the prime field of the circuit in order to avoid overflow during the polynomial multiplication.
        // // (Q-1) * (Q-1) * DEG < p according to the assumption of the circuit.

        // let pk0_u =
        //     poly_mul_equal_deg::<{ DEG - 1 }, F>(ctx_gate, pk0.clone(), u.clone(), range.gate());

        // Enforce pk0_u = pk0 * u using `constrain_poly_mul` gate
        constrain_poly_mul(
            pk0,
            pk0_len,
            u.clone(),
            u_len,
            pk0_u.clone(),
            pk0_u_len,
            ctx_gate,
            ctx_rlc,
            rlc,
            range.gate(),
        );

        // pk0_u is a polynomial of degree (DEG - 1) * 2 = 2*DEG - 2
        // pk0_u has coefficients in the [0, (Q-1) * (Q-1) * DEG] range

        // Reduce the coefficients by modulo `Q`

        // get the number of bits needed to represent the value of (Q-1) * (Q-1) * DEG

        let binary_representation = format!("{:b}", ((Q - 1) * (Q - 1) * (DEG as u64)));
        let num_bits_1 = binary_representation.len();

        // The coefficients of pk0_u are in the range [0, (Q-1) * (Q-1) * DEG] according to the polynomial multiplication constraint set above.
        // Therefore the coefficients of pk0_u are known to have <= `num_bits_1` bits, therefore they satisfy the assumption of the `poly_reduce` chip

        let pk0_u = poly_reduce::<{ 2 * DEG - 2 }, Q, F>(ctx_gate, &pk0_u, range, num_bits_1);

        // pk0_u is a polynomial of degree (DEG - 1) * 2 = 2*DEG - 2
        // pk0_u has coefficients in the [0, Q-1] range
        // cyclo is a polynomial of degree DEG
        // Reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1 to get a polynomial of degree DEG - 1

        // Dealing with the assumption of the `poly_divide_by_cyclo` chip
        // - The degree of dividend (pk0_u) is equal to (2 * DEG) - 2 according to the constraint set above
        // - The coefficients of dividend are in the [0, Q-1] range according to the constraint set above
        // - The divisor is a cyclotomic polynomial of degree DEG with coefficients either 0 or 1
        // - The coefficients of dividend and divisor can be expressed as u64 values as long as Q - 1 is less than 2^64
        // - Q is chosen such that (Q-1) * (2 * DEG - 2 - DEG + 1)] + Q-1 < p. Note that this is a subset of the condition (Q-1) * (Q-1) * DEG < p which is an assumption of the circuit.

        let pk0_u =
            poly_divide_by_cyclo::<{ 2 * DEG - 2 }, DEG, Q, F>(ctx_gate, &pk0_u, &cyclo, range);

        // assert that the degree of pk0_u is 2*DEG - 2

        assert_eq!(pk0_u.len() - 1, 2 * DEG - 2);

        // But actually, the degree of pk0_u is DEG - 1, the first DEG - 1 coefficients are just zeroes

        // Enforce that the first DEG - 1 coefficients of pk0_u are zeroes

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

        // m * delta

        // Perform the polynomial scalar multiplication between m and delta.

        // DEGREE ANALYSIS
        // The degree of m is DEG - 1 according to the constraint set above
        // Delta is a scalar constant
        // The degree of m_delta is constrained to be DEG - 1 according to the logic of the `poly_scalar_mul` chip

        // COEFFICIENTS OVERFLOW ANALYSIS
        // The coefficients of m are in the range [0, T/2] OR [Q - T/2, Q - 1] according to the constaints set above.
        // Delta is a constant equal to Q/T (integer division) where T < Q according to the assumption of the circuit.
        // The maximum value of a coefficient of m_delta is (Q-1) * (Q/T)
        // If the condition (Q-1) * (Q/T) < p is satisfied there is no risk of overflow during the scalar multiplication.
        // Note that this condition is a subset of the condition (Q-1) * (Q-1) * DEG < p which is an assumption of the circuit.

        let m_delta = poly_scalar_mul::<{ DEG - 1 }, F>(
            ctx_gate,
            &m,
            &Constant(F::from(DELTA)),
            range.gate(),
        );

        // m_delta is a polynomial of degree DEG - 1
        // Coefficients of m_delta are in the [0, (Q-1) * (Q/T)] range

        // Reduce the coefficients of `m_delta` by modulo `Q`

        // get the number of bits needed to represent the value of (Q-1) * (Q/T)

        let binary_representation = format!("{:b}", ((Q - 1) * (Q / T)));
        let num_bits_2 = binary_representation.len();

        // The coefficients of m_delta are in the range [0,  (Q-1) * (Q/T)] according to the polynomial scalar multiplication constraint set above.
        // Therefore the coefficients of m_delta are known to have <= `num_bits_2` bits, therefore they satisfy the assumption of the `poly_reduce` chip

        let m_delta = poly_reduce::<{ DEG - 1 }, Q, F>(ctx_gate, &m_delta, range, num_bits_2);

        // Note: Scalar multiplication does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1
        // m_delta is a polynomial in the R_q ring

        // pk0_u_trimmed + m_delta

        // Perform the polynomial addition between pk0_u_trimmed and m_delta.

        // DEGREE ANALYSIS
        // The degree of pk0_u_trimmed is DEG - 1 according to the constraint set above
        // The degree of m_delta is DEG - 1 according to the constraint set above
        // The degree of pk0_u_trimmed_plus_m_delta is constrained to be DEG - 1 according to the logic of the `poly_add` chip

        // COEFFICIENTS OVERFLOW ANALYSIS
        // The coefficients of pk0_u_trimmed and m_delta are in the [0, Q -1] range according to the constraint set above.
        // The coefficients of m_delta are in the [0, Q -1] range according to the constraint set above.
        // The maximum value of the coefficient of pk0_u_trimmed_plus_m_delta is (Q-1) + (Q-1) = 2Q - 2.
        // If the condition  (Q-1) + (Q-1) < p is satisfied there is no risk of overflow during the polynomial addition.
        // Note that this condition is a subset of the condition (Q-1) * (Q-1) * DEG < p which is an assumption of the circuit.

        let pk0_u_trimmed_plus_m_delta =
            poly_add::<{ DEG - 1 }, F>(ctx_gate, &pk0_u_trimmed, &m_delta, range.gate());

        // Reduce the coefficients of `m_delta` by modulo `Q`
        // Coefficients of pk0_u_trimmed_plus_m_delta are in the [0, 2Q - 2] range

        // get the number of bits needed to represent the value of 2Q - 2

        let binary_representation = format!("{:b}", (2 * Q - 2));
        let num_bits_3 = binary_representation.len();

        // The coefficients of pk0_u_trimmed_plus_m_delta are in the range [0, 2Q - 2] according to the polynomial addition constraint set above.
        // Therefore the coefficients of m_delta are known to have <= `num_bits_3` bits, therefore they satisfy the assumption of the `poly_reduce` chip

        // Reduce the coefficients of `pk0_u_trimmed_plus_m_delta` by modulo `Q`

        let pk0_u_trimmed_plus_m_delta = poly_reduce::<{ DEG - 1 }, Q, F>(
            ctx_gate,
            &pk0_u_trimmed_plus_m_delta,
            range,
            num_bits_3,
        );

        // Note: Addition does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1
        // pk0_u_trimmed_plus_m_delta is a polynomial in the R_q ring

        // c0 = pk0_u_trimmed_plus_m_delta + e0

        // Perform the polynomial addition between pk0_u_trimmed_plus_m_delta and e0.

        // DEGREE ANALYSIS
        // The degree of pk0_u_trimmed_plus_m_delta is DEG - 1 according to the constraint set above
        // The degree of e0 is DEG - 1 according to the constraint set above
        // The degree of c0 is constrained to be DEG - 1 according to the logic of the `poly_add` chip

        // COEFFICIENTS OVERFLOW ANALYSIS
        // The coefficients of pk0_u_trimmed_plus_m_delta and m_delta are in the [0, Q -1] range according to the constraint set above.
        // The cofficients of e0 are in the range [0, b] OR [q-b, q-1] according to the constraint set above.
        // The maximum value of the coefficient of c0 is (Q-1) + (Q-1) = 2Q - 2.
        // If the condition (Q-1) + (Q-1) < p is satisfied there is no risk of overflow during the polynomial addition.
        // Note that this condition is a subset of the condition (Q-1) * (Q-1) * DEG < p which is an assumption of the circuit.

        let c0 =
            poly_add::<{ DEG - 1 }, F>(ctx_gate, &pk0_u_trimmed_plus_m_delta, &e0, range.gate());

        // The coefficients of c0 are in the range [0, 2Q - 2] according to the polynomial addition constraint set above.
        // Therefore the coefficients of c0 are known to have <= `num_bits_3` bits, therefore they satisfy the assumption of the `poly_reduce` chip

        // Reduce the coefficients of `pk0_u_trimmed_plus_m_delta` by modulo `Q`

        let c0 = poly_reduce::<{ DEG - 1 }, Q, F>(ctx_gate, &c0, range, num_bits_3);

        // Note: Addition does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1
        // c0 is a polynomial in the R_q ring!

        // 1. COMPUTE C1 (c1 is the second ciphertext component)

        // pk1 * u

        // Perform the polynomial multiplication between pk1 and u.

        // DEGREE ANALYSIS
        // The degree of pk1 is DEG - 1 according to the constraint set above
        // The degree of u is DEG - 1 according to the constraint set above
        // The degree of pk1_u is constrained to be DEG - 1 + DEG - 1 = 2*DEG - 2 according to the logic of the `poly_mul_equal_deg` chip

        // COEFFICIENTS OVERFLOW ANALYSIS
        // The coefficients of pk1 are in the range [0, Q-1] according to the assumption of the circuit
        // The coefficients of u are either 0, 1 or Q-1 according to the constraint set above.
        // The coefficients of pk1_u are calculated as $c_{k} = \sum_{i=0}^{k} pk1[i] * u[k - i]$. Where k is the index of the coefficient c of pk1_u.
        // For two polynomials of the same degree n, the maximum number of multiplications in the sum is for k = n. Namely for the coefficient c_n.
        // The number of multiplications in the sum for the coefficient c_n is n + 1.
        // Given that the input polynomials are of degree DEG - 1, the maximum number of multiplications in the sum is for k = DEG - 1.
        // In that case there are max DEG multiplications in the sum.
        // It follows that the maximum value that a coefficient of pk1_u can have is (Q-1) * (Q-1) * DEG.
        // Q needs to be chosen such that (Q-1) * (Q-1) * DEG < p where p is the prime field of the circuit in order to avoid overflow during the polynomial multiplication.
        // (Q-1) * (Q-1) * DEG < p according to the assumption of the circuit.

        // let pk1_u =
        //     poly_mul_equal_deg::<{ DEG - 1 }, F>(ctx_gate, pk1.clone(), u.clone(), range.gate());

        // Enforce pk1_u = pk1 * u using `constrain_poly_mul` gate
        constrain_poly_mul(
            pk1,
            pk1_len,
            u.clone(),
            u_len,
            pk1_u.clone(),
            pk1_u_len,
            ctx_gate,
            ctx_rlc,
            rlc,
            range.gate(),
        );

        // pk1_u is a polynomial of degree (DEG - 1) * 2 = 2*DEG - 2
        // pk1_u has coefficients in the [0, (Q-1) * (Q-1) * DEG] range

        // Reduce the coefficients by modulo `Q`

        // The coefficients of pk1_u are in the range [0, (Q-1) * (Q-1) * DEG] according to the polynomial multiplication constraint set above.
        // Therefore the coefficients of pk1_u are known to have <= `num_bits_1` bits, therefore they satisfy the assumption of the `poly_reduce` chip

        let pk1_u = poly_reduce::<{ 2 * DEG - 2 }, Q, F>(ctx_gate, &pk1_u, range, num_bits_1);

        // pk1_u is a polynomial of degree (DEG - 1) * 2 = 2*DEG - 2
        // pk1_u has coefficients in the [0, Q-1] range
        // cyclo is a polynomial of degree DEG
        // Reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1 to get a polynomial of degree DEG - 1

        // Dealing with the assumption of the `poly_divide_by_cyclo` chip
        // - The degree of dividend (pk0_1) is equal to (2 * DEG) - 2 according to the constraint set above
        // - The coefficients of dividend are in the [0, Q-1] range according to the constraint set above
        // - The divisor is a cyclotomic polynomial of degree DEG with coefficients either 0 or 1
        // - The coefficients of dividend and divisor can be expressed as u64 values as long as Q - 1 is less than 2^64
        // - Q is chosen such that (Q-1) * (2 * DEG - 2 - DEG + 1)] + Q-1 < p. Note that this is a subset of the condition (Q-1) * (Q-1) * DEG < p which is an assumption of the circuit.

        let pk1_u =
            poly_divide_by_cyclo::<{ 2 * DEG - 2 }, DEG, Q, F>(ctx_gate, &pk1_u, &cyclo, range);

        // assert that the degree of pk1_u is 2*DEG - 2

        assert_eq!(pk1_u.len() - 1, 2 * DEG - 2);

        // But actually, the degree of pk1_u is DEG - 1, the first DEG - 1 coefficients are just zeroes

        // Enforce that the first DEG - 1 coefficients of pk1_u are zeroes

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

        // c1 = pk1_u_trimmed + e1

        // Perform the polynomial addition between pk1_u_trimmed and e1.

        // DEGREE ANALYSIS
        // FIX The degree of pk1_u_trimmed is DEG - 1 according to the constraint set above
        // The degree of e1 is DEG - 1 according to the constraint set above
        // The degree of c1 is constrained to be DEG - 1 according to the logic of the `poly_add` chip

        // COEFFICIENTS OVERFLOW ANALYSIS
        // The coefficients of pk1_u_trimmed are in the [0, Q-1] range according to the constraint set above.
        // The cofficients of e1 are in the range [0, b] OR [q-b, q-1] according to the constraint set above.
        // The maximum value of the coefficient of c1 is (Q-1) + (Q-1) = 2Q - 2.
        // If the condition (Q-1) + (Q-1) < p is satisfied there is no risk of overflow during the polynomial addition.
        // Note that this condition is a subset of the condition (Q-1) * (Q-1) * DEG < p which is an assumption of the circuit.

        // Perform the polynomial addition between pk1_u_trimmed and e1.

        let c1 = poly_add::<{ DEG - 1 }, F>(ctx_gate, &pk1_u_trimmed, &e1, range.gate());

        // The coefficients of c1 are in the range [0, 2Q - 2] according to the polynomial addition constraint set above.
        // Therefore the coefficients of c1 are known to have <= `num_bits_3` bits, therefore they satisfy the assumption of the `poly_reduce` chip

        // Reduce the coefficients of `c1` by modulo `Q`

        let c1 = poly_reduce::<{ DEG - 1 }, Q, F>(ctx_gate, &c1, range, num_bits_3);

        // Note: Addition does not change the degree of the polynomial, therefore we do not need to reduce the coefficients by the cyclotomic polynomial of degree `DEG` => x^DEG + 1
        // c1 is a polynomial in the R_q ring

        // That that c0 and c1 computed inside the circuit are equal to the ciphertexts provided as input in the test vector json file
        // Check outside the circuit that the remainder matches the expected one
        for i in 0..DEG {
            assert_eq!(*c0[i].value(), F::from(input.c0[i]));
            assert_eq!(*c1[i].value(), F::from(input.c1[i]));
        }
    };

    callback
}

fn main() {
    env_logger::init();

    let args = Cli::parse();

    run_eth(bfv_encryption_circuit, args);
}
