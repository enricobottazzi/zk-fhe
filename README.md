# zk-fhe
Zk proving the correct execution of encryption operation under BFV Fully Homomorphic Encryption scheme

Implementation based on [Revisiting Homomorphic Encryption Schemes for Finite Fields](https://eprint.iacr.org/2021/204.pdf)

The application is not production ready and is only meant to be used for educational purposes.

## Disclaimer

This is a research project and is not meant to be used in production. The code is not audited.

## Guide

Many polynomial operations performed inside the [circuit](./examples/bfv.rs) involve careful handling of coefficients in order to avoid overflows. A thorough explanation of the range of expansion of polynomial coefficients during this operations can be found in [Polynomial Coefficients and Degree Analysis](https://hackmd.io/@letargicus/Bk4KtYkSp)

## Quick Start

**Mock Prover**

`cargo run --example bfv -- --name bfv -k 9 --input bfv.in mock`

The `MockProver` does not run the cryptographic prover on your circuit, but instead directly checks if constraints are satisfied. This is useful for testing purposes, and runs faster than the actual prover.

- `LOOKUP_BITS`, in the backend build a lookup table filled with value in the range [0, 2**LOOKUP_BITS)
- `bfv` is the name of the circuit located in `examples/bfv.rs` 
- `bfv.in` is the input file for the circuit located in `data/bfv.in`. This test vector file can be generated for different encryption using [bfv-py](https://github.com/yuriko627/bfv-py)
- `-k` is the DEGREE of the circuit as you specify to set the circuit to have `2^k` number of rows. The number of rows is determined by the number of constraints in the circuit. Working with larger data inputs will require a larger degree.

**Key Generation**

`cargo run --example bfv -- --name bfv -k 9 --input bfv.in keygen`

To generate a random universal trusted setup (for testing only!) and the proving and verifying keys for your circuit.

For technical reasons (related to [axiom Halo2-scaffold](https://github.com/axiom-crypto/halo2-scaffold)), keygen still requires an input file of the correct format. You can use the same input file as for the prover. But be aware that the actual input data are not encoded in the key generation. 

This will generate a proving key `data/bfv.pk` and a verifying key `data/bfv.vk`. It will also generate a file `configs/bfv.json` which describes (and pins down) the configuration of the circuit. This configuration file is later read by the prover.

**Proof Generation**

`cargo run --example bfv -- --name bfv -k 9  --input bfv.in prove`

This creates a SNARK proof, stored as a binary file `data/bfv.snark`, using the inputs read (by default) from `data/halbfvo2_lib.in``. You can specify a different input file with the option `--input filename.in`, which would look for a file at `data/filename.in``.

Using the same proving key, you can generate proofs for the same ZK circuit on different inputs using this command.

**Proof Verification**

`cargo run --example bfv -- --name bfv -k 9 --input bfv.in verify`

Verify the proof generated above

## Chips 

- `check_poly_coefficients_in_range` - Enforces polynomial coefficients to be within a specified range
- `check_poly_from_distribution_chi_key` - Enforces polynomial to be sampled from the chi key
- `poly_add` - Enforces polynomial addition
- `constrain_poly_mul` - Constrain polynomial multiplication
- `poly_scalar_mul` - Enforces scalar multiplication of a polynomial
- `poly_reduce_by_modulo_q` - Enforces reduction of polynomial coefficients by a modulus
- `constraint_poly_reduction_by_cyclo` - Constrain polynomial reduction by a cyclotomic polynomial

## Benchmarks

- **Proving time: 9.5s** 
- **Verification time: 312ms**

Benches run using `bfv_2` run on M2 Macbook Pro with 12 cores and 32GB of RAM.

DEG and Q Parameters of the BFV encryption scheme should be chosen according to TABLES of RECOMMENDED PARAMETERS for 128-bits security level => https://homomorphicencryption.org/wp-content/uploads/2018/11/HomomorphicEncryptionStandardv1.1.pdf. 

In order to reproduce the benchmark modify the following parameters in `examples/bfv.rs`:

```rust
const DEG: usize = 1024;
const Q: u64 = 536870909;
```

Then run the following commands:

```bash
cargo run --example bfv -- --name bfv -k 13 --input bfv_2.in keygen
cargo run --example bfv -- --name bfv -k 13 --input bfv_2.in prove
cargo run --example bfv -- --name bfv -k 13 --input bfv_2.in verify
```