# zk-fhe
Zk proving the correct execution of encryption operation under BFV Fully Homomorphic Encryption scheme

Implementation based on [Revisiting Homomorphic Encryption Schemes for Finite Fields](https://eprint.iacr.org/2021/204.pdf)

## Disclaimer

This is a research project and is not meant to be used in production. The code is not audited.

## Quick Start

**Mock Prover**

```bash
cargo run --example bfv -- --name bfv -k 9 --input bfv_test/bfv.in mock
```

The `MockProver` does not run the cryptographic prover on your circuit, but instead directly checks if constraints are satisfied. This is useful for testing purposes, and runs faster than the actual prover.

- `bfv` is the name of the circuit located in `examples/bfv.rs` 
- `bfv_test/bfv.in` is the input file for the circuit located in `data/bfv_test/bfv.in`. A different test vector file can be generated using [bfv-py](https://github.com/yuriko627/bfv-py)
- `-k` is the DEGREE of the circuit as you specify to set the circuit to have `2^k` number of rows. The number of rows is determined by the number of constraints in the circuit. Working with larger data inputs will require a larger degree.

**Key Generation**

```bash
cargo run --example bfv -- --name bfv -k 9 --input bfv_test/bfv_empty.in keygen
```

To generate a random universal trusted setup (for testing only!) and the proving and verifying keys for your circuit.

For technical reasons (related to [axiom Halo2-scaffold](https://github.com/axiom-crypto/halo2-scaffold)), keygen still requires an input file of the correct format. In this case, the input file is empty as the the actual input data are not encoded in the key generation. 

This will generate a proving key `data/bfv.pk` and a verifying key `data/bfv.vk`. It will also generate a file `configs/bfv.json` which describes (and pins down) the configuration of the circuit. This configuration file is later read by the prover.

**Proof Generation**

```
cargo run --example bfv -- --name bfv -k 9 --input bfv_test/bfv.in prove
```

Note: during proof generation we must pass an input file containing the actual input data. 

**Proof Verification**

```
cargo run --example bfv -- --name bfv -k 9 --input bfv_test/bfv_empty.in verify
```

Note: during proof verification we can pass an empty input file.

## Benchmarks

- **Proving time: 8.5s** 
- **Verification time: 299ms**

Benches run using `bfv_2` run on M2 Macbook Pro with 12 cores and 32GB of RAM.

N and Q Parameters of the BFV encryption scheme should be chosen according to TABLES of RECOMMENDED PARAMETERS for 128-bits security level => https://homomorphicencryption.org/wp-content/uploads/2018/11/HomomorphicEncryptionStandardv1.1.pdf. 

In order to reproduce the benchmark modify the following parameters in `examples/bfv.rs`:

```rust
const N: usize = 1024;
const Q: u64 = 536870909;
```

Then run the following commands:

```bash
cargo run --example bfv -- --name bfv -k 13 --input bfv/bfv.in mock
cargo run --example bfv -- --name bfv -k 13 --input bfv/bfv_empty.in keygen
cargo run --example bfv -- --name bfv -k 13 --input bfv/bfv.in prove
cargo run --example bfv -- --name bfv -k 13 --input bfv/bfv_empty.in verify
```

## Warning: Overflow 

Many polynomial operations performed inside the [circuit](./examples/bfv.rs) involve careful handling of coefficients in order to avoid overflows on the prime field. This [guide](https://zipcpu.com/dsp/2017/07/21/bit-growth.html) is recommended to understand the bit growth of the coefficients when performing polynomial operations. `N` and `DEG` must be provided at `keygen` time. Certain combinations of `N` and `DEG` can potentially lead to the risk of overflow during proof generation, which is something that can be maliciously exploited by the prover. `keygen` will fail if this is these cases.