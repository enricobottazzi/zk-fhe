# zk-fhe
Zk proving the correct execution of encryption operation under BFV Fully Homomorphic Encryption scheme

Implementation based on [Revisiting Homomorphic Encryption Schemes for Finite Fields](https://eprint.iacr.org/2021/204.pdf)

The application is not production ready and is only meant to be used for educational purposes.

`LOOKUP_BITS=8 cargo run --example bfv -- --name bfv -k 14  mock`

The input data is located in the `data` folder. This file can be generated using [rlwe-py](https://github.com/yuriko627/rlwe-py)