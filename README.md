# zk-fhe
Zk proving the correct execution of encryption operation under BFV Fully Homomorphic Encryption scheme

Implementation based on [Revisiting Homomorphic Encryption Schemes for Finite Fields](https://eprint.iacr.org/2021/204.pdf)

The application is not production ready and is only meant to be used for educational purposes.

`cargo run --example bfv -- --name bfv -k 9  mock`

The input data is located in the `data` folder. This test vector file can be generated using [rlwe-py](https://github.com/yuriko627/rlwe-py)

### Chips 

- `check_poly_from_distribution_chi_error` - Enforces polynomial to be sampled from the chi distribution
- `check_poly_from_distribution_chi_key` - Enforces polynomial to be sampled from the chi key
- `poly_add` - Enforces polynomial addition
- `poly_mul_equal_deg` - Enforces polynomial multiplication between polynomials of equal degree
- `poly_mul_diff_deg` - Enforces polynomial multiplication between polynomials of different degree
- `poly_scalar_mul` - Enforces scalar multiplication of a polynomial
- `poly_reduce` - Enforces reduction of polynomial coefficients by a modulus
- `poly_divide_by_cyclo` - Enforces the reduction of a polynomial by a cyclotomic polynomial
