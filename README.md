# zk-fhe
Zk proving the correct execution of encryption operation under BFV Fully Homomorphic Encryption scheme

### TO DOs

// - [ ] Add assumptions for each chip
// - [ ] Add master example to test the whole encryption/decryption flow
// - [ ] Avoid overflow in poly addition and poly multiplication in the circuit level 
// - [ ] Add function to reduce by cyclotomic polynomial in poly_utils
// - [ ] Pass inputs as variable instead of JSON

`cargo run --example bfv -- --name poly_input -k 11  mock`