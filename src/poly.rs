use halo2_base::utils::log2_ceil;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::Zero;

/// Struct to perform polynomial operations (outside of the circuit).
/// `max_bits` keeps track of the maximum number of bits that a coefficient can have.
#[derive(Clone, Debug)]
pub struct Poly {
    pub coefficients: Vec<BigInt>,
    pub degree: usize,
    pub max_bits: u64,
}

impl Poly {
    /// Create a new polynomial from a vector of strings. `modulus` is the modulus in which the coefficients should live.
    /// Coefficients must be passed as `[a_deg, a_deg-1, ..., a_1, a_0]` where `a_0` is the constant term
    ///
    /// # Assumptions
    /// * All coefficients are less than or equal to `modulus`
    pub fn from_string(coefficients: Vec<String>, modulus: u64) -> Self {
        let modulus_bigint = BigInt::from(modulus);

        let coefficients: Vec<BigInt> = coefficients
            .iter()
            .map(|x| {
                let coeff = BigInt::parse_bytes(x.as_bytes(), 10).unwrap();
                assert!(coeff <= modulus_bigint);
                coeff
            })
            .collect();

        let degree = coefficients.len() - 1;

        Self {
            coefficients,
            degree,
            max_bits: modulus_bigint.bits(),
        }
    }

    /// Create a new polynomial from a vector of big ints.
    /// Coefficients must be passed as [a_deg, a_deg-1, ..., a_1, a_0] where a_0 is the constant term
    ///
    /// # Assumptions
    /// * All coefficients are at most `max_bits` bits long
    fn from_big_int(coefficients: Vec<BigInt>, max_bits: u64) -> Self {
        let degree = coefficients.len() - 1;

        for coeff in &coefficients {
            assert!(coeff.bits() <= max_bits);
        }

        Self {
            coefficients,
            degree,
            max_bits,
        }
    }

    /// Get the degree.
    pub fn deg(&self) -> usize {
        self.degree
    }

    /// Get the max number of bits.
    fn max_bits(&self) -> u64 {
        self.max_bits
    }

    /// Perform polynomial multiplication between `self` and `other` and return the result.
    ///
    /// # Assumptions
    /// * self and other are the same degree
    pub fn mul(&self, other: &Self) -> Self {
        let deg_a = self.deg();
        let deg_b = other.deg();
        assert_eq!(deg_a, deg_b);

        let deg_c = deg_a + deg_b;

        // initialize the output polynomial with zeroes
        let mut c = vec![BigInt::from(0_u64); deg_c + 1];

        // perform polynomial multiplication
        for i in 0..=deg_a {
            for j in 0..=deg_b {
                c[i + j] += self.coefficients[i].clone() * other.coefficients[j].clone();
            }
        }

        assert!(c.len() == deg_c + 1);

        // Each coefficient in the resulting polynomial is a sum of products of coefficients from the input polynomials.
        // Thus, the maximum number of bits in any coefficient of the result is influenced by:
        // * Number of bits in a Single Product: this is the sum of the max_bits of the two input polynomials = a.max_bits() + b.max_bits()
        // * Number of terms contributing to a single coefficient:
        // In the worst case scenario, the number of terms contributing to a single coefficient is deg_a (=deg_b) + 1.
        // Considering the sum of deg_a + 1 elements of a.max_bits() + b.max_bits() bits each, the maximum number of bits in a single coefficient can be defined as follows:
        // max_bits = deg_a.max_bits() + deg_b.max_bits() + ceil(log2(deg_a + 1))
        let max_bits = self.max_bits() + other.max_bits() + log2_ceil(deg_a as u64 + 1) as u64;
        Self::from_big_int(c, max_bits)
    }

    /// Perform polynomial division between `self` and cyclotomic polynomial `cyclo`, where self is the dividend and `cyclo` is the divisor.
    /// Returns `quotient` and `remainder`.
    /// * `quotient` is padded with zeroes at the beginning to make its degree equal to `cyclo.deg()`
    /// * `remainder` is padded with zeroes at the beginning to make its degree equal to `2 * cyclo.deg()`
    /// * `remainder` is reduced by `modulus`
    ///
    /// # Assumptions
    /// * `cyclo` is a cyclotomic polynomial of form `x^n + 1`
    pub fn divide_by_cyclo(&self, cyclo: &Poly, modulus: u64) -> (Self, Self) {
        // In division by a cyclotomic polynomial, the quotient polynomial have at most modulus.bits() bits. => https://hackmd.io/EcLZAJiyQRShalL3BLZJ7Q?view#Polynomial-Division
        // After reduction by modulus, the remainder polynomial have at most modulus.bits() bits.
        let modulus_bits = BigInt::from(modulus).bits();

        if self.coefficients.is_empty() || self.coefficients.iter().all(BigInt::is_zero) {
            return (
                Self::from_big_int(vec![BigInt::zero(); cyclo.deg() + 1], modulus_bits), // We set `max_bits` to `modulus_bits` for correct `keygen` phase
                Self::from_big_int(vec![BigInt::zero(); 2 * cyclo.deg() + 1], modulus_bits), // We set `max_bits` to `modulus_bits` for correct `keygen` phase
            );
        }

        assert_eq!(self.coefficients.len() - 1, self.deg());
        assert_eq!(cyclo.coefficients.len() - 1, cyclo.deg());

        let mut dividend = self.coefficients.to_vec();
        let divisor = cyclo.coefficients.to_vec();
        let mut quotient: Vec<BigInt> = Vec::new();
        let mut remainder: Vec<BigInt> = Vec::new();

        while dividend.len() > divisor.len() - 1 {
            let leading_coefficient_ratio = &dividend[0] / &divisor[0];
            quotient.push(leading_coefficient_ratio.clone());

            for (i, coeff) in divisor.iter().enumerate() {
                dividend[i] = &dividend[i] - &(&leading_coefficient_ratio * coeff);
            }

            dividend.remove(0);
        }

        for coeff in &dividend {
            remainder.push(coeff.clone());
        }

        while !quotient.is_empty() && quotient[0].is_zero() {
            quotient.remove(0);
        }

        while !remainder.is_empty() && remainder[0].is_zero() {
            remainder.remove(0);
        }

        // pad quotient with zeroes at the beginning to make its degree equal to self.deg() - cyclo.deg()
        let mut quotient = quotient;
        while quotient.len() - 1 < cyclo.deg() {
            quotient.insert(0, BigInt::from(0u32));
        }

        // pad remainder with zeroes at the beginning to make its degree equal to 2 * cyclo.deg() - 2
        let mut remainder = remainder;
        while remainder.len() - 1 < 2 * cyclo.deg() {
            remainder.insert(0, BigInt::from(0u32));
        }

        // reduce remainder by modulus
        remainder = remainder
            .iter()
            .map(|x| x.mod_floor(&BigInt::from(modulus)))
            .collect();
        (
            Self::from_big_int(quotient, modulus_bits),
            Self::from_big_int(remainder, modulus_bits),
        )
    }

    /// Reduce coefficients by `modulus` and return the result.
    pub fn reduce_by_modulus(&mut self, modulus: u64) -> Poly {
        let modulus_bigint = BigInt::from(modulus);
        let coefficients = self
            .coefficients
            .iter()
            .map(|x| x.mod_floor(&modulus_bigint))
            .collect();

        let modulus_bits = BigInt::from(modulus).bits();

        Self::from_big_int(coefficients, modulus_bits)
    }
}
