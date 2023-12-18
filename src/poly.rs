use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::Zero;

/// Struct to perform polynomial operations (outside of the circuit).
#[derive(Clone, Debug)]
pub struct Poly {
    pub coefficients: Vec<BigInt>,
    pub degree: usize,
}

impl Poly {
    /// Create a new polynomial from a vector of strings.
    /// Coefficients are parsed as [a_deg, a_deg-1, ..., a_1, a_0] where a_0 is the constant term
    pub fn from_string(coefficients: Vec<String>) -> Self {
        let coefficients: Vec<BigInt> = coefficients
            .iter()
            .map(|x| BigInt::parse_bytes(x.as_bytes(), 10).unwrap())
            .collect();

        let degree = coefficients.len() - 1;

        Self {
            coefficients,
            degree,
        }
    }

    /// Create a new polynomial from a vector of big ints.
    /// Coefficients are parsed as [a_deg, a_deg-1, ..., a_1, a_0] where a_0 is the constant term
    pub fn from_big_int(coefficients: Vec<BigInt>) -> Self {
        let degree = coefficients.len() - 1;

        Self {
            coefficients,
            degree,
        }
    }

    /// Get the degree.
    pub fn deg(&self) -> usize {
        self.degree
    }

    /// Reduce coefficients by modulus
    pub fn reduce_by_modulus(&mut self, modulus: u64) {
        for coeff in self.coefficients.iter_mut() {
            *coeff = coeff.mod_floor(&BigInt::from(modulus));
        }
    }

    /// Perform polynomial multiplication.
    pub fn mul(&self, other: &Self) -> Self {
        let deg_a = self.deg();
        let deg_b = other.deg();
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

        Self::from_big_int(c)
    }

    /// Perform polynomial division between self and other, where self is the dividend and other is the divisor.
    /// Returns the quotient and remainder.
    /// The quotient is padded with zeroes at the beginning to make its degree equal to other.deg()
    /// The remainder is padded with zeroes at the beginning to make its degree equal to 2 * other.deg()
    pub fn div_euclid(&self, other: &Poly) -> (Self, Self) {
        if other.coefficients.is_empty() || other.coefficients.iter().all(BigInt::is_zero) {
            panic!("Cannot divide by a zero polynomial!");
        }
        if self.coefficients.is_empty() || self.coefficients.iter().all(BigInt::is_zero) {
            return (
                Self::from_big_int(vec![BigInt::zero(); self.deg() - other.deg() + 1]),
                Self::from_big_int(vec![BigInt::zero(); other.deg()]),
            );
        }

        assert_eq!(self.coefficients.len() - 1, self.deg());
        assert_eq!(other.coefficients.len() - 1, other.deg());

        let mut dividend = self.coefficients.to_vec();
        let divisor = other.coefficients.to_vec();
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

        // pad quotient with zeroes at the beginning to make its degree equal to self.deg() - other.deg()
        let mut quotient = quotient;
        while quotient.len() - 1 < other.deg() {
            quotient.insert(0, BigInt::from(0u32));
        }

        // pad remainder with zeroes at the beginning to make its degree equal to 2 * other.deg() - 2
        let mut remainder = remainder;
        while remainder.len() - 1 < 2 * other.deg() {
            remainder.insert(0, BigInt::from(0u32));
        }

        (Self::from_big_int(quotient), Self::from_big_int(remainder))
    }
}
