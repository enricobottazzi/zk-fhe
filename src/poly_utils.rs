use ark_bn254::Fr;
use ark_ff::fields::PrimeField;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
use halo2_base::utils::ScalarField;

// Process polynomial coefficients to fit in the [0, q) range
// q is the modulus of the cipher text field
// If a coefficient is negative, take (q - |c|)
// If a coefficient is positive and greater than q, take c % q
fn process_coefficients(coefficients: Vec<i64>, q: u64) -> Vec<u64> {
    let mut processed_coefficients: Vec<u64> = Vec::new();

    for c in coefficients {
        if c < 0 {
            processed_coefficients.push(q - (c.unsigned_abs()));
        } else {
            processed_coefficients.push((c as u64) % q);
        }
    }

    for c in &processed_coefficients {
        assert!(*c < q);
    }

    processed_coefficients
}

fn reduce_poly_by_modulus(poly: DensePolynomial<Fr>, q: u64) -> DensePolynomial<Fr> {
    let unreduced_poly_coefficients = poly
        .coeffs
        .iter()
        .map(|x| x.into_bigint().to_string())
        .collect::<Vec<String>>();

    // convert to u64
    let unreduced_poly_coefficients = unreduced_poly_coefficients
        .iter()
        .map(|x| x.parse::<u64>().unwrap())
        .collect::<Vec<u64>>();

    // reduce coefficients by looping through them
    let mut reduced_poly_coefficients: Vec<u64> = Vec::new();

    for c in unreduced_poly_coefficients {
        reduced_poly_coefficients.push(c % q);
    }

    DensePolynomial::<Fr>::from_coefficients_vec(
        reduced_poly_coefficients
            .iter()
            .map(|x| Fr::from(*x))
            .collect::<Vec<Fr>>(),
    )
}

pub fn poly_to_scalar_field_coeffs<F: ScalarField>(poly: DensePolynomial<Fr>) -> Vec<F> {
    poly.coeffs
        .iter()
        .map(|x| {
            let coeff_string = x.into_bigint().to_string();
            F::from_str_vartime(&coeff_string).unwrap()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_process_coefficients() {
        let coefficients = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let q = 11;
        let processed_coefficients = process_coefficients(coefficients, q);

        assert_eq!(processed_coefficients, vec![1, 2, 3, 4, 5, 6, 7, 8]);

        let coefficients = vec![-1, -2, -3, -4, -5, -6, -7, -8];
        let q = 11;
        let processed_coefficients = process_coefficients(coefficients, q);

        assert_eq!(processed_coefficients, vec![10, 9, 8, 7, 6, 5, 4, 3]);

        let coefficients = vec![1, 2, 3, 4, 5, 6, 11, 23];

        let q = 11;
        let processed_coefficients = process_coefficients(coefficients, q);

        assert_eq!(processed_coefficients, vec![1, 2, 3, 4, 5, 6, 0, 1]);
    }

    #[test]
    fn test_reduce_poly_by_modulus() {
        let poly = DensePolynomial::<Fr>::from_coefficients_vec(vec![
            Fr::from(1u64),
            Fr::from(2u64),
            Fr::from(3u64),
            Fr::from(4u64),
            Fr::from(5u64),
            Fr::from(13u64),
            Fr::from(12u64),
            Fr::from(11u64),
        ]);

        let q = 11;

        let reduced_poly = reduce_poly_by_modulus(poly, q);

        assert_eq!(
            reduced_poly,
            DensePolynomial::<Fr>::from_coefficients_vec(vec![
                Fr::from(1u64),
                Fr::from(2u64),
                Fr::from(3u64),
                Fr::from(4u64),
                Fr::from(5u64),
                Fr::from(2u64),
                Fr::from(1u64),
                Fr::from(0u64),
            ])
        );
    }
}
