use halo2_base::{utils::ScalarField, AssignedValue};

/// Performs long polynomial division on two polynomials
/// Returns the quotient and remainder
/// 
/// * Input polynomials are parsed as a vector of assigned coefficients [a_DEG, a_DEG-1, ..., a_1, a_0] where a_0 is the constant term
/// * DEG_DVD is the degree of the dividend
/// * DEG_DVS is the degree of the divisor
/// * Q is the modulus of the Ring. All the coefficients will be in the range [0, Q-1]
/// * Assumes that coefficients of the dividend and divisor are u64 values
pub fn div_euclid<const DEG_DVD: usize, const DEG_DVS: usize, const Q: u64>(
    dividend: &Vec<u64>,
    divisor: &Vec<u64>,
) -> (Vec<u64>, Vec<u64>) {
    if divisor.is_empty() || divisor.iter().all(|&x| x == 0) {
        panic!("Cannot divide by a zero polynomial!");
    }
    if dividend.is_empty() || dividend.iter().all(|&x| x == 0) {
        let quotient = vec![0; DEG_DVD - DEG_DVS + 1];
        let remainder = vec![0; DEG_DVS];

        // turn quotient and remainder into u64
        let quotient = quotient.iter().map(|&x| x as u64).collect::<Vec<u64>>();
        let remainder = remainder.iter().map(|&x| x as u64).collect::<Vec<u64>>();

        return (quotient, remainder);
    }

    // assert that the degree of the dividend is equal to DEG_DVD
    assert_eq!(dividend.len() - 1, DEG_DVD);

    // assert that the degree of the divisor is equal to DEG_DVS
    assert_eq!(divisor.len() - 1, DEG_DVS);

    // transform the dividend and divisor into a vector of i64
    let mut dividend = dividend.iter().map(|&x| x as i64).collect::<Vec<i64>>();
    let divisor = divisor.iter().map(|&x| x as i64).collect::<Vec<i64>>();

    let mut quotient = Vec::new();
    let mut remainder = Vec::new();

    while dividend.len() > divisor.len() - 1 {
        let leading_coefficient_ratio = dividend[0] / divisor[0];
        quotient.push(leading_coefficient_ratio);

        for (i, coeff) in divisor.iter().enumerate() {
            let diff = dividend[i] - leading_coefficient_ratio * *coeff;
            dividend[i] = diff;
        }

        dividend.remove(0);
    }

    for coeff in &dividend {
        remainder.push(*coeff);
    }

    // Trim the leading zeros from quotient and remainder
    while !quotient.is_empty() && quotient[0] == 0 {
        quotient.remove(0);
    }

    while !remainder.is_empty() && remainder[0] == 0 {
        remainder.remove(0);
    }

    // Range over remainder. If any element is negative, add Q to it
    for coeff in &mut remainder {
        if *coeff < 0 {
            *coeff += Q as i64;
        }
    }

    // Convert remainder back to u64
    let remainder = remainder.iter().map(|&x| x as u64).collect::<Vec<u64>>();

    // Convert quotient back to u64
    let quotient = quotient.iter().map(|&x| x as u64).collect::<Vec<u64>>();

    // pad quotient with zeroes at the beginning to make its degree equal to DEG_DVD - DEG_DVS
    let mut quotient = quotient;
    while quotient.len() - 1 < DEG_DVD - DEG_DVS {
        quotient.insert(0, 0);
    }

    (quotient, remainder)
}

/// Convert a vector of AssignedValue to a vector of u64
/// 
/// * Assumes that each element of AssignedValue can be represented in 8 bytes
pub fn vec_assigned_to_vec_u64<F: ScalarField>(vec: &Vec<AssignedValue<F>>) -> Vec<u64> {
    let mut vec_u64 = Vec::new();

    for i in 0..vec.len() {
        let value_bytes_le = vec[i].value().to_bytes_le();

        // slice value_to_bytes_le the first 8 bytes
        let value_8_bytes_le = &value_bytes_le[..8];
        let mut array_value_8_bytes_le = [0u8; 8];
        array_value_8_bytes_le.copy_from_slice(value_8_bytes_le);
        let num = u64::from_le_bytes(array_value_8_bytes_le);
        vec_u64.push(num);
    }
    vec_u64
}
