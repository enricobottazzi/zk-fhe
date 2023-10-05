


pub fn div_euclid(dividend: &Vec<u64>, divisor: &Vec<u64>) -> (Vec<u64>, Vec<u64>) {
    if divisor.is_empty() || divisor.iter().all(|&x| x == 0) {
        panic!("Cannot divide by a zero polynomial!");
    }

    // transform the dividend and divisor into a vector of i64
    let mut dividend = dividend.iter().map(|&x| x as i64).collect::<Vec<i64>>();
    let mut divisor = divisor.iter().map(|&x| x as i64).collect::<Vec<i64>>();

    // reverse the dividend and divisor
    dividend.reverse();
    divisor.reverse();

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
    let mut remainder = remainder.iter().map(|&x| x as u64).collect::<Vec<u64>>();

    // Convert quotient back to u64
    let mut quotient = quotient.iter().map(|&x| x as u64).collect::<Vec<u64>>();

    // reverse the quotient and remainder
    quotient.reverse();
    remainder.reverse();

    (quotient, remainder)
}