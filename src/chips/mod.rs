pub mod poly_big_int_chip;
pub mod poly_distribution;
pub mod poly_operations;
pub mod utils;

use axiom_eth::Field;
use halo2_base::AssignedValue;

/// A struct that holds a polynomial and its length.
#[derive(Clone, Debug)]
pub struct PolyWithLength<F: Field> {
    assigned_poly: Vec<AssignedValue<F>>,
    assigned_length: AssignedValue<F>,
}

impl<F: Field> PolyWithLength<F> {
    /// Create a new polynomial with length.
    pub fn new(assigned_poly: Vec<AssignedValue<F>>, assigned_length: AssignedValue<F>) -> Self {
        assert!(assigned_length.value() == &F::from(assigned_poly.len() as u64));
        Self {
            assigned_poly,
            assigned_length,
        }
    }

    /// Get the polynomial.
    pub fn get_poly(&self) -> &Vec<AssignedValue<F>> {
        &self.assigned_poly
    }

    /// Get the length.
    pub fn get_length(&self) -> &AssignedValue<F> {
        &self.assigned_length
    }
}
