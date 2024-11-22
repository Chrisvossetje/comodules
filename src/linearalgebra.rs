use crate::fp::PrimeField;

pub trait Matrix<F: PrimeField> {
    fn get_matrix(&self) -> Vec<Vec<F>>;
}


