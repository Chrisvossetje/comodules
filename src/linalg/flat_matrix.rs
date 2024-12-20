use serde::{Deserialize, Serialize};

use super::{field::Field, matrix::Matrix};

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct FlatMatrix<F: Field> {
    pub data: Vec<F>,
    pub domain: usize,
    pub codomain: usize,
}

impl<F: Field> Matrix<F> for FlatMatrix<F> {
    fn zero(domain: usize, codomain: usize) -> Self {
        todo!()
    }

    fn identity(d: usize) -> Self {
        todo!()
    }

    fn get(&self, domain: usize, codomain: usize) -> F {
        todo!()
    }

    fn set(&mut self, domain: usize, codomain: usize, f: F) {
        todo!()
    }

    fn add_at(&mut self, domain: usize, codomain: usize, f: F) {
        todo!()
    }

    fn get_row(&self, codomain: usize) -> &[F] {
        todo!()
    }

    fn set_row(&mut self, codomain: usize, row: &[F]) {
        todo!()
    }

    fn compose(&self, rhs: &Self) -> Self {
        todo!()
    }

    fn transpose(&self) -> Self {
        todo!()
    }

    fn domain(&self) -> usize {
        todo!()
    }

    fn codomain(&self) -> usize {
        todo!()
    }

    fn pivots(&self) -> Vec<(usize, usize)> {
        todo!()
    }

    fn vstack(&mut self, other: &mut Self) {
        todo!()
    }

    fn block_sum(&mut self, other: &mut Self) {
        todo!()
    }

    fn rref(&mut self) {
        todo!()
    }

    fn kernel(&self) -> Self {
        todo!()
    }

    fn first_non_zero_entry(&self) -> Option<(usize, usize)> {
        todo!()
    }
}
