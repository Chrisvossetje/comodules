use crate::linalg::field::Field;
use std::fmt::Debug;

use super::field::CRing;

pub trait RModMorphism<R: CRing> {
    fn zero(domain: usize, codomain: usize) -> Self;
    fn identity(d: usize) -> Self;

    fn get(&self, domain: usize, codomain: usize) -> R;
    fn set(&mut self, domain: usize, codomain: usize, f: R);
    fn add_at(&mut self, domain: usize, codomain: usize, f: R);

    fn get_row(&self, codomain: usize) -> &[R];
    fn set_row(&mut self, codomain: usize, row: &[R]);

    fn is_row_non_zero(&self, codomain: usize) -> bool {
        (0..self.domain()).any(|domain| !self.get(domain, codomain).is_zero())
    }

    // domain l == codomain r, l \circ r
    fn compose(&self, rhs: &Self) -> Self;

    fn transpose(&self) -> Self;

    fn domain(&self) -> usize;
    fn codomain(&self) -> usize;

    fn vstack(&mut self, other: &mut Self);
    fn block_sum(&mut self, other: &Self);
}

pub trait Matrix<F: Field>: RModMorphism<F> + Clone + Send + Sync + PartialEq + Debug {
    // OR JUST A USIZE, NOT KNOWING WHICH COLUMN IT ORIGINATED FROM ?
    // nah, both
    fn pivots(&self) -> Vec<(usize, usize)>;

    // Returns the change of basis matrix for the domain!
    fn rref(&mut self);

    fn kernel(&self) -> Self;
    fn cokernel(&self) -> Self {
        self.transpose().kernel()
    }

    fn first_non_zero_entry(&self) -> Option<(usize, usize)>;
}
