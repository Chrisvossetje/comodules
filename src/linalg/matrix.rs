use crate::linalg::field::Field;

pub trait Matrix<F: Field>: Clone + Send + Sync + PartialEq {
    fn zero(domain: usize, codomain: usize) -> Self;
    fn identity(d: usize) -> Self;

    fn get(&self, domain: usize, codomain: usize) -> F;
    fn set(&mut self, domain: usize, codomain: usize, f: F);
    fn add_at(&mut self, domain: usize, codomain: usize, f: F);

    fn get_row(&self, codomain: usize) -> &[F];
    fn set_row(&mut self, codomain: usize, row: &[F]);

    // domain l == codomain r, l \circ r
    fn compose(&self, rhs: &Self) -> Self;

    fn transpose(&self) -> Self;

    fn domain(&self) -> usize;
    fn codomain(&self) -> usize;

    // OR JUST A USIZE, NOT KNOWING WHICH COLUMN IT ORIGINATED FROM ?
    // nah, both
    fn pivots(&self) -> Vec<(usize, usize)>;

    fn vstack(&self, other: &Self) -> Self;
    fn block_sum(&self, other: &Self) -> Self;

    // Returns the change of basis matrix for the domain!
    fn rref(&mut self);

    fn kernel(&self) -> Self;
    fn cokernel(&self) -> Self {
        self.transpose().kernel()
    }

    fn first_non_zero_entry(&self) -> Option<(usize, usize)>;
}
