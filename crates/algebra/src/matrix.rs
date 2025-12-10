use deepsize::DeepSizeOf;

use crate::ring::CRing;
use std::fmt::Debug;


/// This represents a map between (free) R-modules
/// (or represents vectors inside an R-module)
pub trait Matrix<R: CRing>: Clone + Send + Sync + PartialEq + Debug + DeepSizeOf {
    fn zero(domain: usize, codomain: usize) -> Self;
    fn identity(d: usize) -> Self;

    fn get(&self, domain: usize, codomain: usize) -> R;
    fn set(&mut self, domain: usize, codomain: usize, r: R);
    fn add_at(&mut self, domain: usize, codomain: usize, r: R);

    fn get_row(&self, codomain: usize) -> &[R];
    fn set_row(&mut self, codomain: usize, row: &[R]);
    
    fn scalar_multiply_row(&mut self, codomain: usize, r: R) {
        for domain in 0..self.domain() {
            let el = self.get(domain, codomain);
            self.set(domain, codomain, r * el);
        }
    }

    fn eval_vector(&self, vector: &[R]) -> Vec<R> {
        debug_assert_eq!(self.domain(), vector.len());
        (0..self.codomain()).map(|y| {
            let r_1 = self.get_row(y);
            let el: R = r_1.iter().zip(vector).map(|(x, y)| {
                *x * *y
            }).sum();
            el
        }).collect()
    }

    fn scalar_multiply_column(&mut self, domain: usize, r: R) {
        for codomain in 0..self.codomain() {
            let el = self.get(domain, codomain);
            self.set(domain, codomain, r * el);
        }
    }
    
    fn get_column(&self, domain: usize) -> Vec<R> {
        let mut r = vec![];
        for i in 0..self.codomain() {
            let el = self.get(domain, i);
            r.push(el);
        }
        r
    }


    fn set_column(&mut self, domain: usize, row: &[R]) {
        for i in 0..self.codomain() {
            self.set(domain, i, row[i]);
        }
    }

    fn set_row_zero(&mut self, codomain: usize) {
        for d in 0..self.domain() {
            self.set(d, codomain, R::zero());
        }
    }

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

    fn extend_one_row(&mut self); 

    /// Swap two rows
    fn swap_rows(&mut self, row1: usize, row2: usize) {
        if row1 == row2 {
            return;
        }
        for j in 0..self.domain() {
            let temp = self.get(j, row1);
            self.set(j, row1, self.get(j, row2));
            self.set(j, row2, temp);
        }
    }

    fn is_unit(&self) -> Result<(),()> {
        for i in 0..self.domain() {
            for j in 0..self.codomain() {
                if i == j {
                    if self.get(i, j) != R::one(){
                        return Err(())
                    }
                } else {
                    if !self.get(i, j).is_zero(){
                        return Err(())
                    }

                }
            }
        }
        Ok(())
    }

    /// Swap two columns  
    fn swap_cols(&mut self, col1: usize, col2: usize) {
        if col1 == col2 {
            return;
        }
        for i in 0..self.codomain() {
            let temp = self.get(col1, i);
            self.set(col1, i, self.get(col2, i));
            self.set(col2, i, temp);
        }
    }

    /// Add a multiple of one row to another: row[target] += factor * row[source]
    fn add_row_multiple(&mut self, target: usize, source: usize, factor: R) {
        if factor.is_zero() {
            return;
        }
        for j in 0..self.domain() {
            let addition = factor.clone() * self.get(j, source);
            self.add_at(j, target, addition);
        }
    }

    /// Add a multiple of one column to another: col[target] += factor * col[source]
    fn add_col_multiple(&mut self, target: usize, source: usize, factor: R) {
        if factor.is_zero() {
            return;
        }
        for i in 0..self.codomain() {
            let addition = factor.clone() * self.get(source, i);
            self.add_at(target, i, addition);
        }
    }
}

