use core::panic;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

use crate::{linalg::ring::CRing};
use super::{
    field::Field,
    matrix::{Matrix, RModMorphism},
};

#[derive(Clone, PartialEq, Deserialize, Serialize)]
pub struct FlatMatrix<R: CRing> {
    pub data: Vec<R>,
    pub domain: usize,
    pub codomain: usize,
}


impl<R: CRing> Debug for FlatMatrix<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for m in 0..self.codomain {
            for n in 0..self.domain {
                let el = self.get(n, m);
                let s = format!("{:?}", el);
                write!(f, "{:<4}", s)?;
            }
            writeln!(f)?;
        }
        writeln!(f)
    }
}

impl<R: CRing> FlatMatrix<R> {
    fn get_element(&self, row: usize, col: usize) -> R {
        self.data[row * self.domain + col].clone()
    }

    fn set_element(&mut self, row: usize, col: usize, value: R) {
        self.data[row * self.domain + col] = value;
    }

    pub fn extensive_pivots(&self) -> Vec<usize> {
        let mut v = vec![None; self.codomain];
        for codom_id in 0..self.domain {
            let column = self.get_column(codom_id);
            let mut units = 0;
            let mut final_id = usize::MAX;
            for (id,a) in column.iter().enumerate() {
                if !a.is_zero() {
                    if a.is_unit() {
                        units += 1;
                        final_id = id;
                    } else {
                        units += 2;
                    }
                }
            }
            if units == 1 {
                v[final_id] = Some(codom_id);
            }
        }

        if cfg!(debug_assertions) {
            for a in &v {
                if a.is_none() {
                    println!("{:?}", self);
                    panic!("Matrix does not contain all pivots");
                }
            }
        }

        v.iter().map(|f|
            f.unwrap()
        ).collect()
    }
}

impl<F: Field> FlatMatrix<F> {
    fn rref_kernel(&self) -> Self {
        let mut free_vars = Vec::new();
        let pivot_cols: Vec<usize> = self.pivots().iter().map(|x| x.0).collect();

        for j in 0..self.domain {
            if !pivot_cols.contains(&j) {
                free_vars.push(j);
            }
        }

        let mut kernel = vec![F::zero(); free_vars.len() * self.domain];

        for (i, &free_var) in free_vars.iter().enumerate() {
            let row_idx = i * self.domain;
            kernel[row_idx + free_var] = F::one();

            for (row, &pivot_col) in pivot_cols.iter().enumerate() {
                kernel[row_idx + pivot_col] = -self.get_element(row, free_var);
            }
        }

        Self {
            data: kernel,
            domain: self.domain,
            codomain: free_vars.len(),
        }
    }
}

impl<F: Field> Matrix<F> for FlatMatrix<F> {
    fn kernel(&self) -> Self {
        let mut clone = self.clone();
        clone.rref();
        let mut kernel = clone.rref_kernel();
        kernel.rref();
        kernel
    }

    fn rref(&mut self) {
        let mut lead = 0;

        for r in 0..self.codomain {
            if lead >= self.domain {
                break;
            }

            let mut i = r;
            while self.get_element(i, lead).is_zero() {
                i += 1;
                if i == self.codomain {
                    i = r;
                    lead += 1;
                    if lead == self.domain {
                        return;
                    }
                }
            }

            for j in 0..self.domain {
                self.data.swap(r * self.domain + j, i * self.domain + j);
            }

            let pivot = self.get_element(r, lead);
            if !pivot.is_zero() {
                let pivot_inv = pivot.inv().expect("Pivot should be invertible");
                for j in 0..self.domain {
                    let idx = r * self.domain + j;
                    self.data[idx] *= pivot_inv;
                }
            }

            for i in 0..self.codomain {
                if i != r {
                    let factor = self.get_element(i, lead);
                    for j in 0..self.domain {
                        let idx = i * self.domain + j;
                        let el = self.get_element(r, j);
                        self.data[idx] -= factor * el;
                    }
                }
            }

            lead += 1;
        }
    }

    fn pivots(&self) -> Vec<(usize, usize)> {
        let mut col = 0;
        let mut pivots = vec![];
        for row in 0..self.codomain {
            while col < self.domain {
                if !self.get_element(row, col).is_zero() {
                    pivots.push((col, row));
                    col += 1;
                    break;
                }
                col += 1;
            }
        }
        pivots
    }

    fn first_non_zero_entry(&self) -> Option<(usize, usize)> {
        for codom_id in 0..self.codomain {
            for dom_id in 0..self.domain {
                if !self.get_element(codom_id, dom_id).is_zero() {
                    return Some((codom_id, dom_id));
                }
            }
        }
        None
    }
}

impl<R: CRing> RModMorphism<R> for FlatMatrix<R> {
    fn transpose(&self) -> Self {
        let mut new_matrix = vec![R::zero(); self.data.len()];
        for i in 0..self.codomain {
            for j in 0..self.domain {
                new_matrix[j * self.codomain + i] = self.get_element(i, j);
            }
        }

        Self {
            data: new_matrix,
            domain: self.codomain,
            codomain: self.domain,
        }
    }

    fn vstack(&mut self, other: &mut Self) {
        debug_assert_eq!(
            self.domain(),
            other.domain(),
            "Domains of the two matrices do not have the same dimension"
        );

        self.data.extend_from_slice(&other.data);
        self.codomain += other.codomain;
    }

    fn block_sum(&mut self, other: &Self) {
        let new_domain = self.domain + other.domain;
        let new_codomain = self.codomain + other.codomain;
        let mut new = Self::zero(new_domain, new_codomain);

        for i in 0..self.codomain {
            let start = i * new_domain;
            new.data[start..(start + self.domain)].copy_from_slice(self.get_row(i));
        }

        for i in 0..other.codomain {
            let start = (self.codomain + i) * new_domain + self.domain;
            new.data[start..(start + other.domain)].copy_from_slice(other.get_row(i));
        }

        *self = new
    }

    fn identity(d: usize) -> Self {
        let mut data = vec![R::zero(); d * d];
        for i in 0..d {
            data[i * d + i] = R::one();
        }
        Self {
            data,
            domain: d,
            codomain: d,
        }
    }

    fn get(&self, domain: usize, codomain: usize) -> R {
        self.get_element(codomain, domain)
    }

    fn set(&mut self, domain: usize, codomain: usize, r: R) {
        self.set_element(codomain, domain, r);
    }

    fn add_at(&mut self, domain: usize, codomain: usize, r: R) {
        let idx = codomain * self.domain + domain;
        self.data[idx] += r;
    }

    fn get_row(&self, codomain: usize) -> &[R] {
        let start = codomain * self.domain;
        let end = start + self.domain;
        &self.data[start..end]
    }

    fn set_row(&mut self, codomain: usize, row: &[R]) {
        debug_assert!(row.len() <= self.domain);
        let start = codomain * self.domain;
        let end = start + row.len();
        self.data[start..end].copy_from_slice(row);
    }

    fn domain(&self) -> usize {
        self.domain
    }

    fn codomain(&self) -> usize {
        self.codomain
    }

    fn compose(&self, rhs: &Self) -> Self {
        debug_assert_eq!(
            self.domain, rhs.codomain,
            "Matrix domain not equal to codomain"
        );

        let mut compose = vec![R::zero(); self.codomain * rhs.domain];

        for x in 0..self.codomain {
            for y in 0..rhs.domain {
                let mut sum = R::zero();
                for k in 0..self.domain {
                    sum += self.get_element(x, k) * rhs.get_element(k, y);
                }
                compose[x * rhs.domain + y] = sum;
            }
        }

        Self {
            data: compose,
            domain: rhs.domain,
            codomain: self.codomain,
        }
    }

    fn zero(domain: usize, codomain: usize) -> Self {
        Self {
            data: vec![R::zero(); domain * codomain],
            domain,
            codomain,
        }
    }
    
    fn extend_one_row(&mut self) {
        self.codomain += 1;
        self.data.extend(vec![R::zero(); self.domain]);
    }
}

#[test]
fn test_rref_kernel() {
    use crate::linalg::field::Fp;

    type TestField = Fp<23>;

    let matrix = FlatMatrix {
        data: vec![
            TestField { 0: 1 },
            TestField { 0: 0 },
            TestField { 0: 22 },
            TestField { 0: 0 },
            TestField { 0: 1 },
            TestField { 0: 1 },
        ],
        domain: 3,
        codomain: 2,
    };
    let kernel = matrix.rref_kernel();
    let expected = FlatMatrix {
        data: vec![TestField { 0: 1 }, TestField { 0: 22 }, TestField { 0: 1 }],
        domain: 3,
        codomain: 1,
    };
    assert_eq!(kernel, expected);
}
