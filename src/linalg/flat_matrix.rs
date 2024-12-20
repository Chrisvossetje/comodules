use serde::{Deserialize, Serialize};

use super::{field::Field, matrix::Matrix};

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct FlatMatrix<F: Field> {
    pub data: Vec<F>,
    pub domain: usize,
    pub codomain: usize,
}

impl<F: Field> FlatMatrix<F> {
    fn get_element(&self, row: usize, col: usize) -> F {
        self.data[row * self.domain + col].clone()
    }

    fn set_element(&mut self, row: usize, col: usize, value: F) {
        self.data[row * self.domain + col] = value;
    }

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

    fn transpose(&self) -> Self {
        let mut new_matrix = vec![F::zero(); self.data.len()];
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

    fn vstack(&self, other: &Self) -> Self {
        assert_eq!(
            self.domain(),
            other.domain(),
            "Domains of the two matrices do not have the same dimension"
        );
        let mut new = self.clone();

        new.data.extend_from_slice(&other.data);
        new.codomain += other.codomain;
        new
    }

    fn block_sum(&self, other: &Self) -> Self {
        let new_domain = self.domain + other.domain;
        let new_codomain = self.codomain + other.codomain;
        let mut new = Self::zero(new_domain, new_codomain);

        for i in 0..self.codomain {
            for j in 0..self.domain {
                new.data[i * (new_domain) + j] = self.get_element(i, j);
            }
        }

        for i in 0..other.codomain {
            for j in 0..other.domain {
                new.data[(self.codomain + i) * (new_domain) + (self.domain + j)] =
                    other.get_element(i, j);
            }
        }

        new
    }

    fn identity(d: usize) -> Self {
        let mut data = vec![F::zero(); d * d];
        for i in 0..d {
            data[i * d + i] = F::one();
        }
        Self {
            data,
            domain: d,
            codomain: d,
        }
    }

    fn get(&self, domain: usize, codomain: usize) -> F {
        self.get_element(codomain, domain)
    }

    fn set(&mut self, domain: usize, codomain: usize, f: F) {
        self.set_element(codomain, domain, f);
    }

    fn add_at(&mut self, domain: usize, codomain: usize, f: F) {
        let idx = codomain * self.domain + domain;
        self.data[idx] += f;
    }

    fn get_row(&self, codomain: usize) -> &[F] {
        let start = codomain * self.domain;
        let end = start + self.domain;
        &self.data[start..end]
    }

    fn set_row(&mut self, codomain: usize, row: &[F]) {
        let start = codomain * self.domain;
        let end = start + self.domain;
        self.data[start..end].copy_from_slice(row);
    }

    fn domain(&self) -> usize {
        self.domain
    }

    fn codomain(&self) -> usize {
        self.codomain
    }

    fn compose(&self, rhs: &Self) -> Self {
        assert_eq!(
            self.domain, rhs.codomain,
            "Matrix domain not equal to codomain"
        );

        let mut compose = vec![F::zero(); self.codomain * rhs.domain];

        for x in 0..self.codomain {
            for y in 0..rhs.domain {
                let mut sum = F::zero();
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
            data: vec![F::zero(); domain * codomain],
            domain,
            codomain,
        }
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
