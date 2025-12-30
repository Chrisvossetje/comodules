use serde::{Deserialize, Serialize};
use std::fmt::Debug;

use crate::{matrix::Matrix, ring::CRing};
use deepsize::DeepSizeOf;



#[derive(Clone, PartialEq, Deserialize, Serialize, DeepSizeOf)]
pub struct FlatMatrix<R: CRing> {
    pub(crate) data: Vec<R>,
    pub(crate) domain: usize,
    pub(crate) codomain: usize,
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
    fn get_index(&self, domain: usize, codomain: usize) -> usize {
        codomain * self.domain + domain
    }

    pub(crate) fn get_element(&self, domain: usize, codomain: usize) -> R {
        self.data[codomain * self.domain + domain].clone()
    }

    pub(crate) fn set_element(&mut self, domain: usize, codomain: usize, value: R) {
        self.data[codomain * self.domain + domain] = value;
    }
}

impl<R: CRing> Matrix<R> for FlatMatrix<R> {
    type UnderlyingRowType = R;
    
    fn transpose(&self) -> Self {
        let mut new_matrix = vec![R::zero(); self.data.len()];
        for i in 0..self.codomain {
            for j in 0..self.domain {
                new_matrix[j * self.codomain + i] = self.get_element(j, i);
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
        self.get_element(domain, codomain)
    }

    fn set(&mut self, domain: usize, codomain: usize, r: R) {
        self.set_element(domain, codomain, r);
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
                    sum += self.get_element(k, x) * rhs.get_element(y, k);
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

    fn swap_rows(&mut self, codom1: usize, codom2: usize) {
        if codom1 == codom2 {
            return;
        }
        for j in 0..self.domain() {
            let a = self.get_index(j, codom1);
            let b = self.get_index(j, codom2);
            self.data.swap(a, b);
        }
    }

    fn swap_cols(&mut self, domain1: usize, domain2: usize) {
        if domain1 == domain2 {
            return;
        }
        for i in 0..self.codomain() {
            let a = self.get_index(domain1, i);
            let b = self.get_index(domain2, i);
            self.data.swap(a, b);
        }
    }
    
    fn eval_vector(&self, vector: &[R]) -> Vec<R> {
        let mut m = FlatMatrix::zero(1, vector.len());
        for (id,v) in vector.iter().enumerate() {
            m.set(0, id, *v);
        }
        self.compose(&m).data
    }
}
