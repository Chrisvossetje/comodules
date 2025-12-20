// TODO : implement an F2_matrix thing, where we compress the data of our matrices

use crate::{matrix::Matrix, rings::finite_fields::F2, ring::CRing};
use serde::{Deserialize, Serialize};
use deepsize::DeepSizeOf;


#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, DeepSizeOf)]
pub struct F2Matrix {
    pub data: Vec<u64>,
    pub(crate) domain: usize,
    pub(crate) codomain: usize,
}

impl F2Matrix {
    pub(crate) fn index(&self, domain: usize, codomain: usize) -> (usize, usize) { // index + shift
        let word_idx = domain >> 6; // Which 64-bit word
        let bit_idx = domain & 0b111111; // Which bit in that word (0-63)
        
        let words_per_row = (self.domain + 63) >> 6; // Ceiling division by 64
        let base_word = codomain * words_per_row;
         
        (base_word + word_idx, bit_idx)
    }

    pub(crate) fn get_element(&self, domain: usize, codomain: usize) -> F2 {
        let (idx, shift) = self.index(domain, codomain);        
        F2(((self.data[idx] >> shift) & 0b1) as u8)
    }
    
    pub(crate) fn set_element(&mut self, domain: usize, codomain: usize, el: F2) {
        let (idx, shift) = self.index(domain, codomain);
        self.data[idx] = (self.data[idx] & (!(1u64 << shift))) | ((el.0 as u64) << shift);
    }
}

impl Matrix<F2> for F2Matrix {
    type UnderlyingRowType = u64;

    fn zero(domain: usize, codomain: usize) -> Self {
        let words_per_row = (domain + 63) >> 6; // Ceiling division by 64
        Self { data: vec![0; words_per_row * codomain], domain, codomain }
    }

    fn get(&self, domain: usize, codomain: usize) -> F2 {
        self.get_element(domain, codomain)
    }

    fn set(&mut self, domain: usize, codomain: usize, r: F2) {
        self.set_element(domain, codomain, r);
    }

    fn add_at(&mut self, domain: usize, codomain: usize, r: F2) {
        let (idx, shift) = self.index(domain, codomain);
        self.data[idx] ^= (r.0 as u64) << shift;
    }

    fn get_row(&self, codomain: usize) -> &[u64] {
        let words_per_row = (self.domain + 63) >> 6; // Ceiling division by 64
        let start = words_per_row * codomain;
        let end = start + words_per_row;

        &self.data[start..end]
    }
    
    fn set_row(&mut self, codomain: usize, row: &[u64]) {
        let words_per_row = (self.domain + 63) >> 6; // Ceiling division by 64
        let start = words_per_row * codomain;
        let end = start + words_per_row;

        self.data[start..end].copy_from_slice(row);
    }

    fn compose(&self, rhs: &Self) -> Self {
        assert_eq!(self.domain, rhs.codomain);
        let mut result = Self::zero(rhs.domain, self.codomain);
        
        for i in 0..rhs.domain {
            for j in 0..self.codomain {
                let mut sum = F2::zero();
                for k in 0..self.domain {
                    sum = sum + rhs.get_element(i, k) * self.get_element(k, j);
                }
                result.set_element(i, j, sum);
            }
        }
        
        result
    }

    fn transpose(&self) -> Self {
        let mut result = Self::zero(self.codomain, self.domain);
        
        for i in 0..self.domain {
            for j in 0..self.codomain {
                result.set_element(j, i, self.get_element(i, j));
            }
        }
        
        result
    }

    fn domain(&self) -> usize {
        self.domain
    }

    fn codomain(&self) -> usize {
        self.codomain
    }

    fn vstack(&mut self, other: &mut Self) {
        assert_eq!(self.domain, other.domain);
        
        let new_codomain = self.codomain + other.codomain;
        
        // Extend our data with other's data
        self.data.extend_from_slice(&other.data);
        self.codomain = new_codomain;
    }

    fn block_sum(&mut self, other: &Self) {
        let old_domain = self.domain;
        let old_codomain = self.codomain;
        let new_domain = self.domain + other.domain;
        let new_codomain = self.codomain + other.codomain;
        
        let words_per_old_row = (old_domain + 63) >> 6;
        let words_per_new_row = (new_domain + 63) >> 6;
        
        let mut new_data = vec![0u64; words_per_new_row * new_codomain];
        
        // Copy self to top-left block
        for row in 0..old_codomain {
            for word in 0..words_per_old_row {
                let old_idx = row * words_per_old_row + word;
                let new_idx = row * words_per_new_row + word;
                new_data[new_idx] = self.data[old_idx];
            }
        }
        
        // Copy other to bottom-right block
        for row in 0..other.codomain {
            for col in 0..other.domain {
                let other_value = other.get_element(col, row);
                if other_value != F2::zero() {
                    let new_row = old_codomain + row;
                    let new_col = old_domain + col;
                    let word_idx = new_col >> 6;
                    let bit_idx = new_col & 0b111111;
                    let data_idx = new_row * words_per_new_row + word_idx;
                    new_data[data_idx] |= 1u64 << bit_idx;
                }
            }
        }
        
        self.data = new_data;
        self.domain = new_domain;
        self.codomain = new_codomain;
    }

    fn extend_one_row(&mut self) {
        let words_per_row = (self.domain + 63) >> 6;

        for _ in 0..words_per_row {
            self.data.push(0);
        }
        
        self.codomain += 1;
    }
    
    fn eval_vector(&self, _vector: &[F2]) -> Vec<F2> {
        todo!()
    }
}