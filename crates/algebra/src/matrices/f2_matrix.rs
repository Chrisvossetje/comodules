use crate::{matrices::flat_matrix::FlatMatrix, matrix::Matrix, ring::CRing, rings::finite_fields::F2};
use serde::{Deserialize, Serialize};
use deepsize::DeepSizeOf;


#[derive(Clone, PartialEq, Deserialize, Serialize, DeepSizeOf)]
pub struct F2Matrix {
    pub data: Vec<u64>,
    pub(crate) domain: usize,
    pub(crate) codomain: usize,
    pub(crate) words_per_row: usize,

    // For each domain, find if it contains a pivot at column i
    pub(crate) pivots: Vec<Option<u32>>,
}

impl std::fmt::Debug for F2Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for m in 0..self.codomain {
            for n in 0..self.domain {
                let el = self.get(n, m);
                let s = format!("{:?}", el);
                write!(f, "{}", s)?;
            }
            writeln!(f)?;
        }
        writeln!(f)
    }
}

impl F2Matrix {
    pub(crate) fn index(&self, domain: usize, codomain: usize) -> (usize, usize) { // index + shift
        let word_idx = domain >> 6; // Which 64-bit word
        let bit_idx = domain & 0b111111; // Which bit in that word (0-63)
        
        let base_word = codomain * self.words_per_row;
         
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

    pub fn from_flat(flat: FlatMatrix<F2>) -> Self {
        let mut z = Self::zero(flat.domain, flat.codomain);
        for i in 0..flat.domain() {
            for j in 0..flat.codomain {
                z.set(i,j,flat.get(i, j));
            } 
        }
        z
    }

    pub fn to_flat(self) -> FlatMatrix<F2> {
        let mut z = FlatMatrix::zero(self.domain, self.codomain);
        for i in 0..self.domain() {
            for j in 0..self.codomain {
                z.set(i,j,self.get(i, j));
            } 
        }
        z
    }


    pub fn get_row_mut(&mut self, codomain: usize) -> &mut [u64] {
        let start = self.words_per_row * codomain;
        let end = start + self.words_per_row;
        &mut self.data[start..end]
    }
}

impl Matrix<F2> for F2Matrix {
    type UnderlyingRowType = u64;

    fn zero(domain: usize, codomain: usize) -> Self {
        let words_per_row = (domain + 63) >> 6; // Ceiling division by 64
        Self { data: vec![0; words_per_row * codomain], domain, codomain, words_per_row, pivots: vec![] }
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
        let start = self.words_per_row * codomain;
        let end = start + self.words_per_row;

        &self.data[start..end]
    }
    
    fn set_row(&mut self, codomain: usize, row: &[u64]) {
        let start = self.words_per_row * codomain;
        let end = start + self.words_per_row;

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

        /// Evaluate `self * vector` (vector length = domain), returning length = codomain.
    /// Uses packed parity computation.
    fn eval_vector(&self, vector: &[F2]) -> Vec<F2> {
        assert_eq!(vector.len(), self.domain);

        // Pack vector bits into u64 words matching matrix layout.
        let v_words = (self.domain + 63) >> 6;
        let mut v = vec![0u64; v_words];
        for (i, bit) in vector.iter().enumerate() {
            if bit.0 & 1 == 1 {
                v[i >> 6] |= 1u64 << (i & 63);
            }
        }

        let mut out = vec![F2(0); self.codomain];
        for row in 0..self.codomain {
            let r = self.get_row(row);
            let mut acc = 0u64;
            for w in 0..self.words_per_row {
                acc ^= r[w] & v[w];
            }
            // parity of acc
            let parity = (acc.count_ones() & 1) as u8;
            out[row] = F2(parity);
        }
        out
    }


    /// [ self  0 ]
    /// [  0  other ]
    fn block_sum(&mut self, other: &Self) {
        let old_domain = self.domain;
        let old_codomain = self.codomain;

        let new_domain = self.domain + other.domain;
        let new_codomain = self.codomain + other.codomain;

        let new_wpr = (new_domain + 63) >> 6;
        let mut new_data = vec![0u64; new_wpr * new_codomain];

        // Copy self rows into top-left.
        for row in 0..old_codomain {
            let src = self.get_row(row);
            let dst_start = row * new_wpr;
            // domain offset = 0, so aligned copy with possible truncation
            new_data[dst_start..dst_start + self.words_per_row].copy_from_slice(src);
        }

        // Copy other rows into bottom-right with bit offset.
        let bit_off = old_domain & 63;
        let word_off = old_domain >> 6;

        for row in 0..other.codomain {
            let src = other.get_row(row);
            let dst_row = old_codomain + row;
            let dst_start = dst_row * new_wpr;

            if bit_off == 0 {
                // aligned at word boundary
                for w in 0..other.words_per_row {
                    new_data[dst_start + word_off + w] |= src[w];
                }
            } else {
                // needs shifting across word boundary
                let sh = bit_off as u32;
                let inv_sh = 64 - sh;

                for w in 0..other.words_per_row {
                    let val = src[w];
                    let low = val << sh;
                    new_data[dst_start + word_off + w] |= low;

                    // spill into next word if needed
                    if dst_start + word_off + w + 1 < dst_start + new_wpr {
                        let high = val >> inv_sh;
                        new_data[dst_start + word_off + w + 1] |= high;
                    }
                }
            }
        }

        self.data = new_data;
        self.domain = new_domain;
        self.codomain = new_codomain;
        self.words_per_row = new_wpr;
    }

    fn extend_one_row(&mut self) {
        for _ in 0..self.words_per_row {
            self.data.push(0);
        }
        
        self.codomain += 1;
    }

    fn swap_rows(&mut self, row1: usize, row2: usize) {
        if row1 == row2 {
            return;
        }
        
        let words_per_row = (self.domain() + 63) >> 6;
        let start1 = row1 * words_per_row;
        let start2 = row2 * words_per_row;
        
        for i in 0..words_per_row {
            self.data.swap(start1 + i, start2 + i);
        }
    }
}