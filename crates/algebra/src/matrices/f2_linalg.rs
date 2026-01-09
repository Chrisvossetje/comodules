use crate::{abelian::Abelian, matrices::{f2_matrix::F2Matrix}, matrix::Matrix, ring::CRing, rings::finite_fields::F2};

impl Abelian<F2> for F2Matrix {
    type Generator = ();

    fn kernel(&self, _domain: &Vec<Self::Generator>, _codomain: &Vec<Self::Generator>) -> (Self, Vec<Self::Generator>) {
        
        let mut clone = self.clone();
        clone.echelonize();

        let mut kernel = clone.rref_kernel();
        kernel.echelonize();
        let len = kernel.codomain();
        (kernel, vec![(); len])
    }

    fn cokernel(&self, codomain: &Vec<Self::Generator>) -> (Self, Self, Vec<Self::Generator>) {                
        return self.transpose().transposed_cokernel(&codomain);
    }
    
    fn transposed_cokernel(&self, _codomain: &Vec<Self::Generator>) -> (Self, Self, Vec<Self::Generator>) {
        let (coker, module) = self.kernel(&vec![], &vec![]);
        
        let mut repr_vecs = F2Matrix::zero(coker.codomain(), coker.domain());
        for (domain, codomain) in coker.pivots.iter().enumerate().filter_map(|x| {if x.1.is_some() {Some((x.0, x.1.unwrap() as usize))} else {None}}) {
            repr_vecs.set(codomain, domain, F2::one());
        }
        
        debug_assert!(coker.compose(&repr_vecs).is_unit().is_ok());

        (coker, repr_vecs, module)
    }
    
    fn cohomology(_f: &Self, _g: &Self, _n: &Vec<Self::Generator>, _q: &Vec<Self::Generator>) -> (Self, Vec<Self::Generator>) {
        unimplemented!()
    }
    
    fn compose(f: &Self, g: &Self, _g_codomain: &Vec<Self::Generator>) -> Self {
        f.compose(g)
    }
    
    fn kernel_destroyers(&self, _domain: &Vec<Self::Generator>, _codomain: &Vec<Self::Generator>) -> Vec<usize> {
        let mut pivots = vec![];
        let (mut kernel, _) = self.kernel(&vec![], &vec![]);

        while let Some((pivot_domain, pivot_codomain)) = kernel.first_non_zero_entry() {
            pivots.push(pivot_domain);
            
            for domain in 0..kernel.domain() {
                kernel.set(domain, pivot_codomain, F2::zero());
            }
            for codom in 0..kernel.codomain() {
                kernel.set(pivot_domain, codom, F2::zero());
            }
        }
        pivots
    }
}





impl F2Matrix {
    pub(crate) fn echelonize(&mut self) {
        if self.codomain() < 64 {
            return self.echelonize_naive(true);
        } else {
            return self.echelonize_m4ri();
        }
    }

    pub(crate) fn xor_row_from_word(&mut self, dst_row: usize, src_row: usize, start_word: usize) {
        if dst_row == src_row {
            return;
        }
        let wpr = self.words_per_row;
        debug_assert!(start_word <= wpr);

        let dst_start = dst_row * wpr;
        let src_start = src_row * wpr;

        // Borrow two disjoint slices from self.data using split_at_mut
        if dst_start < src_start {
            let (left, right) = self.data.split_at_mut(src_start);
            let dst = &mut left[dst_start..dst_start + wpr];
            let src = &right[0..wpr];
            for w in start_word..wpr {
                dst[w] ^= src[w];
            }
        } else {
            let (left, right) = self.data.split_at_mut(dst_start);
            let src = &left[src_start..src_start + wpr];
            let dst = &mut right[0..wpr];
            for w in start_word..wpr {
                dst[w] ^= src[w];
            }
        }
    }

    /// In-place row echelonization (GF(2)).
    /// `reduced = false` => REF
    /// `reduced = true`  => RREF (Gauss-Jordan)
    pub(crate) fn echelonize_naive(&mut self, reduced: bool) {
        let rows = self.codomain;
        let cols = self.domain;

        self.pivots = vec![None; cols];

        let mut rank = 0usize;

        // Forward elimination
        for col in 0..cols {
            if rank >= rows {
                break;
            }
            let word = col >> 6;
            let bit = col & 63;
            let mask = 1u64 << bit;

            // Find pivot row
            let mut pivot_row = None;
            for r in rank..rows {
                if (self.get_row(r)[word] & mask) != 0 {
                    pivot_row = Some(r);
                    break;
                }
            }
            let Some(p) = pivot_row else { continue };

            // Swap pivot into position `rank`
            if p != rank {
                self.swap_rows(p, rank);
            }

            // Eliminate above ?
            if reduced {
                for r in 0..rank {
                    if (self.get_row(r)[word] & mask) != 0 {
                        self.xor_row_from_word(r, rank, word);
                    }
                }
            }

            // Eliminate below
            for r in (rank + 1)..rows {
                if (self.get_row(r)[word] & mask) != 0 {
                    self.xor_row_from_word(r, rank, word);
                }
            }

            self.pivots[col] = Some(rank as u32);
            rank += 1;
        }
    }

    /// Get the pivot positions (column, row) after rref
    pub fn pivots(&self) -> Vec<(usize, usize)> {
        let mut domain = 0;
        let mut pivots = vec![];
        for codomain in 0..self.codomain() {
            while domain < self.domain() {
                if self.get_element(domain, codomain) == F2::one() {
                    pivots.push((domain, codomain));
                    domain += 1;
                    break;
                }
                domain += 1;
            }
        }
        pivots
    }

    /// Find the first non-zero entry in the matrix
    pub fn first_non_zero_entry(&self) -> Option<(usize, usize)> {
        for codom_id in 0..self.codomain() {
            for dom_id in 0..self.domain() {
                if self.get_element(dom_id, codom_id) == F2::one() {
                    return Some((dom_id, codom_id));
                }
            }
        }
        None
    }   

    /// Compute the kernel (null space) of the matrix after rref
    pub fn rref_kernel(&self) -> Self {
        let mut free_vars = Vec::new();
        
        // Find free variables (columns without pivots)
        for j in 0..self.domain() {
            if !self.pivots[j].is_some() {
                free_vars.push(j);
            }
        }

        let mut kernel = Self::zero(self.domain(), free_vars.len());

        for (i, &free_var_domain) in free_vars.iter().enumerate() {
            // Set the free variable to 1
            kernel.set_element(free_var_domain, i, F2::one());

            // Set the dependent variables based on the rref form
            for (pivot_domain, &pivot_codomain) in self.pivots.iter().enumerate() {
                if let Some(pivot_codomain) = pivot_codomain {
                    // In F2, negation is the same as the value itself
                    kernel.set_element(pivot_domain, i, self.get_element(free_var_domain, pivot_codomain as usize));
                }
            }
        }

        kernel
    }

    /// Get the rank of the matrix (number of pivots)
    pub fn rank(&self) -> usize {
        let mut clone = self.clone();
        clone.echelonize();
        clone.pivots().len()
    }

    /// Get the nullity (dimension of kernel) of the matrix
    pub fn nullity(&self) -> usize {
        self.domain() - self.rank()
    }

    /// Check if the matrix is in reduced row echelon form
    pub fn is_rref(&self) -> bool {
        // Check that pivot columns are increasing
        for i in 1..self.pivots.len() { 
            if self.pivots[i].is_some() && self.pivots[i] <= self.pivots[i-1] {
                return false;
            }
        }

        // Check that each pivot is the only non-zero entry in its column
        for (pivot_row, pivot_col) in self.pivots.iter().enumerate() {
            if let Some(pivot_col) = pivot_col {
                for row in 0..self.codomain() {
                    if row == pivot_row && self.get_element(row, *pivot_col as usize) != F2::one() {
                        return false;
                    }
                    if row != pivot_row && self.get_element(row, *pivot_col as usize) != F2::zero() {
                        return false;
                    }
                }
            }
        }

        true
    }
}