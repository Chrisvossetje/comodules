use crate::{abelian::Abelian, matrices::f2_matrix::F2Matrix, matrix::Matrix, rings::finite_fields::F2, ring::CRing};

impl Abelian<F2> for F2Matrix {
    type Generator = ();

    // TODO: WTF am i doing here ?
    fn kernel(&self, _domain: &Vec<Self::Generator>, _codomain: &Vec<Self::Generator>) -> (Self, Vec<Self::Generator>) {
        let mut clone = self.clone();
        clone.rref();
        let mut kernel = clone.rref_kernel();
        kernel.rref();
        let len = kernel.codomain();
        (kernel, vec![(); len])
    }
    
    fn cokernel(&self, _codomain: &Vec<Self::Generator>) -> (Self, Self, Vec<Self::Generator>) {
        let (coker, module) = self.transpose().kernel(&vec![], &vec![]);
        let p = coker.pivots(); // TODO: Make it clear wtf pivots returns
        
        let mut repr_vecs = F2Matrix::zero(coker.codomain(), coker.domain());
        for (domain, codomain) in p {
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
        // TODO : This could prob be smarter 
        let mut pivots = vec![];
        let mut mat = self.clone();
        while let Some(pivot) = mat.kernel_find_single_generator() {
            pivots.push(pivot);
            let codom = mat.codomain();
            mat.extend_one_row();
            mat.set(pivot, codom, F2::one());
        }
        pivots
    }
}

impl F2Matrix {
    pub(crate) fn kernel_find_single_generator(&self) -> Option<usize> {
        let (kernel, _) = self.kernel(&vec![], &vec![]);
        kernel.first_non_zero_entry().map(|(x, _)| x)
    } 

    /// Perform Row Echelon Form (rref) on the matrix in-place
    /// Optimized for F2 fields where addition is XOR and multiplication is AND
    pub fn rref(&mut self) {
        let mut lead = 0;

        for r in 0..self.codomain() {
            if lead >= self.domain() {
                break;
            }

            // Find pivot row
            let mut i = r;
            while self.get_element(lead, i) == F2::zero() {
                i += 1;
                if i == self.codomain() {
                    // panic!("WTF");
                    i = r;
                    lead += 1;
                    if lead == self.domain() {
                        return;
                    }
                }
            }

            // Swap rows if needed
            if i != r {
                self.swap_rows(r, i);
            }

            // Since we're working over F2, the pivot is always 1 (no need to normalize)
            // Eliminate other entries in this column
            for i in 0..self.codomain() {
                if i != r && self.get_element(lead, i) == F2::one() {
                    // Add row r to row i (XOR operation in F2)
                    self.add_row_to_row(r, i);
                }
            }

            lead += 1;
        }
    }

    /// Add row `from` to row `to` (XOR operation in F2)
    pub(crate) fn add_row_to_row(&mut self, from: usize, to: usize) {
        let words_per_row = (self.domain() + 63) >> 6;
        let start_from = from * words_per_row;
        let start_to = to * words_per_row;
        
        for i in 0..words_per_row {
            self.data[start_to + i] ^= self.data[start_from + i];
        }
    }

    /// Get the pivot positions (column, row) after rref
    pub fn pivots(&self) -> Vec<(usize, usize)> {
        let mut col = 0;
        let mut pivots = vec![];
        for row in 0..self.codomain() {
            while col < self.domain() {
                if self.get_element(col, row) == F2::one() {
                    pivots.push((col, row));
                    col += 1;
                    break;
                }
                col += 1;
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
        let pivot_cols: Vec<usize> = self.pivots().iter().map(|x| x.0).collect();

        // Find free variables (columns without pivots)
        for j in 0..self.domain() {
            if !pivot_cols.contains(&j) {
                free_vars.push(j);
            }
        }

        let mut kernel = Self::zero(self.domain(), free_vars.len());

        for (i, &free_var) in free_vars.iter().enumerate() {
            // Set the free variable to 1
            kernel.set_element(free_var, i, F2::one());

            // Set the dependent variables based on the rref form
            for (row, &pivot_col) in pivot_cols.iter().enumerate() {
                if row < self.codomain() && free_var < self.domain() {
                    // In F2, negation is the same as the value itself
                    kernel.set_element(pivot_col, i, self.get_element(free_var, row));
                }
            }
        }

        kernel
    }

    /// Get the rank of the matrix (number of pivots)
    pub fn rank(&self) -> usize {
        let mut clone = self.clone();
        clone.rref();
        clone.pivots().len()
    }

    /// Get the nullity (dimension of kernel) of the matrix
    pub fn nullity(&self) -> usize {
        self.domain() - self.rank()
    }

    /// Check if the matrix is in reduced row echelon form
    pub fn is_rref(&self) -> bool {
        let pivots = self.pivots();
        
        // Check that pivot columns are increasing
        for i in 1..pivots.len() {
            if pivots[i].0 <= pivots[i-1].0 {
                return false;
            }
        }

        // Check that each pivot is the only non-zero entry in its column
        for (pivot_col, pivot_row) in pivots {
            for row in 0..self.codomain() {
                if row != pivot_row && self.get_element(row, pivot_col) != F2::zero() {
                    return false;
                }
            }
        }

        true
    }
}