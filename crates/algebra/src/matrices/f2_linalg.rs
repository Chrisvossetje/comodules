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

        println!("{:?}\n{:?}\n{:?}", &self, coker, p);
        
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
        kernel.first_non_zero_entry().map(|(_, y)| y)
    } 

    /// Perform Reduced Row Echelon Form (RREF) on the matrix in-place
    /// Optimized for F2 fields where addition is XOR and multiplication is AND
    pub fn rref(&mut self) {
        let mut lead = 0;

        for r in 0..self.codomain() {
            if lead >= self.domain() {
                break;
            }

            // Find pivot row
            let mut i = r;
            while self.get_element(i, lead) == F2::zero() {
                i += 1;
                if i == self.codomain() {
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
                if i != r && self.get_element(i, lead) == F2::one() {
                    // Add row r to row i (XOR operation in F2)
                    self.add_row_to_row(r, i);
                }
            }

            lead += 1;
        }
    }

    /// Swap two rows in the matrix
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

    /// Add row `from` to row `to` (XOR operation in F2)
    fn add_row_to_row(&mut self, from: usize, to: usize) {
        let words_per_row = (self.domain() + 63) >> 6;
        let start_from = from * words_per_row;
        let start_to = to * words_per_row;
        
        for i in 0..words_per_row {
            self.data[start_to + i] ^= self.data[start_from + i];
        }
    }

    /// Get the pivot positions (column, row) after RREF
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
                if self.get_element(codom_id, dom_id) == F2::one() {
                    return Some((codom_id, dom_id));
                }
            }
        }
        None
    }   

    /// Compute the kernel (null space) of the matrix after RREF
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

            // Set the dependent variables based on the RREF form
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

    /// Solve the linear system Ax = b for x, returns None if no solution exists
    pub fn solve(&self, b: &Self) -> Option<Self> {
        if self.codomain() != b.codomain() {
            return None; // Incompatible dimensions
        }

        // Augment matrix [A|b]
        let mut augmented = Self::zero(self.domain() + b.domain(), self.codomain());
        
        // Copy A
        for i in 0..self.domain() {
            for j in 0..self.codomain() {
                augmented.set_element(i, j, self.get_element(i, j));
            }
        }
        
        // Copy b
        for i in 0..b.domain() {
            for j in 0..b.codomain() {
                augmented.set_element(self.domain() + i, j, b.get_element(i, j));
            }
        }

        // Perform RREF on augmented matrix
        augmented.rref();

        // Check for inconsistency (pivot in augmented part)
        let a_pivots = self.rank();
        let aug_pivots = augmented.rank();
        
        if aug_pivots > a_pivots {
            return None; // No solution
        }

        // Extract solution
        let mut solution = Self::zero(self.domain(), b.domain());
        let pivots = augmented.pivots();
        
        for (pivot_col, pivot_row) in pivots {
            if pivot_col < self.domain() {
                for sol_col in 0..b.domain() {
                    solution.set_element(pivot_col, sol_col, 
                        augmented.get_element(pivot_row, self.domain() + sol_col));
                }
            }
        }

        Some(solution)
    }
}