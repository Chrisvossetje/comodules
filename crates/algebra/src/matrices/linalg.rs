use crate::{abelian::Abelian, field::Field, matrices::flat_matrix::FlatMatrix, matrix::Matrix};


impl<F: Field> Abelian<F> for FlatMatrix<F> {
    // TODO: Should i make this usize ???
    type Module = ();

    fn kernel(&self, _domain: &Self::Module, _codomain: &Self::Module) -> (Self, Self::Module) {
        let mut clone = self.clone();
        clone.rref();
        let mut kernel = clone.rref_kernel();
        kernel.rref();
        (kernel, ())
    }
    
    fn cokernel(&self, _codomain: &Self::Module) -> (Self, Self, Self::Module) {
        // TODO : This cokernel should do something with pivots
        let (coker, module) = self.transpose().kernel(&(), &());
        (coker, FlatMatrix::zero(0, 0), module)
    }
    
    fn cohomology(_f: &Self, _g: &Self, _n: &Self::Module, _q: &Self::Module) -> (Self, Self::Module) {
        unimplemented!()
    }
    
    fn compose(f: &Self, g: &Self, _g_codomain: &Self::Module) -> Self {
        f.compose(g)
    }
    
    fn kernel_generators(&self, _domain: &Self::Module, _codomain: &Self::Module) -> Vec<usize> {
        // TODO : This could prob be smarter 
        let mut pivots = vec![];
        let mut mat = self.clone();
        while let Some(pivot) = mat.kernel_find_single_generator() {
            pivots.push(pivot);
            let codom = mat.codomain;
            mat.extend_one_row();
            mat.set(pivot, codom, F::one());
        }
        pivots
    }
}

impl<F: Field> FlatMatrix<F> {
    pub(crate) fn kernel_find_single_generator(&self) -> Option<usize> {
        let (kernel, _) = self.kernel(&(), &());
        kernel.first_non_zero_entry().map(|(_, y)| y)
    } 

    pub(crate) fn rref(&mut self) {
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

    // TODO : pub(crate)
    pub(crate) fn pivots(&self) -> Vec<(usize, usize)> {
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

    // TODO : Remove
    pub(crate) fn first_non_zero_entry(&self) -> Option<(usize, usize)> {
        for codom_id in 0..self.codomain {
            for dom_id in 0..self.domain {
                if !self.get_element(codom_id, dom_id).is_zero() {
                    return Some((codom_id, dom_id));
                }
            }
        }
        None
    }   

    pub(crate) fn rref_kernel(&self) -> Self {
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
