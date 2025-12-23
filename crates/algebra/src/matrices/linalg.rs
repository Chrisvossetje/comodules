use crate::{abelian::Abelian, field::Field, matrices::{f2_matrix::F2Matrix, flat_matrix::FlatMatrix}, matrix::Matrix, ring::CRing, rings::finite_fields::F2};

// TODO:
// impl<F: Field> Abelian<F> for FlatMatrix<F> {
impl Abelian<F2> for FlatMatrix<F2> {
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

        let mut repr_vecs = FlatMatrix::zero(coker.codomain, coker.domain);
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
            let codom = mat.codomain;
            mat.extend_one_row();
            mat.set(pivot, codom, F2::one());
        }
        pivots
    }
}



// impl<F: Field> FlatMatrix<F> {
impl FlatMatrix<F2> {
    pub(crate) fn kernel_find_single_generator(&self) -> Option<usize> {
        let (kernel, _) = self.kernel(&vec![], &vec![]);
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

        let mut kernel = vec![F2::zero(); free_vars.len() * self.domain];

        for (i, &free_var) in free_vars.iter().enumerate() {
            let row_idx = i * self.domain;
            kernel[row_idx + free_var] = F2::one();

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
