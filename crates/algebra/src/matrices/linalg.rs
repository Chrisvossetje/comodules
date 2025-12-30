use crate::{abelian::Abelian, field::Field, matrices::{flat_matrix::FlatMatrix}, matrix::Matrix, ring::CRing};

impl<F: Field> Abelian<F> for FlatMatrix<F> {
    type Generator = ();

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
        let p = coker.pivots();

        let mut repr_vecs = FlatMatrix::zero(coker.codomain, coker.domain);
        for (domain, codomain) in p {
            repr_vecs.set(codomain, domain, F::one());
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
            mat.set(pivot, codom, F::one());
        }
        pivots
    }
}



impl<F: Field> FlatMatrix<F> {
    pub(crate) fn kernel_find_single_generator(&self) -> Option<usize> {
        let (kernel, _) = self.kernel(&vec![], &vec![]);
        kernel.first_non_zero_entry().map(|(x, _)| x)
    } 

    pub(crate) fn rref(&mut self) {
        let mut lead = 0;

        for r in 0..self.codomain {
            if lead >= self.domain {
                break;
            }

            let mut i = r;
            while self.get_element(lead, i).is_zero() {
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

            let pivot = self.get_element(lead, r);
            if !pivot.is_zero() {
                let pivot_inv = pivot.inv().expect("Pivot should be invertible");
                for j in 0..self.domain {
                    let idx = r * self.domain + j;
                    self.data[idx] *= pivot_inv;
                }
            }

            for i in 0..self.codomain {
                if i != r {
                    let factor = self.get_element(lead, i);
                    for j in 0..self.domain {
                        let idx = i * self.domain + j;
                        let el = self.get_element(j, r);
                        self.data[idx] -= factor * el;
                    }
                }
            }

            lead += 1;
        }
    }

    pub(crate) fn pivots(&self) -> Vec<(usize, usize)> {
        let mut domain = 0;
        let mut pivots = vec![];
        for codomain in 0..self.codomain {
            while domain < self.domain {
                if !self.get_element(domain, codomain).is_zero() {
                    pivots.push((domain, codomain));
                    domain += 1;
                    break;
                }
                domain += 1;
            }
        }
        pivots
    }

    pub(crate) fn first_non_zero_entry(&self) -> Option<(usize, usize)> {
        for codom_id in 0..self.codomain {
            for dom_id in 0..self.domain {
                if !self.get_element(dom_id, codom_id).is_zero() {
                    return Some((dom_id, codom_id));
                }
            }
        }
        None
    }   

    pub(crate) fn rref_kernel(&self) -> Self {
        let mut free_vars = Vec::new();
        let pivot_doms: Vec<usize> = self.pivots().iter().map(|x| x.0).collect();

        for j in 0..self.domain {
            if !pivot_doms.contains(&j) {
                free_vars.push(j);
            }
        }

        let mut kernel = vec![F::zero(); free_vars.len() * self.domain];

        for (i, &free_var) in free_vars.iter().enumerate() {
            let codom_idx = i * self.domain;
            kernel[codom_idx + free_var] = F::one();

            for (codom, &pivot_dom) in pivot_doms.iter().enumerate() {
                kernel[codom_idx + pivot_dom] = -self.get_element(free_var, codom);
            }
        }

        Self {
            data: kernel,
            domain: self.domain,
            codomain: free_vars.len(),
        }
    }
}
