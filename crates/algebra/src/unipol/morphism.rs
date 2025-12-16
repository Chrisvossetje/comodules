
use std::cmp::Ordering;

use itertools::Itertools;

use crate::{abelian::Abelian, field::Field, matrices::flat_matrix::FlatMatrix, matrix::Matrix, ring::CRing, rings::univariate_polynomial_ring::UniPolRing, snf::SmithNormalForm, unipol::{UniPolMap, UniPolModule, cohomology::internal_cohomology}};

impl<F: Field> Abelian<UniPolRing<F>> for UniPolMap<F> {
    type Generator = Option<u16>;
    
    fn kernel(&self, _domain: &Vec<Self::Generator>, _codomain: &Vec<Self::Generator>) -> (Self, Vec<Self::Generator>) {
        todo!()
    }
    
    fn cokernel(&self, codomain: &Vec<Self::Generator>) -> (Self, Self, Vec<Self::Generator>) {
        self.internal_cokernel(codomain)
    }
    
    fn cohomology(f: &Self, g: &Self, n: &Vec<Self::Generator>, q: &Vec<Self::Generator>) -> (Self, Vec<Self::Generator>) {
        internal_cohomology(f,g,n,q)
    }


    fn compose(f: &Self, g: &Self, g_codomain: &Vec<Self::Generator>) -> Self {
        let mut comp = f.compose(g);
        comp.reduce(g_codomain);
        comp
    }
    
    fn kernel_destroyers(&self, domain: &Vec<Self::Generator>, codomain: &Vec<Self::Generator>) -> Vec<usize> {
        let mut pivots = vec![];
        let mut mat = self.clone();
        let mut codomain = codomain.clone();
        while let Some(pivot) = mat.kernel_find_single_generator(domain, &codomain) {
            pivots.push(pivot);
            let codom = mat.codomain;
            mat.extend_one_row();
            mat.set(pivot, codom, UniPolRing::one());

            let pivot_structure = domain[pivot];
            codomain.push(pivot_structure);
        }
        pivots
    }
}

impl<F: Field> UniPolMap<F> {
    pub fn verify(&self, domain: &UniPolModule, codomain: &UniPolModule) -> Result<(), String> {
        if self.domain != domain.len() {
            return Err(format!("The domain size of the map does not equal that of the module"));
        }
        if self.codomain != codomain.len() {
            return Err(format!("The codomain size of the map does not equal that of the module"));
        }

        for dom in 0..self.domain {
            for codom in 0..self.codomain {
                let el = self.get(dom, codom);
                if !el.is_zero() {
                    let power = el.1;
                    let dom_module = domain[dom];
                    let codom_module = codomain[codom];

                    if let Some(codom_module_power) = codom_module {
                        if power >= codom_module_power {
                            return Err(format!("There is an element mapping to some power of a codomain, but the codomain has quotient which divides out this power. Meaning the matrix should have been reduced and this should have shown zero. Domain {dom} id, codomain {codom} id, power {power}."));
                        }
                    }

                    match dom_module {
                        Some(dom_cycle) => {
                            match codom_module {
                                Some(codom_cycle) => {
                                    if dom_cycle + power < codom_cycle {
                                        return Err(format!("Maps a cyclic t^{dom_cycle} module to a cyclic t^{codom_cycle} module, which is illegal with power between those being {power}"));   
                                    }
                                },
                                None => {
                                    return Err(format!("Cannot have a (non-zero) map from a cyclic module to a free module."));
                                },
                            }
                        },
                        None => {},
                    }
                }
            }
        }
        Ok(())
    }
    
    
    pub(super) fn reduce(&mut self, codomain: &UniPolModule) {
        for (y, pow) in codomain.iter().enumerate() {
            match pow {
                Some(pow) => {
                    for x in 0..self.domain {
                        if self.get(x, y).1 >= *pow {
                            self.set(x, y, UniPolRing::zero());
                        }
                    }
                },
                None => {},
            }
        }
    }


    pub(crate) fn kernel_find_single_generator(&self, domain: &UniPolModule, codomain: &UniPolModule) -> Option<usize> {
        debug_assert_eq!(self.codomain, codomain.len());
        
        if self.codomain == 0 {
            if self.domain == 0 {
                return None;
            }
            else {
                return Some(0);
            }
        }    
    
        let new_domain = self.domain + self.codomain;
        let mut new_map = FlatMatrix::zero(new_domain, self.codomain);
        
        for n in  0..new_map.codomain {
            new_map.set_row(n, self.get_row(n));
            
            
            if let Some(power) = codomain[n] {
                new_map.set(self.domain + n, n, UniPolRing(F::one().neg(),power));
            }
        }
        
        let (_, s, v, _, _) = new_map.full_snf();

        let mut possible_id = None;
        let mut possible_pow = u16::MAX;

        // first verify if the i'th vector in s is in the kernel.
        // And if it is check that vector is non zero
        for r in 0..s.codomain {
            let el = s.get(r,r);
            
            if el.is_zero() {
                // consider this vector
                for y in 0..self.domain {
                    let domain_pow = domain[y];
                    let super_el = v.get(r, y);
                    if !super_el.is_zero() {
                        if super_el.1 < possible_pow {
                            if let Some(domain_pow) = domain_pow {
                                if domain_pow <= super_el.1 {
                                    continue;
                                } 
                            }
                            possible_id = Some(y);
                            possible_pow = super_el.1;
                        } 
                    }
                }
            }            
        }
        // The vectors outside of s.codomain will always map to zero
        // verify if there truncation is still a non-zero vector
        // If it is choose one original basis element as a pivot 
        for r in s.codomain..s.domain {
            for y in 0..self.domain {
                let domain_pow = domain[y];
                let super_el = v.get(r, y);
                if !super_el.is_zero() {

                    if let Some(domain_pow) = domain_pow {
                        if domain_pow <= super_el.1 {
                            continue;
                        } 
                    }

                    if super_el.1 < possible_pow {
                        possible_id = Some(y);
                        possible_pow = super_el.1;
                    } 
                }
            } 
        }

        possible_id
    }

} 


pub fn order_maps<F: Field>(f: &UniPolMap<F>, g: &UniPolMap<F>, n: &UniPolModule, q: &UniPolModule) -> (UniPolMap<F>, UniPolMap<F>, UniPolModule, UniPolModule, FlatMatrix<UniPolRing<F>>) {
    let n_sorted: Vec<_> = n.iter().enumerate().sorted_by(|a, b|{
        match a.1 {
            Some(pow) => {
                match b.1 {
                    Some(b_pow) => {
                        if pow < b_pow {
                            Ordering::Greater
                        } else if pow == b_pow {
                            Ordering::Equal
                        } else {
                            Ordering::Less
                        }
                    },
                    None => {
                        Ordering::Greater
                    },
                } 
            },
            None => {
                match b.1 {
                    Some(_) => {
                        Ordering::Less
                    },
                    None => {
                        Ordering::Equal
                    },
                }
            },
        }
    }).collect();

    let q_sorted: Vec<_> = q.iter().enumerate().sorted_by(|a, b|{
        match a.1 {
            Some(pow) => {
                match b.1 {
                    Some(b_pow) => {
                        if pow < b_pow {
                            Ordering::Greater
                        } else if pow == b_pow {
                            Ordering::Equal
                        } else {
                            Ordering::Less
                        }
                    },
                    None => {
                        Ordering::Greater
                    },
                } 
            },
            None => {
                match b.1 {
                    Some(_) => {
                        Ordering::Less
                    },
                    None => {
                        Ordering::Equal
                    },
                }
            },
        }
    }).collect();


    let mut n_trans_map = FlatMatrix::zero(f.codomain, f.codomain);    
    let mut n_trans_map_inv = FlatMatrix::zero(f.codomain, f.codomain);    
    
    for (target, (origin, _)) in n_sorted.iter().enumerate() {
        n_trans_map.set(*origin, target, UniPolRing::one());
        n_trans_map_inv.set(target, *origin, UniPolRing::one());
    }
    
    let mut q_trans_map = FlatMatrix::zero(g.codomain, g.codomain);   

    for (target, (origin, _)) in q_sorted.iter().enumerate() {
        q_trans_map.set(*origin, target, UniPolRing::one());
    }
    
    let new_f = n_trans_map.compose(&f);
    let new_g = q_trans_map.compose(&g.compose(&n_trans_map_inv));

    let new_n = n_sorted.iter().map(|x| x.1.clone()).collect();
    let new_q = q_sorted.iter().map(|x| x.1.clone()).collect();
    (new_f, new_g, new_n, new_q, n_trans_map_inv)
}

