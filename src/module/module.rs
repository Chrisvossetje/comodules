use std::cmp::Ordering;

use ahash::HashMap;
use itertools::Itertools;
use rayon::Yield;
use serde::{Deserialize, Serialize};

use crate::{
    basiselement::{BasisElement, kBasisElement}, grading::{Grading, UniGrading}, linalg::{field::{F2, Field}, flat_matrix::FlatMatrix, matrix::{RModMorphism, SmithNormalForm}, ring::{CRing, UniPolRing, ValuationRing}}, tensor::Tensor
};


#[allow(type_alias_bounds)]
pub type Module<B: BasisElement> = Vec<(B, UniGrading, Option<u16>)>;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Default)]
pub struct GradedModule<G: Grading, B: BasisElement>(pub HashMap<G, Module<B>>);



impl<G: Grading, B: BasisElement> GradedModule<G, B> {
    pub fn dimensions(&self) -> HashMap<G, usize> {
        self.0.iter().map(|(g, v)| (*g,v.len())).collect()
    }

    pub fn dimension_in_grade(&self, grade: &G) -> usize {
        self.0.get(grade).map(|x| x.len()).unwrap_or(0)
    }

    pub fn generate_tensor_as_module(&self, coalgebra: &Self, tensor: &Tensor<G>) -> GradedModule<G, B> {
        let mut map = HashMap::default();
        for (&gr, dim) in &tensor.dimensions {
            map.insert(gr, vec![(B::default(), UniGrading(0), None); *dim]);
        }
        for (&(t_gr, t_id), &((a_gr, a_id),(m_gr, m_id))) in &tensor.deconstruct {
            let module = &mut map.get_mut(&t_gr).unwrap()[t_id];
            let a_unigrade = coalgebra.0.get(&a_gr).unwrap()[a_id].1;
            let (_, m_unigrade, m_quotient) = self.0.get(&m_gr).unwrap()[m_id];
            module.1 = a_unigrade + m_unigrade;
            module.2 = m_quotient;
        }

        GradedModule(map)
    }
}



fn reduce<F: Field, B: BasisElement>(map: &mut Map<F>, codomain: &Module<B>) {
    for (y, (_,_,pow)) in codomain.iter().enumerate() {
        match pow {
            Some(pow) => {
                for x in 0..map.domain {
                    if map.get(x, y).1 >= *pow {
                        map.set(x, y, UniPolRing::zero());
                    }
                }
            },
            None => {},
        }
    }
}


type Map<F> = FlatMatrix<UniPolRing<F>>;

pub fn order_maps<F: Field, B: BasisElement>(f: &Map<F>, g: &Map<F>, n: &Module<B>, q: &Module<B>) -> (Map<F>, Map<F>, Module<B>, Module<B>, FlatMatrix<UniPolRing<F>>) {
    let n_sorted: Vec<_> = n.iter().enumerate().sorted_by(|a, b|{
        match a.1.2 {
            Some(pow) => {
                match b.1.2 {
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
                match b.1.2 {
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
        match a.1.2 {
            Some(pow) => {
                match b.1.2 {
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
                match b.1.2 {
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

fn multiply_vector<R: CRing>(vec: Vec<R>, el: R) -> Vec<R> {
    vec.into_iter().map(|x| el * x).collect()
}

/// M -f> N -g> Q
pub fn cohomology<F: Field, B: BasisElement>(f: &Map<F>, g: &Map<F>, n: &Module<B>, q: &Module<B>) -> (Module<B>, FlatMatrix<UniPolRing<F>>) {
    debug_assert_eq!(f.codomain, g.domain);

    if cfg!(debug_assertions) {
        let mut comp = g.compose(&f);
        reduce(&mut comp, &q);
        for x in 0..comp.domain {
            for y in 0..comp.codomain {
                assert!(comp.get(x, y).is_zero());
            }
        }
    }

    println!("f:\n{:?}g:\n{:?}n: {:?}\nq: {:?}\n", f, g, n,q);

    let (new_f,new_g, new_n, new_q, trans_map_inv) = order_maps(f, g, n, q);
    let (f,g,n,q) = (&new_f, &new_g, &new_n, &new_q);
    
    println!("f:\n{:?}g:\n{:?}n: {:?}\nq: {:?}\n", f, g, n,q);


    let new_domain = g.domain + g.codomain;
    let mut g_aug = FlatMatrix::zero(new_domain, g.codomain);

    for i in 0..g_aug.codomain {
        g_aug.set_row(i, g.get_row(i));
        
        
        if let Some(power) = q[i].2 {
            g_aug.set(g.domain + i, i, UniPolRing(F::one().neg(),power));
        }
    }


    let (_,s_aug,v_aug) = g_aug.snf();

    let mut start_zeros = s_aug.codomain;
    // Assume that s is ordered in some way (nonzero first)
    for r in 0..s_aug.codomain {
        let el = s_aug.get(r, r);
        if el.is_zero() {
            start_zeros = r;
            break;
        }
    }

    println!("g_aug:\n{:?}s_aug:\n{:?}v_aug:\n{:?}", g_aug, s_aug, v_aug);

    let g_ker_size = s_aug.domain - start_zeros;
    let mut g_ker = FlatMatrix::zero(g_ker_size, g.domain);

    
    // g_ker has as columns the vectors which generate the kernel
    for r in 0..g.domain {
        let row = &v_aug.get_row(r)[start_zeros..s_aug.domain];
        g_ker.set_row(r, row);
    }

    
    // TODO : THIS IS OPTIoNAL ??
    // reduce(&mut g_ker, &n);
    // println!("g_ker:\n{:?}", g_ker);
    
    
    let (_,s_ker,v_ker, uinv_ker, _) = g_ker.full_snf();
    
    let mut non_zero_els = 0;
    for r in 0..s_ker.codomain {
        let el = s_ker.get(r, r);
        if el.is_zero() {
            break;
        }
        non_zero_els += 1;
    }

    println!("s_ker:\n{:?}uinv_ker:\n{:?}v_ker:\n{:?}", s_ker, uinv_ker, v_ker);
    

    let mut vecs = vec![];
    for r in 0..non_zero_els {
        let vector = uinv_ker.get_column(r);
        vecs.push(multiply_vector(vector, s_ker.get(r, r)));
    }



    let mut real_ker_module_structure = vec![];
    vecs = vecs.into_iter().filter(|v| {
        let mut mod_str = Some(0);
        let mut mod_gr = UniGrading(0);
        for y in 0..v.len() {
            let el = v[y];
            if el.is_zero() {
                continue;
            }

            let gr = n[y].1;
            let n_x_structure = n[y].2;
            mod_gr = gr - UniGrading(el.1 as i32);
            match n_x_structure {
                Some(n_x_power) => {
                    let new_power = n_x_power - el.1; // If non-negative then the map was unreduced, probably ?

                    if let Some(p) = mod_str {
                        if p < new_power {
                            mod_str = Some(new_power);
                        }
                    }
                },
                None => {
                    mod_str = None;
                    break;
                },
            }
        }
        
        if mod_str == Some(0) {
            false 
        } else {
            real_ker_module_structure.push((kBasisElement::default(), mod_gr, mod_str));
            true
        }
    }).collect();

    let mut real_g_ker = FlatMatrix::zero(vecs.len(), g.domain);
    
    
    for (id, column) in vecs.into_iter().enumerate() {
        real_g_ker.set_column(id, &column[..]);
    }

    // TODO : OPTIONAL ?
    // reduce(&mut real_g_ker, n);
    println!("real_g_ker:\n{:?}", real_g_ker);


    

    // Check if vectors are REALLY zero
    if cfg!(debug_assertions) {
        for a in 0..real_g_ker.domain {
            let column = real_g_ker.get_column(a);
            let eval = g.eval_vector(&column);
            
            for (id,r) in eval.iter().enumerate() {
                let power = match q[id].2 {
                    Some(p) => p,
                    None => u16::MAX,
                };
                if !(r.is_zero() || r.1 >= power) {
                    panic!("OH OH, kernel is NOT zero")
                }
            }

        }
    }

    println!("real_ker_module: {:?}\n\n", real_ker_module_structure);
    
    let (u_real_ker,s_real_ker,_) = real_g_ker.snf();
    let mut Y = u_real_ker.compose(&f);
    
    // reduce(Y, codomain);
    
    println!("u_real_ker:\n{:?}s_real_ker:\n{:?}Y:\n{:?}", u_real_ker,s_real_ker, Y);

    let mut sol = FlatMatrix::zero(f.domain + real_ker_module_structure.len(), real_g_ker.domain);
    for x in 0..f.domain {
        for y in 0..sol.codomain {
            let y_i = Y.get(x, y);
            if y_i.is_zero() {
                continue;
            }
            let d_i = s_real_ker.get(y, y);
            let z_i = y_i.unsafe_divide(d_i);
            sol.set(x, y, z_i);
        }
    }

    for (r,el) in real_ker_module_structure.iter().enumerate() {
        match el.2 {
            Some(power) => {
                sol.set(f.domain + r, r, UniPolRing(F::one().neg(), power));
                
            },
            None => {},
        }
    }

    
    let (_,sol_s,_, sol_uinv, _) = sol.full_snf();

    println!("sol_s:\n{:?}sol_uinv:\n{:?}sol:\n{:?}", sol_s, sol_uinv, sol);

    let mut module = Module::default();
    let mut columns = vec![];
    for r in 0..sol_s.codomain {
        let el = sol_s.get(r, r);
        if !el.is_unit() {
            let structure = {
                if el.is_zero() { None } else { Some(el.1) }
            };
            let col = sol_uinv.get_column(r);

            let mut grade = UniGrading::infty();
            for c in 0..col.len() {
                if !col[c].is_zero() {
                    let gr = real_ker_module_structure[c].1;
                    grade = gr - (UniGrading(col[c].1 as i32));
                    break;
                }
            }

            if cfg!(debug_assertions) {
                if grade == UniGrading::infty() {
                    panic!("THIS CAN NOT HAPPEN :(")
                }
            }

            module.push((B::default(), grade, structure));
            columns.push(col);
        }
    }
    
    let mut cohom_to_ker = FlatMatrix::zero(module.len(), real_g_ker.domain);
    for (id, column) in columns.into_iter().enumerate() {
        cohom_to_ker.set_column(id, &column[..]);
    }

    let ker_to_n = trans_map_inv.compose(&real_g_ker);

    println!("ker_to_n:\n{:?}\n\nker_to_cohom:\n{:?}", ker_to_n, cohom_to_ker);

    (module, ker_to_n.compose(&cohom_to_ker))
}
