
use std::{cmp::Ordering, ops::Neg};

use ahash::HashMap;
use itertools::Itertools;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use crate::{
    basiselement::BasisElement,
    grading::{Grading, UniGrading},
    linalg::{
        field::Field, flat_matrix::FlatMatrix, graded::BasisIndex, matrix::{RModMorphism, SmithNormalForm}, ring::{CRing, UniPolRing}
    }, module::{module::GradedModule}
};

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Default)]
pub struct GradedModuleMap<G: Grading, F: Field> {
    // (Grade, (BasisIndex, t multiplication))
    pub maps: HashMap<G, FlatMatrix<UniPolRing<F>>>,
}


impl<G: Grading, F: Field> GradedModuleMap<G, F> {    
    pub fn zero_codomain<B: BasisElement>(domain: &GradedModule<G,B>) -> Self {
        let maps = domain.0.iter().map(|(g,v)| {
            (*g,FlatMatrix::zero(v.len(), 0))
        }).collect();

        Self {
            maps
        }
    }

    pub fn vstack(&mut self, other: &mut Self) {
        self.maps
            .iter_mut()
            .for_each(|(grade, self_mat)| match other.maps.get_mut(grade) {
                Some(other_mat) => {
                    self_mat.vstack(other_mat);
                }
                None => {}
            });
        other.maps.drain().for_each(|(g, map)| {
            if !self.maps.contains_key(&g) {
                self.maps.insert(g, map);
            }
        });
    }

    /// Note that calling this function invalidates the domains / codomains it references unless these are also summed
    pub fn block_sum(&mut self, other: &mut Self) {
        self.maps
            .iter_mut()
            .for_each(|(grade, self_mat)| match other.maps.get_mut(grade) {
                Some(other_mat) => {
                    self_mat.block_sum(other_mat);
                }
                None => {}
            });
        other.maps.drain().for_each(|(g, map)| {
            if !self.maps.contains_key(&g) {
                self.maps.insert(g, map);
            }
        });
    }

    pub fn verify<B: BasisElement>(&self, domain: &GradedModule<G,B>, codomain: &GradedModule<G,B>) -> Result<(), String> {
        for (grade, map) in &self.maps {
            if map.domain != domain.dimension_in_grade(grade) {
                return Err(format!("In {grade}, the domain size of the map does not equal that of the module"));
            }
            if map.codomain != codomain.dimension_in_grade(grade) {
                return Err(format!("In {grade}, the codomain size of the map does not equal that of the module"));
            }

            for dom in 0..map.domain {
                for codom in 0..map.codomain {
                    let el = map.get(dom, codom);
                    if !el.is_zero() {
                        let power = el.1;
                        let (_, dom_grade, dom_module) = domain.0.get(&grade).unwrap()[dom];
                        let (_, codom_grade, codom_module) = codomain.0.get(&grade).unwrap()[codom];

                        if let Some(codom_module_power) = codom_module {
                            if power >= codom_module_power {
                                return Err(format!("There is an element mapping to some power of a codomain, but the codomain has quotient which divides out this power. Domain {dom} id, codomain {codom} id, power {power}."));
                            }
                        }

                        if dom_grade + (UniGrading(power as i32)) != codom_grade  {
                            return Err(format!("Grade of domain does not map to correct grade of codomain. Domain {dom_grade}, codomain {codom_grade}, power {power}."));
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
                                        return Err(format!("Cannot have a (non-zero) map from a cyclic module to a free module. In grade {grade}."));
                                    },
                                }
                            },
                            None => {},
                        }
                    }
                }
            }
        }
        Ok(())
    }

   

    // pub fn reduce(&self) -> Self {
    //     let mut total_map = HashMap::default();
    //     let mut total_expl = HashMap::default();

    //     for (g, map) in &self.maps {
    //         let expl = self.explained.get(g).expect("The codomain should have been explained for the map in this grade");

    //         let mut new_codoms = vec![];
    //         for codom in 0..map.codomain {
    //             let zero = map.get_row(codom).iter().fold(true, |b, r| r.is_zero() && b);
    //             if !zero {
    //                 new_codoms.push(codom);
    //             }
    //         }

    //         let mut new_map = FlatMatrix::zero(map.domain, new_codoms.len());
    //         let mut new_expl = vec![];

    //         for (index, codom) in new_codoms.iter().enumerate() {
    //             let row = map.get_row(*codom);
    //             new_map.set_row(index, row);

    //             let exp = expl[*codom];
    //             new_expl.push(exp);
    //         }

    //         total_map.insert(*g, new_map);
    //         total_expl.insert(*g, new_expl);
    //     }

    //     Self {
    //         maps: total_map,
    //         explained: total_expl
    //     }
    // }

    pub fn compose(&self, rhs: &Self) -> Self {
        let mut compose: HashMap<G, FlatMatrix<UniPolRing<F>>> = self
            .maps
            .par_iter()
            .filter_map(|(k, v)| match rhs.maps.get(&k) {
                None => None,
                Some(t) => Some((*k, v.compose(t))),
            })
            .collect();

        for (self_gr, val) in self.maps.iter() {
            if !compose.contains_key(self_gr) {
                compose.insert(*self_gr, FlatMatrix::zero(0, val.codomain()));
            }
        }

        for (rhs_gr, val) in rhs.maps.iter() {
            if !compose.contains_key(rhs_gr) {
                compose.insert(*rhs_gr, FlatMatrix::zero(val.domain(), 0));
            }
        }

        // TODO ! check if powers are now too high
        
        Self { 
            maps: compose 
        }
    }
} 


impl<G: Grading, F: Field> GradedModuleMap<G, F> {
    pub fn zero<B: BasisElement>(domain: &GradedModule<G, B>, codomain: &GradedModule<G, B>) -> Self {
        let mut maps: HashMap<G, FlatMatrix<UniPolRing<F>>> = domain
            .0
            .iter()
            .map(|(g, els)| {
                let codom_len = codomain.dimension_in_grade(g);
                (*g, RModMorphism::zero(els.len(), codom_len))
            })
            .collect();
        codomain.0.iter().for_each(|(g, v)| {
            if !maps.contains_key(g) {
                maps.insert(*g, RModMorphism::zero(0, v.len()));
            }
        });
        Self {
            maps,
        }
    }

    /// If self is a map from A -> B, let Q be the cokernel. 
    /// Then we return the map from (B -> Q, Q -> B, Q). 
    /// "The map from B to Q, an 'inverse' map from Q to B and the cokernel object Q. 
    pub fn cokernel<B: BasisElement>(&self, codomain: &GradedModule<G,B>) -> (Self, Self, GradedModule<G, B>) {        
        if cfg!(debug_assertions) {
            for g in codomain.0.keys() {
                assert!(self.maps.contains_key(g));
            }
        }
        
        let mut coker = HashMap::default();
        let mut coker_map = HashMap::default();        
        let mut coker_inv_map = HashMap::default(); 

        for (&g, mat) in &self.maps { // TODO, parallelization

            let codom_g = match codomain.0.get(&g) {
                Some(e) => e,
                None => {
                    continue;
                },
            };

            let sorted: Vec<_> = codom_g.iter().enumerate().sorted_by(|a, b|{
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
            }).map(|x| (x.0, {
                match x.1.2 {
                    Some(pow) => {
                        UniPolRing(F::one().neg(), pow)
                    },
                    None => {
                        UniPolRing::zero()
                    },
                }})).collect();


            let mut trans_map = FlatMatrix::zero(mat.codomain, mat.codomain);    
            let mut trans_map_inv = FlatMatrix::zero(mat.codomain, mat.codomain);    
            
            for (target, (origin, _)) in sorted.iter().enumerate() {
                trans_map.set(*origin, target, UniPolRing::one());
                trans_map_inv.set(target, *origin, UniPolRing::one());
            }

            let mut mat_relations = FlatMatrix::zero(mat.domain + codom_g.len(), mat.codomain);

            // println!("{:?}", sorted);

            // TODO : Optimize by considering less relations ? (i.e. the free k[t] modules)
            for (target, (origin, el)) in sorted.clone().into_iter().enumerate() {
                
                mat_relations.set(mat.domain + target, target, el);
                
                let row = mat.get_row(origin);
                mat_relations.set_row(target, row);
            }

             
            let (u, s, _, uinv, _) = mat_relations.full_snf();


            
            let mut module: Vec<_> = (0..s.codomain).filter_map(|r| {
                let el = s.get(r, r);
                
                if el.is_unit() { None }
                else if !el.is_zero() {
                    // The grading is set correctly later !
                    Some((B::default(), UniGrading(0), Some(el.1)))
                } else { 
                    // The grading is set correctly later !
                    Some((B::default(), UniGrading(0), None))
                }
            }).collect();

            // println!("{:?}", g);
            // println!("{:?}", trans_map);
            // let v: Vec<_> = codom_g.iter().enumerate().map(|(id,x)| (id, x.2)).collect();
            // println!("{:?}", sorted);
            // let v_alt: Vec<_> = module.iter().enumerate().map(|(id,x)| (id, x.2)).collect();
            // println!("{:?}", v_alt);
            // println!("{:?}{:?}{:?}{:?}", mat, mat_relations, u, s);
        

            let mut g_map = FlatMatrix::zero(s.codomain, module.len());
            let mut g_inv_map = FlatMatrix::zero(module.len(), s.codomain);

            let diff = mat.codomain - module.len();

            for coker_id in 0..module.len() {

                let row = u.get_row(coker_id + diff);
                g_map.set_row(coker_id, &row);
                
                for (id,el) in row.iter().enumerate() {
                    if el.is_unit() {
                        // This is incorrect ??
                        module[coker_id].1 = codomain.0.get(&g).unwrap()[id].1; 
                    }
                }

                // This is correct!
                // Note that our SNF is ordered, thus all units appear in the topleft
                // Thusss we just truncate the row
                let uinv_column = &uinv.get_column(coker_id + diff)[0..s.codomain];
                g_inv_map.set_column(coker_id, &uinv_column);
            }


            let mut g_map = g_map.compose(&trans_map);
            let mut g_inv_map = trans_map_inv.compose(&g_inv_map);

            // Reduce g_map, (if something maps to t^k which is zero in some thing, set it to zero)
            for coker_id in 0..g_map.codomain {
                if let Some(power) = module[coker_id].2 {
                    for codom_id in 0..g_map.domain {
                        let el = g_map.get(codom_id, coker_id);
                        if el.1 >= power && !el.is_zero() {
                            g_map.set(codom_id, coker_id, UniPolRing::zero());
                        }
                    }
                }
            }

            // Reduce g_map, (if something maps to t^k which is zero in some thing, set it to zero)
            for codom_id in 0..g_inv_map.codomain {
                if let Some(power) = codom_g[codom_id].2 {
                    for coker_id in 0..g_inv_map.domain {
                        let el = g_inv_map.get(coker_id, codom_id);
                        if el.1 >= power && !el.is_zero() {
                            g_inv_map.set(coker_id, codom_id, UniPolRing::zero());
                        }
                    }
                }
            }


            for coker_id in 0..module.len() {
                for y in 0..g_map.domain {
                    if g_map.get(y, coker_id).is_unit() {
                        module[coker_id].1 = codomain.0.get(&g).unwrap()[y].1; 
                    }
                }
            }

            // println!("{:?}", codom_g);
            // println!("{:?}", module);

            // let mut divide_in_chunks = vec![];
            // let mut count = 0;
            // let mut latest = Some(0);
            // for el in &module {
            //     if el.2 == latest {
            //         count += 1;
            //     } else {
            //         divide_in_chunks.push(count);
            //         count = 1;
            //         latest = el.2
            //     }
            // }
            // divide_in_chunks.push(count);


            // println!("{:?}", divide_in_chunks);
            // println!("{:?}", g_map);

            // let mut total = 0;
            // for a in divide_in_chunks {
            //     if a > 1 {
            //         let mut possible = vec![true; g_map.domain];
            //         'outer : for codom in 0..g_map.domain {
            //             for coker_id in 0..total {
            //                 if !g_map.get(codom, coker_id).is_zero() {
            //                     possible[codom] = false;
            //                     continue 'outer; 
            //                 }
            //             }
                        
            //             for coker_id in total+a..g_map.codomain {
            //                 if !g_map.get(codom, coker_id).is_zero() {
            //                     possible[codom] = false;
            //                     continue 'outer; 
            //                 }
            //             }
            //         }
            //         let count = possible.iter().fold(0, |x, y| if *y { x + 1} else { x });
            //         println!("{:?}", possible);
            //         println!("{count}");
            //         let k: Vec<_> = module.iter().map(|x| x.2).enumerate().collect();
            //         println!("{:?}", k);
            //         println!("{:?}", g_map);
                    
            //         let mut fixed_pivots = vec![];
            //         if count < a { 
            //             println!("OH NO, NO POSSIBLE PIVOTS HERE :("); 
            //         };

            //         for (codom_id, a) in possible.iter().enumerate() {
            //             if *a {
            //                 for coker_id in 0..g_map.codomain {
            //                     if !fixed_pivots.contains(&coker_id) {
            //                         let el = g_map.get(codom_id, coker_id);
            //                         if let Some(inv) = el.try_inverse() {
            //                             g_map.multiply_row(coker_id, inv);
            //                             // This iterator can be with less items ?
            //                             for alt_coker_id in (0..coker_id).chain((coker_id+1)..g_map.codomain) {
            //                                 // Should check if ot_el is invertible ?
            //                                 let ot_el = g_map.get(codom_id, alt_coker_id);
            //                                 if !ot_el.is_zero() {
            //                                     g_map.add_row_multiple(alt_coker_id, coker_id, ot_el.neg());
            //                                 }
            //                             }
            //                             fixed_pivots.push(coker_id);
            //                             break;
            //                         }
            //                     }
            //                 }
                            
            //                 // println!("{:?}", g_map);
            //             }
            //         }
            //         if count < a { 
            //             println!("{:?}", g_map);

            //             panic!("OH NO, NO POSSIBLE PIVOTS HERE :(") };

            //     }
            //     total += a;
            // }


            // println!("{:?}", g);
            

            // TODO: MAKE SURE THE PIVOTS ARE ACTUALLY 1 INSTEAD OF JUST A UNIT
            // let g_inv_map = g_map.extensive_pivots();
            // let mut g_inv_map = FlatMatrix::zero(g_map.codomain, g_map.domain);


            // for (coker_id, &codom_id) in pivots.iter().enumerate() {
            //     g_inv_map.set(coker_id, codom_id, UniPolRing::one());
            // }


            // println!("{:?}{:?}", g_map, g_inv_map);
            
            if cfg!(debug_assertions) {
                let mut comp = g_map.compose(&g_inv_map);
                for (y, (_,_,power)) in module.iter().enumerate() {
                    match power {
                        Some(power) => {
                            for x in 0..comp.domain() {
                                let el = comp.get(x, y);
                                if !el.is_zero() {
                                    if el.1 >= *power {
                                        comp.set(x,y,UniPolRing::zero());
                                    }
                                } 
                            }
                        },
                        None => {},
                    }
                }
                comp.is_unit().unwrap();
            }

            if module.len() > 0 {
                coker.insert(g, module);
            }

            coker_map.insert(g, g_map);
            coker_inv_map.insert(g, g_inv_map);
        }

        let to = GradedModuleMap { maps: coker_map };
        let inv = GradedModuleMap { maps: coker_inv_map };
        let coker = GradedModule(coker);

        if cfg!(debug_assertions) {
            to.verify(codomain, &coker).unwrap();
        }
        
        (to, inv, coker)
    }


    /// Let F: M -> N and g a grade 
    /// this function returns a map k: ker(F_g) -> M_g and a module ker(F)
    pub fn basis_element_kernel_pivot_in_grade<B: BasisElement>(&self, domain: &GradedModule<G, B>, codomain: &GradedModule<G,B>, grade: G) -> Option<BasisIndex<G>> {
        let map = self.maps.get(&grade).unwrap();        
        let real_domain = domain.0.get(&grade).unwrap();

        if map.codomain == 0 {
            if map.domain == 0 {
                return None;
            }
            else {
                return Some((grade, 0));
            }
        }    
    
        let new_domain = map.domain + map.codomain;
        let mut new_map = FlatMatrix::zero(new_domain, map.codomain);
        
        for n in  0..new_map.codomain {
            new_map.set_row(n, map.get_row(n));
            
            
            if let Some(power) = codomain.0.get(&grade).unwrap()[n].2 {
                new_map.set(map.domain + n, n, UniPolRing(F::one().neg(),power));
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
                for y in 0..map.domain {
                    let domain_pow = real_domain[y].2;
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
            for y in 0..map.domain {
                let domain_pow = real_domain[y].2;
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

        match possible_id {
            Some(id) => {
                Some((grade, id))
            },
            None => None,
        }
    }
}