use ahash::{HashMap, HashMapExt};
use serde::{Deserialize, Serialize};

use crate::linalg::{field::Field, flat_matrix::FlatMatrix, graded::{BasisElement, BasisIndex}, grading::{BiGrading, Grading}, matrix::{RModMorphism, SmithNormalForm}, ring::{CRing, UniPolRing, ValuationRing}};



pub fn print_matrix<F: Field>(mat: &FlatMatrix<UniPolRing<F>>) {
    for m in 0..mat.codomain {
        for n in 0..mat.domain {
            if mat.get(n, m).0.is_zero() {
                print!("0     ");
            } else {
                if mat.get(n, m).1 == 0 {
                    print!("1     ");

                } else {
                    print!("t^{:?}   ", mat.get(n, m).1);
                }
            }
        }
        println!();
    }
    println!();

}
    

// impl<F: Field> SmithNormalForm<UniPolRing<F>> for FlatMatrix<UniPolRing<F>> {
impl<R: ValuationRing> SmithNormalForm<R> for FlatMatrix<R> {
    fn full_snf(&self) -> (Self, Self, Self, Self, Self) {

        enum Action<R: ValuationRing> {
            Swap(usize, usize),
            Add(usize, usize, R)
        }

        
        let mut s = self.clone();
        let mut u = Self::identity(self.codomain());
        let mut v = Self::identity(self.domain());

        let mut u_inv_actions = vec![];
        let mut v_inv_actions = vec![];
        
        let m = s.codomain();
        let n = s.domain();
        
        let min_r = m.min(n);

        for r in 0..min_r {

            // Find element in (sub)matrix which divides all others
            let mut candidate = (r,r);
            let mut candidate_val = s.get(r, r); 
            for x in r..n {
                for y in r..m {
                    let el = s.get(x, y);
                    if !candidate_val.divides(&el) {
                        candidate_val = el;
                        candidate = (x,y);
                    } 
                }
            }

            // swap rows and colums to get candidate in correct position
            u.swap_rows(r, candidate.1);
            u_inv_actions.push(Action::Swap(r, candidate.1));
            s.swap_rows(r, candidate.1);

            s.swap_cols(r, candidate.0);
            v.swap_cols(r, candidate.0);
            v_inv_actions.push(Action::Swap(r, candidate.0));


            let pivot = candidate_val;

            // reduce rows below (r,r)
            for k in (r + 1)..m {
                let entry = s.get(r, k);
                if !entry.is_zero() {
                    let factor = -entry.unsafe_divide(pivot);

                    u.add_row_multiple(k, r, factor);
                    u_inv_actions.push(Action::Add(k,r,factor));
                    s.add_row_multiple(k, r, factor);
                }
            }

            // reduce columns next to (r,r)
            for l in (r + 1)..n {
                let entry = s.get(l, r);
                if !entry.is_zero() {
                    let factor = -entry.unsafe_divide(pivot);
                    s.add_col_multiple(l, r, factor);
                    v.add_col_multiple(l, r, factor);
                    v_inv_actions.push(Action::Add(l,r,factor));
                }
            }
        }

        
        // construct u and v inverse
        let mut u_inv = Self::identity(self.codomain());
        let mut v_inv = Self::identity(self.domain());
        
        for a in u_inv_actions.iter().rev() {
            match a {
                Action::Swap(b, c) => {
                    u_inv.swap_rows(*b, *c);
                },
                Action::Add(b, c, factor) => {
                    u_inv.add_row_multiple(*b, *c, factor.neg());
                },
            }
        }
        for a in v_inv_actions.iter().rev() {
            match a {
                Action::Swap(b, c) => {
                    v_inv.swap_cols(*b, *c);
                },
                Action::Add(b, c, factor) => {
                    v_inv.add_col_multiple(*b, *c, factor.neg());
                },
            }
        }

        (u, s, v, u_inv, v_inv)
    }
}

// TODO: THIS CAN BE REMOVED / IMPROVED TO WORK WITH PID's
// PROBABLY ONLY RELEVANT FOR THE PID Z (the integers) (or maybe a generic k[x])
// k[x,y] will require a different method (groebner :( )

// impl<F: Field> SmithNormalForm<UniPolRing<F>> for FlatMatrix<UniPolRing<F>> {    
//     fn snf(&self) -> (Self, Self, Self) 
//     where 
//         Self: Clone,
//     {
//         let mut s = self.clone();
//         let mut u = Self::identity(self.codomain());
//         let mut v = Self::identity(self.domain());
        
//         let m = s.codomain();
//         let n = s.domain();
//         // let mut column_indices = Vec::new();
//         let mut last_j = 0;
        
//         // Main algorithm: iterate t from 1 to m
//         for t in 0..m {
//             // println!("{}:::!!!\n",t);
//             // print_matrix(&s);

//             // Step I: Choosing a pivot
//             // Find smallest column index with non-zero entry starting from last_j            
//             let mut search = None;
//             'outer: for j in last_j..n { // column is choice of domain
//                 for a in t..m { // row is choice of codomain
//                     if !s.get(j, a).is_zero() {
//                         search = Some((a, j));
//                         break 'outer;
//                     }
//                 }
//             }

//             let (a, j_t) = match search {
//                 Some((a, j)) => (a, j),
//                 None => break, // No more non-zero entries
//             };
//             last_j = j_t + 1;
            
//             u.swap_rows(t, a);
//             s.swap_rows(t, a);

//             // As all "prime" factors are finite, this loop is finite
//             loop {
//                 // println!("Loop");
//                 // print_matrix(&s);

//                 let mut improved = false;
                
                
//                 // FIRST DO ROWS
                
//                 // Step II: Improving the pivot
//                 // Check if pivot divides all entries in column j_t
//                 let mut smallest_index = t;
//                 for k in (t + 1)..m { // rows
//                     if !s.get(j_t, smallest_index).divides(&s.get(j_t, k)) {
//                         smallest_index = k;
//                     }
//                 }

//                 u.swap_rows(t, smallest_index);
//                 s.swap_rows(t, smallest_index);
//                 if t != smallest_index {
//                     improved = true;
//                 }

//                 // println!("Row Step II");
//                 // print_matrix(&s);


//                 // Step III: Eliminating entries
//                 // Clear column j_t below position (t, j_t)
//                 let pivot = s.get(j_t, t);
//                 for k in (t + 1)..m {
//                     let entry = s.get(j_t, k);
//                     if !entry.is_zero() {
//                         let factor = -entry.unsafe_divide(pivot);
//                         // TODO
//                         u.add_row_multiple(k, t, factor);
//                         s.add_row_multiple(k, t, factor);
//                     }
//                 }
//                 // println!("Row Step III");
//                 // print_matrix(&s);


//                 // SECOND DO COLUMNS

//                 // Step II: Improving the pivot
//                 // Check if pivot divides all entries in column j_t
//                 let mut smallest_index = j_t;
//                 for k in (j_t + 1)..n {
//                     if !s.get(smallest_index, t).divides(&s.get(k, t)) {
//                         smallest_index = k;
//                     }
//                 }

                
//                 s.swap_cols(j_t, smallest_index);
//                 v.swap_cols(j_t, smallest_index);
//                 if j_t != smallest_index {
//                     improved = true;
//                 }
                
//                 // println!("Column Step II");
//                 // print_matrix(&s);


//                 // Step III: Eliminating entries
//                 // Clear column j_t below position (t, j_t)

//                 let pivot = s.get(j_t, t);
//                 for k in (j_t + 1)..n {
//                     let entry = s.get(k, t);
//                     if !entry.is_zero() {
//                         let factor = -entry.unsafe_divide(pivot);
//                         s.add_col_multiple(k, j_t, factor);
//                         v.add_col_multiple(k, j_t, factor);
//                     }
//                 }

//                 // println!("Column Step III");
//                 // print_matrix(&s);

//                 if !improved {
//                     break;
//                 }
//             }
            
            
//             // // Clear row t to the right of position (t, j_t) using column operations
//             // loop {
//             //     let mut cleared_row = true;
//             //     for j in (j_t + 1)..n {
//             //         let entry = s.get(j, t);
//             //         if !entry.is_zero() {
//             //             cleared_row = false;
                        
//             //             // Simple column elimination: subtract multiple of pivot column
//             //             let pivot = s.get(j_t, t);
//             //             if !pivot.is_zero() {
//             //                 // Use simple approach: subtract pivot column from current column
//             //                 s.add_col_multiple(j, j_t, -R::one());
//             //                 v.add_col_multiple(j, j_t, -R::one());
//             //             }
//             //             break;
//             //         }
//             //     }
                
//             //     if cleared_row {
//             //         break;
//             //     }
                
//             //     // Re-clear column if needed (entries may have become non-zero again)
//             //     let pivot = s.get(j_t, t);
//             //     for k in (t + 1)..m {
//             //         let entry = s.get(j_t, k);
//             //         if !entry.is_zero() && !pivot.is_zero() {
//             //             // Simple row elimination
//             //             s.add_row_multiple(k, t, -R::one());
//             //             u.add_row_multiple(k, t, -R::one());
//             //         }
//             //     }
//             // }
//         }
        
//         // Final step: ensure diagonal entries satisfy divisibility condition
//         // self.ensure_divisibility(&mut s, &mut u, &mut v);
        
//         (u, s, v)
//     }
// }



// #[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
// pub struct Module<B: BasisElement>(pub Vec<(B, Option<usize>)>);

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Default)]
pub struct GradedModule<G: Grading, B: BasisElement>(pub HashMap<G, Vec<(B, Option<usize>)>>);


impl<G: Grading, B: BasisElement> GradedModule<G, B> {
    pub fn dimension_in_grade(&self, grade: &G) -> usize {
        self.0.get(grade).map(|x| x.len()).unwrap_or(0)
    }
}

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Default)]
/// Note that the Grading is "nice" wrt the codomain, the domain has the more complex part
pub struct GradedModuleMap<G: Grading, F: Field> {
    // (Grade, (BasisIndex, t multiplication))
    // Do we want to assume that this vec is ordered by wrt t ?  Nah
    pub maps_domain: HashMap<G, Vec<(BasisIndex<G>, usize)>>, 
    pub maps_codomain: HashMap<G, Vec<(BasisIndex<G>, usize)>>,
    pub maps: HashMap<G, FlatMatrix<UniPolRing<F>>>,    
}

pub trait PolyGrading: Grading {
    fn gen_diff_power(source: Self, target: Self, gen: Self) -> Option<usize>;
}

impl PolyGrading for BiGrading {
    fn gen_diff_power(source: Self, target: Self, gen: Self) -> Option<usize> {
        let diff = target-source;
        let factor = {
            if diff.0 == 0 && gen.0 == 0 {
                if diff.1 % gen.1 == 0 {
                    if gen.1 == 0 {
                        if diff.1 == 0 {
                            Some(0) 
                        } else {
                            None
                        }
                    } else {
                        Some(diff.1 / gen.1)
                    }
                } else {
                    None
                }
            } else {
                if diff.0 % gen.0 == 0 {
                    if gen.0 == 0 {
                        if diff.0 == 0 {
                            Some(0) 
                        } else {
                            None
                        }
                    } else {
                        let factor = diff.0 / gen.0;
                        if gen.1 * factor == diff.1 {
                            Some(factor)
                        } else {
                            None
                        }
                    }
                } else {
                    None
                }
            }
        };


        match factor {
            Some(t) => {
                if t < 0 {
                    None
                } else {
                    Some(t as usize)
                }
            },
            None => None,
        }
    }
}

impl<G: PolyGrading,F: Field> GradedModuleMap<G, F> {
    pub fn zero<B: BasisElement>(domain: &GradedModule<G, B>, codomain: &GradedModule<G, B>, gen_grading: G) -> Self {
        let mut maps_domain = HashMap::new();
        for (g, _) in &codomain.0 {
            maps_domain.entry(*g).or_insert(vec![]);
            for (gg, mm) in &domain.0 {
                maps_domain.entry(*gg).or_insert(vec![]);
                if let Some(power) = G::gen_diff_power(*g, *gg, gen_grading) {
                    maps_domain.entry(*g).and_modify(|l| {
                        for (index, (_, quotient)) in mm.iter().enumerate() {
                            match quotient {
                                Some(quotient_power) => {
                                    // Here we check if the domain outscales 
                                    //(R/x^4) -x^5> codomain will always be zero !
                                    if *quotient_power > power {
                                        l.push(((*gg, index), power));
                                    }
                                },
                                None => {
                                    l.push(((*gg, index), power));
                                },
                            }

                        }
                    });
                }
            }
        }

        let maps_codomain = codomain.0.iter().map(|(g,v)| {
            (*g, (0..v.len()).map(|n| ((*g,n),0)).collect())
        }).collect();

        let mut maps = HashMap::new();
        for (g, l) in &maps_domain {
            let codom_size = codomain.0.get(&g).map_or(0, |e| e.len());
            let m = FlatMatrix::zero(l.len(), codom_size);
            maps.insert(*g, m);
        }

        Self {
            maps_domain,
            maps_codomain,
            maps,
        }
    }

    pub fn zero_codomain<B: BasisElement>(domain: &GradedModule<G,B>) -> Self {
        let maps_domain = domain.0.iter().map(|(g,v)| {
            (*g,(0..v.len()).map(|n| ((*g, n),0)).collect())
        }).collect();
        let maps_codomain = domain.0.iter().map(|(g, _)| {
            (*g,vec![])
        }).collect();
        let maps = domain.0.iter().map(|(g,v)| {
            (*g,FlatMatrix::zero(v.len(), 0))
        }).collect();
        
        Self {
            maps_domain,
            maps_codomain,
            maps,
        }
    }

    pub fn vstack(&mut self, other: &mut Self) {
        let (other_map, other_codom) = (&mut other.maps, &mut other.maps_codomain);

        self.maps
            .iter_mut()
            .for_each(|(grade, self_mat)| match other_map.get_mut(grade) {
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

        self.maps_codomain
            .iter_mut()
            .for_each(|(grade, self_list)| match other_codom.get_mut(grade) {
                Some(other_list) => {
                    self_list.append(other_list);
                }
                None => {}
            });
        other_codom.drain().for_each(|(g, map)| {
            if !self.maps_codomain.contains_key(&g) {
                self.maps_codomain.insert(g, map);
            }
        });
    }

    pub fn block_sum(&mut self, other: &mut Self) {
        let (other_map, other_dom, other_codom) = (&mut other.maps, &mut other.maps_domain, &mut other.maps_codomain);

        self.maps
            .iter_mut()
            .for_each(|(grade, self_mat)| match other_map.get(grade) {
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


        self.maps_domain
            .iter_mut()
            .for_each(|(grade, self_list)| match other_dom.get_mut(grade) {
                Some(other_list) => {
                    self_list.append(other_list);
                }
                None => {}
            });
        other_dom.drain().for_each(|(g, map)| {
            if !self.maps_domain.contains_key(&g) {
                self.maps_domain.insert(g, map);
            }
        });

        self.maps_codomain
            .iter_mut()
            .for_each(|(grade, self_list)| match other_codom.get_mut(grade) {
                Some(other_list) => {
                    self_list.append(other_list);
                }
                None => {}
            });
        other_codom.drain().for_each(|(g, map)| {
            if !self.maps_codomain.contains_key(&g) {
                self.maps_codomain.insert(g, map);
            }
        });


    }

    pub fn check_valid_quotients() {
        unimplemented!()
    }

    /// If self is a map from A -> B, let Q be the cokernel. 
    /// Then we return the map from (B -> Q, Q -> B, Q). 
    /// "The map from B to Q, an 'inverse' map from Q to B and the cokernel object Q. 
    pub fn cokernel<B: BasisElement>(&self) -> (Self, Self, GradedModule<G, B>) {
        let normalized_self = self.fix_codomain();

        let mut coker = HashMap::new();
        let mut coker_map = HashMap::new();
        let mut coker_map_domain = HashMap::new();
        let mut coker_map_codomain = HashMap::new();
        
        let mut coker_inv_map = HashMap::new(); 
        let mut coker_inv_domain = HashMap::new();
        let mut coker_inv_codomain = HashMap::new();

        for (&g, mat) in &normalized_self.maps {
            let (u, s, _, uinv, _) = mat.full_snf();
            
            let min_r = s.domain.min(s.codomain);
            let mut cokernel_module: Vec<Option<(B, Option<usize>)>> = (0..min_r).map(|r| {
                let el = s.get(r, r);
                if el.is_unit() {None}
                else if !el.is_zero() {
                    Some((B::default(), Some(el.1)))
                } else {
                    Some((B::default(), None))
                }
            }).collect();

            for _ in min_r..s.codomain {
                cokernel_module.push(Some((B::default(), None)));
            }

            let module: Vec<(B, Option<usize>)> = cokernel_module.clone().into_iter().filter_map(|e| e).collect();
            let mut g_map = FlatMatrix::zero(s.codomain, module.len());
            let mut g_inv_map = FlatMatrix::zero(module.len(), s.codomain);

            let mut count = 0;
            for index in 0..module.len() {
                if cokernel_module.get(index).is_some() {
                    g_map.set_row(count, u.get_row(index));
                    g_inv_map.set_row(count, uinv.get_row(index));
                    count += 1;
                }
            }
            
            let l = module.len();
            coker.insert(g, module);

            coker_map.insert(g, g_map);            
            coker_map_domain.insert(g, (0..s.codomain).map(|n| {
                ((g, n), 0)
            }).collect());
            coker_map_codomain.insert(g, (0..l).map(|n| {
                ((g, n), 0)
            }).collect());
            
            coker_inv_map.insert(g, g_inv_map);
            coker_inv_domain.insert(g, (0..l).map(|n| {
                ((g, n), 0)
            }).collect());
            coker_inv_codomain.insert(g, (0..s.codomain).map(|n| {
                ((g, n), 0)
            }).collect());
        }

        (GradedModuleMap {
            maps_domain: coker_map_domain,
            maps_codomain: coker_map_codomain,
            maps: coker_map,
        }, GradedModuleMap {
            maps_domain: coker_inv_domain,
            maps_codomain: coker_inv_codomain,
            maps: coker_inv_map
        },
        GradedModule(coker))
    }


    /// If self is a map from A -> B, let Q be the cokernel. 
    /// Then we return the map from (B -> Q, Q -> B, Q). 
    /// "The map from B to Q, an 'inverse' map from Q to B and the cokernel object Q. 
    pub fn kernel_in_grade<B: BasisElement>(&self, _grade: G) -> (Self, Self, GradedModule<G, B>) {
        todo!()
    }

    pub fn fix_domain(&self) -> Self {
        todo!()
    }

    pub fn fix_codomain(&self) -> Self {
        todo!()
    }

    pub fn compose(&self, other: &Self) -> Self {
        let lhs = self.fix_codomain();
        let rhs = other.fix_domain();

        let compose: HashMap<G, FlatMatrix<UniPolRing<F>>> = lhs.maps.into_iter().filter_map(|(g, mat)| {
            match rhs.maps.get(&g) {
                Some(rhs_map) => Some((g, mat.compose(rhs_map))),
                None => None,
            }
        }).collect();

        let domain = compose.iter().map(|(&g,_)| {
            (g, lhs.maps_domain[&g].clone())
        }).collect();
        let codomain = compose.iter().map(|(&g,_)| {
            (g, lhs.maps_codomain[&g].clone())
        }).collect();

        Self {
            maps_domain: domain,
            maps_codomain: codomain,
            maps: compose,
        }
    }
} 