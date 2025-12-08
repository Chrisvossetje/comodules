use std::cmp::Ordering;

use itertools::Itertools;

use crate::{field::Field, matrices::flat_matrix::FlatMatrix, matrix::Matrix, ring::CRing, rings::univariate_polynomial_ring::UniPolRing, snf::SmithNormalForm, unipol::UniPolModule};

/// If self is a map from A -> B, let Q be the cokernel. 
/// Then we return the map from (B -> Q, Q -> B, Q). 
/// "The map from B to Q, representing vectors of Q in B and the cokernel object Q. 

impl<F: Field> FlatMatrix<UniPolRing<F>> {
    pub(super) fn internal_cokernel(&self, codomain: &UniPolModule) -> (Self, Self, UniPolModule) {  
        // TODO : Force maps to always be sorted      
        let sorted: Vec<_> = codomain.iter().enumerate().sorted_by(|a, b|{
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
        }).map(|x| (x.0, {
            match x.1 {
                Some(pow) => {
                    UniPolRing(F::one().neg(), *pow)
                },
                None => {
                    UniPolRing::zero()
                },
            }})).collect();


        let mut trans_map = FlatMatrix::zero(self.codomain, self.codomain);    
        let mut trans_map_inv = FlatMatrix::zero(self.codomain, self.codomain);    
        
        for (target, (origin, _)) in sorted.iter().enumerate() {
            trans_map.set(*origin, target, UniPolRing::one());
            trans_map_inv.set(target, *origin, UniPolRing::one());
        }

        let mut mat_relations = FlatMatrix::zero(self.domain + codomain.len(), self.codomain);

        // println!("{:?}", sorted);

        // TODO : Optimize by considering less relations ? (i.e. the free k[t] modules)
        for (target, (origin, el)) in sorted.clone().into_iter().enumerate() {
            
            mat_relations.set(self.domain + target, target, el);
            
            let row = self.get_row(origin);
            mat_relations.set_row(target, row);
        }

            
        let (u, s, _, uinv, _) = mat_relations.full_snf();


        
        // TODO : As we force maps to already be sorted, this must change a litlle bit maybs ?
        let module: Vec<_> = (0..s.codomain).filter_map(|r| {
            let el = s.get(r, r);
            
            if el.is_unit() { None }
            else if !el.is_zero() {
                Some(Some(el.1))
            } else { 
                Some(None)
            }
        }).collect();

        // println!("{:?}", g);
        // println!("{:?}", trans_map);
        // let v: Vec<_> = codomain.iter().enumerate().map(|(id,x)| (id, x.2)).collect();
        // println!("{:?}", sorted);
        // let v_alt: Vec<_> = module.iter().enumerate().map(|(id,x)| (id, x.2)).collect();
        // println!("{:?}", v_alt);
        // println!("{:?}{:?}{:?}{:?}", self, mat_relations, u, s);
    

        let mut g_map = FlatMatrix::zero(s.codomain, module.len());
        let mut g_inv_map = FlatMatrix::zero(module.len(), s.codomain);

        let diff = self.codomain - module.len();


        for coker_id in 0..module.len() {

            let row = u.get_row(coker_id + diff);
            g_map.set_row(coker_id, &row);
            
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
            if let Some(power) = module[coker_id] {
                for codom_id in 0..g_map.domain {
                    let el = g_map.get(codom_id, coker_id);
                    if el.1 >= power && !el.is_zero() {
                        g_map.set(codom_id, coker_id, UniPolRing::zero());
                    }
                }
            }
        }

        // Reduce g_inv_map, (if something maps to t^k which is zero in some thing, set it to zero)
        for codom_id in 0..g_inv_map.codomain {
            if let Some(power) = codomain[codom_id] {
                for coker_id in 0..g_inv_map.domain {
                    let el = g_inv_map.get(coker_id, codom_id);
                    if el.1 >= power && !el.is_zero() {
                        g_inv_map.set(coker_id, codom_id, UniPolRing::zero());
                    }
                }
            }
        }

        // println!("{:?}", codomain);
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
            for (y, power) in module.iter().enumerate() {
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

        (g_map, g_inv_map, module)
    }
}