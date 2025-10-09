use std::{sync::Arc};

use ahash::HashMap;
use itertools::Itertools;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::{comodule::{kcomodule::kBasisElement, rcomodule::RComodule, tensor::Tensor, traits::{Comodule, ComoduleMorphism}}, linalg::{field::Field, flat_matrix::FlatMatrix, graded::BasisIndex, grading::{OrderedGrading}, matrix::RModMorphism, module::{GradedModuleMap, PolyGrading}, ring::{CRing, UniPolRing}}};



#[derive(Debug, Clone)]
pub struct RComoduleMorphism<G: PolyGrading, F: Field> {
    pub domain: Arc<RComodule<G, F>>,
    pub codomain: Arc<RComodule<G, F>>,
    pub map: GradedModuleMap<G, F>, 
}




impl<G: PolyGrading, F: Field> ComoduleMorphism<G, RComodule<G,F>>
    for RComoduleMorphism<G, F>
{
    fn cokernel(&self) -> Self {
        let (coker_to, coker_inv, coker) = self.map.cokernel::<kBasisElement>();

        let coalg = self.codomain.coalgebra.as_ref();
        let tensor = Tensor::generate(&coalg.space, &coker);

        let m_lut: HashMap<BasisIndex<G>, Vec<(usize, UniPolRing<F>)>> = self
            .codomain
            .tensor
            .construct
            .par_iter()
            .map(|((f_gr, f_id), _)| {
                // Transfer a specific codomain grade and id (f_gr, f_id) to a list of elements which it maps to in the cokernel
                let v = (0..coker_to
                    .maps
                    .get(f_gr)
                    .map(|map| map.codomain())
                    .unwrap_or(0))
                    .filter_map(|q_index| {
                        let val = coker_to.maps.get(f_gr).unwrap().get(*f_id, q_index);
                        match val.is_zero() {
                            true => None,
                            false => Some((q_index, val)),
                        }
                    })
                    .collect();
                ((*f_gr, *f_id), v)
            })
            .collect();

        let maps = coker.0.par_iter().map(|(g,v)| {
            let g_tensor_dimen = tensor.get_dimension(g);
            let mut g_coaction = FlatMatrix::<UniPolRing::<F>>::zero(v.len(), g_tensor_dimen);

            let map = &coker_inv.maps[g];
            
            for coker_id in 0..map.domain {
                for codom_id in 0..map.codomain {
                    let inv_val = map.get(coker_id, codom_id);
                    if inv_val.is_zero() {
                        continue;
                    }

                    let coact_size = self.codomain.tensor.dimensions[g];
                    for codom_coact_id in 0..coact_size {
                        let coact_val =
                            self.codomain.coaction.maps[g].get(codom_id, codom_coact_id);
                        if !coact_val.is_zero() {
                            let ((alg_gr, alg_id), (mod_gr, mod_id)) =
                                self.codomain.tensor.deconstruct[&(*g, codom_coact_id)];

                            for (target_id, val) in m_lut.get(&(mod_gr, mod_id)).unwrap() {
                                let (_, final_id) =
                                    tensor.construct[&(mod_gr, *target_id)][&(alg_gr, alg_id)];
                                g_coaction.add_at(coker_id, final_id, coact_val * *val);
                            }
                        }
                    }
                }
            }
            (*g, g_coaction)

        }).collect();

        let domain = coker.0.iter().map(|(g,v)| {
            (*g, (0..v.len()).map(|n| {
                ((*g, n), n)
            }).collect())
        }).collect();
        let codomain = tensor.dimensions.iter().map(|(g,v)| {
            (*g, (0..*v).map(|n| {
                ((*g, n), n)
            }).collect())
        }).collect();

        let coaction = GradedModuleMap {
            maps_domain: domain,
            maps_codomain: codomain,
            maps,
        };

        let coker_comod = RComodule {
            coalgebra: self.codomain.coalgebra.clone(),
            space: coker,
            coaction,
            tensor,
        };

        Self {
            domain: self.codomain.clone(),
            codomain: Arc::new(coker_comod),
            map: coker_to,
        }
    }

    fn inject_codomain_to_cofree(&self, limit: G) -> Self 
    where
        G: OrderedGrading {
        
        let mut growing_map: GradedModuleMap<G, F> =
            GradedModuleMap::zero_codomain(&self.codomain.space);
        let mut growing_comodule = RComodule::zero_comodule(self.codomain.coalgebra.clone());
        let mut iteration = 0;

        let grades: Vec<G> = growing_map
            .maps
            .iter()
            .map(|(g, _)| *g)
            .sorted_by(|a,b| {
                a.compare(b)
            })  
            .filter(|&g| g.compare(&limit).is_le())
            .collect();
        let prev_grade = 0;

        let fixed_limit = limit.incr().incr();


        loop {
            // Get lowest graded pivot element
            let pivot = None;
            for grade_id in prev_grade..grades.len() {
                let _grade = grades[grade_id];
                todo!()
                // let kernel = growing_map.kernel_in_grade(&grade);
                
                // match kernel.first_non_zero_entry() {
                //     Some(loc) => {
                //         prev_grade = grade_id;
                //         pivot = Some((loc, grade));
                //         break;
                //     }
                //     None => {}
                // }
            }
            let (pivot, pivot_grade): ((usize, usize), G) = match pivot {
                Some(p) => p,
                _ => {
                    break;
                }
            };

            let alg_to_tens = self
                .codomain
                .tensor
                .construct
                .get(&(pivot_grade, pivot.1))
                .expect("The tensor should exist on the codomain in this grade");

            let coalg_space = &self.codomain.coalgebra.space;

            // TODO: Verify is this parallel iterator is faster or not for big(ger) coalgebras
            

            let cofree_map: HashMap<G,FlatMatrix<UniPolRing<F>>> = coalg_space.0.iter().filter_map(|(alg_gr, alg_gr_space)| {
                let t_gr = *alg_gr + pivot_grade;

                if t_gr > fixed_limit {
                    return None;
                }

                let codomain_len = self.codomain.space.dimension_in_grade(&t_gr);
                let coalg_len = alg_gr_space.len();

                if !alg_to_tens.contains_key(&(*alg_gr, 0)) {
                    let zero_map = FlatMatrix::zero(codomain_len, coalg_len);
                    return Some((t_gr, zero_map));
                };

                let mut map = FlatMatrix::zero(codomain_len, coalg_len);

                for a_id in 0..coalg_len {
                    let (t_gr, t_id) = alg_to_tens.get(&(*alg_gr,a_id)).expect("This BasisIndex should exist on the tensor object in the to inject comodule");
                    let slice = self
                        .codomain
                        .coaction
                        .maps
                        .get(t_gr)
                        .expect("This grade should exist on the coaction of the injecting comodule")
                        .get_row(*t_id);
                    map.set_row(a_id, slice);
                }

                Some((t_gr, map))
            }).collect();

            let mut f = RComodule::cofree_comodule(
                self.codomain.coalgebra.clone(),
                iteration,
                pivot_grade,
                fixed_limit,
            );

            let mut extend_map = GradedModuleMap {
                maps_domain: HashMap::default(),
                maps_codomain: HashMap::default(),
                maps: cofree_map,
            };

            growing_comodule.direct_sum(&mut f);
            growing_map.vstack(&mut extend_map);

            iteration += 1;
        }

        Self {
            domain: self.codomain.clone(),
            codomain: Arc::new(growing_comodule),
            map: growing_map,
        }
    }

    fn zero_morphism(comodule: Arc<RComodule<G,F>>) -> Self {
        let zero = Arc::new(RComodule::zero_comodule(comodule.coalgebra.clone()));

        let zero_map = GradedModuleMap::zero(&zero.space, &comodule.space, G::zero());

        Self {
            domain: zero,
            codomain: comodule,
            map: zero_map,
        }
    }

    fn compose(l: &Self, r: &Self) -> Self {
        debug_assert!(
            l.domain == r.codomain,
            "l-domain and r-codomain should be equal when composing"
        );

        let codomain = l.codomain.clone();
        let domain = r.domain.clone();

        let map = l.map.compose(&r.map);

        Self {
            domain,
            codomain,
            map,
        }
    }

    fn get_codomain(&self) -> Arc<RComodule<G,F>> {
        self.codomain.clone()
    }

    fn get_structure_lines(&self) -> Vec<(usize, usize, usize, String)> {
        todo!()
    }
}