use std::marker::PhantomData;

use ahash::HashMap;
use algebra::{
    abelian::Abelian, field::Field, ring::CRing, rings::univariate_polynomial_ring::UniPolRing,
};
use deepsize::DeepSizeOf;
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{
    grading::{grading::Grading, tensor::TensorMap},
    k_t_comodule::{
        graded_module::GradedktFieldMap,
        k_t_coalgebra::ktCoalgebra,
        k_t_comodule::{ktCofreeComodule, ktComodule},
    },
    traits::{Coalgebra, CofreeComodule, Comodule, ComoduleMorphism},
    types::{CoalgebraIndexType, CofreeBasis, ComoduleIndex, ComoduleIndexType},
};

#[derive(Debug, Clone, DeepSizeOf)]
#[allow(non_camel_case_types)]
pub struct ktComoduleMorphism<G: Grading, F: Field, M: Abelian<UniPolRing<F>>> {
    pub map: GradedktFieldMap<G, F, M>,
}

impl<G: Grading, F: Field, M: Abelian<UniPolRing<F>>> ComoduleMorphism<G, ktCoalgebra<G, F, M>>
    for ktComoduleMorphism<G, F, M>
{
    fn cokernel(
        &self,
        coalgebra: &ktCoalgebra<G, F, M>,
        codomain: &ktCofreeComodule<G, ktCoalgebra<G, F, M>>,
    ) -> (Self, ktComodule<G, F, M>) {
        // TODO :
        // if cfg!(debug_assertions) {
        //     self.map.verify(&self.domain.space, &self.codomain.space).unwrap();
        //     self.codomain.verify().unwrap()
        // }

        let (coker_to, coker_inv, coker) = self.map.cokernel(&codomain.to_module_generators());

        let tensor = TensorMap::generate(&coalgebra.space, &coker);

        let m_lut: HashMap<ComoduleIndex<G>, Vec<(ComoduleIndexType, UniPolRing<F>)>> = codomain
            .space
            .0
            // the choice for the tensor here is not neccessary
            // Could also be self.codomain.space and iterate over len of the module
            .par_iter()
            .flat_map(|(f_gr, m)| {
                // Transfer a specific codomain grade and id (f_gr, f_id) to a list of elements which it maps to in the cokernel
                (0..m.len()).into_par_iter().map(|f_id| {
                    let v = (0..coker_to
                        .maps
                        .get(f_gr)
                        .map(|map| map.codomain())
                        .unwrap_or(0))
                        .filter_map(|q_index| {
                            let val = coker_to.maps.get(f_gr).unwrap().get(f_id, q_index);
                            match val.is_zero() {
                                true => None,
                                false => Some((q_index as u32, val)),
                            }
                        })
                        .collect();
                    ((*f_gr, f_id as u32), v)
                })
            })
            .collect();

        let mut find: HashMap<CofreeBasis<G>, ComoduleIndex<G>> = HashMap::default();
        for (gr, v) in &codomain.space.0 {
            for (id, (el, _)) in v.iter().enumerate() {
                find.insert(((el.0.0, el.0.1 as u16), el.1), (*gr, id as u32));
            }
        }

        // TODO :
        // let coker_tensor_module = coker.generate_tensor_as_module(&self.codomain.coalgebra.space, &tensor);

        let coaction = coker
            .0
            .par_iter()
            .map(|(g, v)| {
                let g_tensor_dimen = tensor.get_dimension(g);
                let mut g_coaction = M::zero(v.len(), g_tensor_dimen);

                let inv_map = coker_inv.maps.get(&g).unwrap();
                let space = &codomain.space.0[g];

                (0..v.len()).for_each(|coker_id| {
                    for codom_id in 0..inv_map.codomain() {
                        let inv_val = inv_map.get(coker_id, codom_id);
                        if inv_val.is_zero() {
                            continue;
                        }

                        let (((alg_gr, alg_id), gen_id), _) = space[codom_id];

                        for ((alg_l_gr, alg_l_id), (alg_r_gr, alg_r_id), coact_val) in
                            coalgebra.coaction.get(&(alg_gr, alg_id)).unwrap()
                        {
                            let (mod_gr, mod_id) = find[&((*alg_r_gr, *alg_r_id), gen_id)];

                            for (target_id, val) in m_lut.get(&(mod_gr, mod_id)).unwrap() {
                                let (_, final_id) = tensor.construct[&(mod_gr, *target_id)]
                                    [&(*alg_l_gr, *alg_l_id)];

                                let final_val = inv_val * *coact_val * *val;

                                // TODO : Coker tensor module thing!
                                // if let Some(tens_el_power) = tensor_el.2 {
                                //     if final_val.1 >= tens_el_power {
                                //         continue;
                                //     }
                                // }
                                g_coaction.add_at(coker_id, final_id as usize, final_val);
                            }
                        }
                    }
                });
                (*g, g_coaction)
            })
            .collect();

        let coaction = GradedktFieldMap {
            maps: coaction,
            _p: PhantomData,
        };

        let coker_comod = ktComodule {
            space: coker,
            coaction,
            tensor,
        };

        (Self { map: coker_to }, coker_comod)
    }

    fn inject_codomain_to_cofree(
        coalgebra: &ktCoalgebra<G, F, M>,
        comodule: &ktComodule<G, F, M>,
        limit: G,
    ) -> (Self, ktCofreeComodule<G, ktCoalgebra<G, F, M>>) {
        let mut growing_map: GradedktFieldMap<G, F, M> =
            GradedktFieldMap::zero_codomain(&comodule.space);
        let mut growing_comodule = ktCofreeComodule::zero_comodule();
        let mut iteration = 0;

        let grades: Vec<G> = growing_map
            .maps
            .iter()
            .map(|(g, _)| *g)
            .sorted_by(|a, b| a.compare(b))
            .filter(|&g| g.compare(&limit).is_le())
            .collect();

        let fixed_limit = limit.incr().incr();

        for pivot_grade in grades {
            // Get lowest graded pivot element
            let map = growing_map.maps.get(&pivot_grade).unwrap();
            let empty = vec![];
            let domain: &Vec<_> = &comodule.space.0.get(&pivot_grade).unwrap_or(&empty);
            let codomain: Vec<_> = growing_comodule
                .space
                .0
                .get(&pivot_grade)
                .map_or(vec![], |x| {
                    x.iter()
                        .map(
                            |y: &(((G, u16), u16), <M as Abelian<UniPolRing<F>>>::Generator)| {
                                y.1.clone()
                            },
                        )
                        .collect()
                });

            let l = map.kernel_destroyers(domain, &codomain);
            for pivot in l {
                let generator = comodule
                    .space
                    .0
                    .get(&pivot_grade)
                    .unwrap()
                    .get(pivot)
                    .unwrap();

                let alg_to_tens = comodule
                    .tensor
                    .construct
                    .get(&(pivot_grade, pivot as u32))
                    .expect("The tensor should exist on the codomain in this grade");

                let coalg_space = &coalgebra.space;

                let cofree_map: HashMap<G,M> = coalg_space.0.iter().filter_map(|(alg_gr, alg_gr_space)| {
                    let t_gr = *alg_gr + pivot_grade;

                    if t_gr > fixed_limit {
                        return None;
                    }

                    let codomain_len = comodule.space.dimension_in_grade(&t_gr);
                    let coalg_len = alg_gr_space.len();

                    if !alg_to_tens.contains_key(&(*alg_gr, 0)) {
                        let zero_map = M::zero(codomain_len, coalg_len);
                        return Some((t_gr, zero_map));
                    };

                    let mut map = M::zero(codomain_len, coalg_len);

                    for a_id in 0..coalg_len {
                        let (t_gr, t_id) = alg_to_tens.get(&(*alg_gr,a_id as CoalgebraIndexType)).expect("This BasisIndex should exist on the tensor object in the to inject comodule");
                        let slice = comodule
                            .coaction
                            .maps
                            .get(t_gr)
                            .expect("This grade should exist on the coaction of the injecting comodule")
                            .get_row(*t_id as usize);

                        map.set_row(a_id, slice);
                    }

                    Some((t_gr, map))
                }).collect();

                let mut f =
                    coalgebra.cofree_comodule(iteration, pivot_grade, fixed_limit, *generator);

                let mut extend_map = GradedktFieldMap {
                    maps: cofree_map,
                    _p: PhantomData,
                };

                growing_comodule.direct_sum(&mut f);
                growing_map.vstack(&mut extend_map);

                // TODO :
                // if cfg!(debug_assertions) {

                //     growing_map.verify(&comodule.space, &growing_comodule.space).unwrap();
                // }

                iteration += 1;
            }
        }

        (Self { map: growing_map }, growing_comodule)
    }

    fn zero_morphism(comodule: &ktComodule<G, F, M>) -> Self {
        let zero: ktComodule<G, F, M> = ktComodule::zero_comodule();

        let zero_map = GradedktFieldMap::zero(&zero.space, &comodule.space);

        Self { map: zero_map }
    }

    fn compose(l: &Self, r: &Self, _codomain: &ktCofreeComodule<G, ktCoalgebra<G, F, M>>) -> Self {
        let map = l.map.compose(&r.map);

        Self { map }
    }
    // fn cokernel(&self) -> Self {
    //
    // }

    // fn inject_codomain_to_cofree(&self, limit: G) -> Self
    // where
    //     G: OrderedGrading {
    //
    // }

    // /// Zero should be the domain, the comodule is the codomain
    // fn zero_morphism(comodule: Arc<ktComodule<G,F>>) -> Self {

    // }

    // TODO : This could be removed but might be useful ?

    // // (generated index , Grade , real index )
    // fn get_structure_lines(&self) -> Vec<((usize, G, usize), (usize, G, usize), UniPolRing<F>, String)> {
    //     let mut lines = vec![];

    //     for (gr, gr_map) in self.map.maps.iter() {
    //         for el_id in 0..self.domain.space.dimension_in_grade(gr) {
    //             let els = self.domain.space.0.get(gr).unwrap();
    //             match els[el_id].0.primitive {
    //                 Some(prim_id) => {
    //                     for t_id in 0..gr_map.codomain() {
    //                         let t_el = &self.codomain.space.0.get(gr).expect("As codomain of the map is non-zero this vector space should contain an element in this grade.")[t_id];
    //                         if t_el.0.generator {
    //                             if !gr_map.get(el_id, t_id).is_zero() {
    //                                 let s_gen_id = els[el_id].0.generated_index;
    //                                 // find source generator
    //                                 let (s_gr, s_id) = {
    //                                     let mut a = None;
    //                                     'outer: for (gr, module) in &self.domain.space.0 {
    //                                         for (id,el) in module.iter().enumerate() {
    //                                             if el.0.generator && el.0.generated_index == s_gen_id {
    //                                                 a = Some((*gr, id));
    //                                                 break 'outer;
    //                                             }
    //                                         }
    //                                     }
    //                                     match a {
    //                                         Some(b) => b,
    //                                         None => {panic!()},
    //                                     }
    //                                 };

    //                                 lines.push((
    //                                     (s_gen_id, s_gr, s_id),
    //                                     (t_el.0.generated_index, *gr, t_id),
    //                                     gr_map.get(el_id, t_id),
    //                                     "h_".to_string() + &prim_id.to_string(),
    //                                 ));
    //                             }
    //                         }
    //                     }
    //                 }
    //                 None => {}
    //             }
    //         }
    //     }

    //     lines
    // }
}
