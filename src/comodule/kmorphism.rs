use std::sync::{Arc, atomic::AtomicPtr};

use ahash::{HashMap, HashMapExt};
use algebra::{abelian::Abelian, field::Field, matrix::Matrix};
use deepsize::DeepSizeOf;
use itertools::Itertools;
use rayon::prelude::*;

use crate::{
    basiselement::kBasisElement,
    comodule::{kcomodule::{Original, kCofreeComodule}, traits::CofreeComodule},
    graded_space::{BasisIndex, GradedLinearMap},
    grading::{Grading, OrderedGrading},
    tensor::TensorMap,
};

use super::{
    kcomodule::kComodule,
    traits::{Comodule, ComoduleMorphism},
};

#[derive(Debug, Clone, DeepSizeOf)]
#[allow(non_camel_case_types)]
pub struct kComoduleMorphism<G: Grading, F: Field, M: Matrix<F>> {
    pub map: GradedLinearMap<G, F, M>, // Question: Shouldn't this be a module morphism?
}

impl<G: Grading, F: Field, M: Matrix<F>> kComoduleMorphism<G, F, M> {
    fn verify_dimensions(
        &self,
        domain: &kComodule<G, F, M>,
        codomain: &kCofreeComodule<G, F, M>,
    ) -> bool {
        for k in domain.space.0.keys() {
            if !self.map.maps.contains_key(k) {
                return false;
            };
        }

        for k in codomain.space.0.keys() {
            if !self.map.maps.contains_key(k) {
                return false;
            };
        }

        for k in self.map.maps.keys() {
            if !codomain.space.0.contains_key(&k) {
                return false;
            }
        }

        for (g, map) in self.map.maps.iter() {
            let dom_dim = domain.space.dimension_in_grade(g);
            let codom_dim = codomain.space.dimension_in_grade(g);
            if dom_dim != map.domain() {
                return false;
            };
            if codom_dim != map.codomain() {
                return false;
            };
        }

        true
    }

    pub fn new(map: GradedLinearMap<G, F, M>) -> Self {
        let new = Self { map };
        new
    }
}

impl<G: Grading, F: Field, M: Abelian<F>> ComoduleMorphism<G> for kComoduleMorphism<G, F, M> {
    type CofreeComodule = kCofreeComodule<G, F, M>;
    type Comodule = kComodule<G, F, M>;
    type BaseRing = F;

    fn cokernel(&self, codomain: &Self::CofreeComodule) -> (Self, Self::Comodule) {
        let cokernel_map = self.map.get_cokernel();

        let coker_space = cokernel_map.codomain_space(());

        let coalg = codomain.coalgebra.as_ref();
        let tensor = TensorMap::generate(&coalg.space, &coker_space);

        let pivots = cokernel_map.pivots();

        // The upcoming should be a "solve commutative square thingy ?"

        let m_lut: HashMap<BasisIndex<G>, Vec<(u32, F)>> = codomain
            .space
            .0
            .par_iter()
            .flat_map(|(f_gr, m)| {
                // Transfer a specific codomain grade and id (f_gr, f_id) to a list of elements which it maps to in the cokernel
                (0..m.len()).into_par_iter().map(|f_id| {
                    let v = (0..cokernel_map
                        .maps
                        .get(f_gr)
                        .map(|map| map.codomain())
                        .unwrap_or(0))
                        .filter_map(|q_index| {
                            let val = cokernel_map.maps.get(f_gr).unwrap().get(f_id, q_index);
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

        let mut find: HashMap<Original<G>, BasisIndex<G>> = HashMap::new();
        for (gr, v) in &codomain.space.0 {
            for (id, (el, _)) in v.iter().enumerate() {
                find.insert(*el, (*gr, id as u32));
            }
        }

        let coaction: HashMap<G, M> = coker_space
            .0
            .par_iter()
            .map(|(g, v)| {
                let g_tensor_dimen = tensor.get_dimension(g);
                let mut g_coaction = M::zero(v.len(), g_tensor_dimen);

                let coact_ref = &mut g_coaction as *mut M;
                let a = AtomicPtr::new(coact_ref);

                let space = &codomain.space.0[g];

                pivots[g].par_iter().for_each(|(codom_id, coker_id)| {
                    let (((alg_gr, alg_id), gen_id), _) = space[*codom_id];

                    let coact_size = codomain.coalgebra.tensor.dimensions[&alg_gr];
                    for codom_coact_id in 0..coact_size as u32 {
                        let coact_val =
                            codomain.coalgebra.coaction.maps[&alg_gr].get(alg_id as usize, codom_coact_id as usize);

                        if !coact_val.is_zero() {
                            let ((alg_l_gr, alg_l_id), (alg_r_gr, alg_r_id)) =
                                codomain.coalgebra.tensor.deconstruct[&(alg_gr, codom_coact_id)];

                            let (mod_gr, mod_id) = find[&((alg_r_gr, alg_r_id), gen_id)];

                            for (target_id, val) in m_lut.get(&(mod_gr, mod_id)).unwrap() {
                                let (_, final_id) =
                                    tensor.construct[&(mod_gr, *target_id)][&(alg_l_gr, alg_l_id)];

                                // As coker_id is seperate across parallel instances
                                // This unsafe code is fine, AS LONG as the matrix is FlatMatrix :)
                                // TODO : This is not reallly generic, and depends on the underlying implementation
                                // This probably breaks for a F2 matrix implementation.
                                unsafe {
                                    (**(a.as_ptr())).add_at(*coker_id, final_id as usize, coact_val * *val);
                                }
                            }
                        }
                    }
                });
                (*g, g_coaction)
            })
            .collect();

        let comodule = kComodule::new(
            codomain.coalgebra.clone(),
            coker_space,
            GradedLinearMap::from(coaction),
            tensor,
        );

        (kComoduleMorphism { map: cokernel_map }, comodule)
    }

    fn inject_codomain_to_cofree(
        comodule: &Self::Comodule,
        limit: G,
    ) -> (Self, Self::CofreeComodule)
    where
        G: OrderedGrading,
    {
        let mut growing_map: GradedLinearMap<G, F, M> =
            GradedLinearMap::zero_codomain(&comodule.space);
        let mut growing_comodule = kCofreeComodule::zero_comodule(comodule.coalgebra.clone());
        let mut iteration = 0;

        let grades: Vec<G> = growing_map
            .maps
            .iter()
            .map(|(g, _)| *g)
            .sorted_by(|a, b| a.compare(b))
            .filter(|&g| g.compare(&limit).is_le())
            .collect();
        let mut prev_grade = 0;

        let fixed_limit = limit.incr().incr();

        loop {
            // Get lowest graded pivot element
            let mut pivot = None;
            for grade_id in prev_grade..grades.len() {
                let grade = grades[grade_id];
                let kernel = growing_map
                    .maps
                    .get(&grade)
                    .expect("This should exist")
                    .kernel_generators(&M::Module::default(), &M::Module::default());

                match kernel.first() {
                    Some(loc) => {
                        prev_grade = grade_id;
                        pivot = Some((*loc, grade));
                        break;
                    }
                    None => {}
                }
            }
            let (pivot, pivot_grade) = match pivot {
                Some(p) => p,
                _ => {
                    break;
                }
            };

            let alg_to_tens = comodule
                .tensor
                .construct
                .get(&(pivot_grade, pivot as u32))
                .expect("The tensor should exist on the codomain in this grade");

            let coalg_space = &comodule.coalgebra.space;

            // TODO: Verify is this parallel iterator is faster or not for big(ger) coalgebras
            let cofree_map: HashMap<G,M> = coalg_space.0.iter().filter_map(|(alg_gr, alg_gr_space)| {
                let t_gr = *alg_gr + pivot_grade;

                if t_gr > fixed_limit {
                    return None;
                }
                // TODO ! UNSTABLE STUFF, THIS IS PROB WRONG
                // if alg_gr > &pivot_grade {
                //     return None;
                // }

                let codomain_len = comodule.space.dimension_in_grade(&t_gr);
                let coalg_len = alg_gr_space.len();

                if !alg_to_tens.contains_key(&(*alg_gr, 0)) {
                    let zero_map = M::zero(codomain_len, coalg_len);
                    return Some((t_gr, zero_map));
                };

                let mut map = M::zero(codomain_len, coalg_len);

                for a_id in 0..coalg_len {
                    let (t_gr, t_id) = alg_to_tens.get(&(*alg_gr, a_id as u32)).expect("This BasisIndex should exist on the tensor object in the to inject comodule");
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

            let mut f = kCofreeComodule::cofree_comodule(
                comodule.coalgebra.clone(),
                iteration,
                pivot_grade,
                fixed_limit,
                (),
            );

            if cfg!(debug_assertions) {}

            growing_comodule.direct_sum(&mut f);
            growing_map.vstack(&mut GradedLinearMap::from(cofree_map));

            iteration += 1;
        }

        let m = kComoduleMorphism { map: growing_map };
        assert!(m.verify_dimensions(&comodule, &growing_comodule));
        (m, growing_comodule)
    }

    fn zero_morphism(comodule: &Self::Comodule) -> Self {
        // Verify how we want to handle this zero map
        let mut zero_map = GradedLinearMap::empty();

        for (gr, elements) in comodule.space.0.iter() {
            zero_map.maps.insert(*gr, M::zero(0, elements.len()));
        }

        Self { map: zero_map }
    }

    // domain l == codomain r, l \circ r
    fn compose(l: &Self, r: &Self, _codomain: &Self::CofreeComodule) -> Self {
        let map = GradedLinearMap::<G, F, M>::compose(&l.map, &r.map);

        Self { map }
    }

    /// (s, gen_index) uniquely defines a generator of Ext
    /// in a specific morphism we only need to know its gen_index
    /// in the resolution we add the s
    /// (from_dot, to_dot, value, line_type)
    fn get_structure_lines(
        &self,
        domain: &Self::CofreeComodule,
        codomain: &Self::CofreeComodule,
    ) -> Vec<((usize, G, usize), (usize, G, usize), Self::BaseRing, String)> {
        let mut lines = vec![];

        let c_space = &domain.coalgebra.space.0;

        for (gr, gr_map) in self.map.maps.iter() {
            for el_id in 0..domain.space.dimension_in_grade(gr) {
                let els = domain.space.0.get(gr).unwrap();
                let (((alg_gr, alg_id), gen_id), _) = els[el_id];
                let alg_el = &c_space[&alg_gr][alg_id as usize];
                match alg_el.primitive {
                    Some(prim_id) => {
                        for t_id in 0..gr_map.codomain() {
                            let (((codom_alg_gr, codom_alg_id), codom_gen_id), _) = codomain.space.0.get(gr).expect("As codomain of the map is non-zero this vector space should contain an element in this grade.")[t_id];

                            let alg_el = &c_space[&codom_alg_gr][codom_alg_id as usize];
                            if alg_el.generator {
                                if !gr_map.get(el_id, t_id).is_zero() {
                                    lines.push((
                                        // TODO
                                        (gen_id as usize, G::zero(), 0),
                                        (codom_gen_id as usize, G::zero(), 0),
                                        gr_map.get(el_id, t_id),
                                        "h_".to_string() + &prim_id.to_string(),
                                    ));
                                }
                            }
                        }
                    }
                    None => {}
                }
            }
        }

        lines
    }

    fn direct_sum(
        a: &mut Self,
        b: &mut Self,
        a_codom: &mut Self::CofreeComodule,
        b_codom: &mut Self::CofreeComodule,
    ) {
        todo!()
    }
}
