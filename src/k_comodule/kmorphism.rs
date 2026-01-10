use std::sync::atomic::AtomicPtr;

use ahash::HashMap;
use algebra::{abelian::Abelian, field::Field};
use deepsize::DeepSizeOf;
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{
    grading::{grading::Grading, tensor::TensorMap},
    k_comodule::{
        graded_space::GradedLinearMap,
        kcoalgebra::kCoalgebra,
        kcomodule::{kCofreeComodule, kComodule},
    },
    traits::{Coalgebra, CofreeComodule, ComoduleMorphism},
    types::{CoalgebraIndexType, CofreeBasis, ComoduleIndex, ComoduleIndexType},
};

#[derive(Debug, Clone, DeepSizeOf)]
#[allow(non_camel_case_types)]
pub struct kComoduleMorphism<G: Grading, F: Field, M: Abelian<F>> {
    pub map: GradedLinearMap<G, F, M>, // Question: Shouldn't this be a module morphism?
}

impl<G: Grading, F: Field, M: Abelian<F>> kComoduleMorphism<G, F, M> {
    fn verify_dimensions(
        &self,
        domain: &kComodule<G, kCoalgebra<G, F, M>>,
        codomain: &kCofreeComodule<G, kCoalgebra<G, F, M>>,
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

impl<G: Grading, F: Field, M: Abelian<F>> ComoduleMorphism<G, kCoalgebra<G, F, M>>
    for kComoduleMorphism<G, F, M>
{
    fn cokernel(
        &self,
        coalgebra: &kCoalgebra<G, F, M>,
        codomain: &kCofreeComodule<G, kCoalgebra<G, F, M>>,
    ) -> (Self, kComodule<G, kCoalgebra<G, F, M>>) {
        let cokernel_map = self.map.get_cokernel();

        let coker_space = cokernel_map.codomain_space(<M as Abelian<F>>::Generator::default());

        let tensor = TensorMap::generate(&coalgebra.space, &coker_space);

        let pivots = cokernel_map.pivots();

        // The upcoming should be a "solve commutative square thingy ?"

        let m_lut: HashMap<ComoduleIndex<G>, Vec<(ComoduleIndexType, F)>> = codomain
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

        let mut find: HashMap<CofreeBasis<G>, ComoduleIndex<G>> = HashMap::default();
        for (gr, v) in &codomain.space.0 {
            for (id, (el, _)) in v.iter().enumerate() {
                find.insert(((el.0.0, el.0.1 as u16), el.1), (*gr, id as u32));
            }
        }

        let coaction: HashMap<G, M> = coker_space
            .0
            .par_iter()
            .map(|(g, v)| {
                let g_tensor_dimen = tensor.get_dimension(g);
                let mut g_coaction_transposed = M::zero(g_tensor_dimen, v.len());

                let coact_ref = &mut g_coaction_transposed as *mut M;
                let map_ref = AtomicPtr::new(coact_ref);

                let space = &codomain.space.0[g];

                pivots[g].par_iter().for_each(|(codom_id, coker_id)| {
                    let (((alg_gr, alg_id), gen_id), _) = space[*codom_id];

                    for ((alg_l_gr, alg_l_id), (alg_r_gr, alg_r_id), coact_val) in
                        coalgebra.coaction.get(&(alg_gr, alg_id)).unwrap()
                    {
                        let (mod_gr, mod_id) = find[&((*alg_r_gr, *alg_r_id), gen_id)];

                        for (target_id, val) in m_lut.get(&(mod_gr, mod_id)).unwrap() {
                            let (_, final_id) =
                                tensor.construct[&(mod_gr, *target_id)][&(*alg_l_gr, *alg_l_id)];

                            // TODO
                            // As coker_id is seperate across parallel instances
                            // SPECIAL care should be taken here
                            unsafe {
                                (**(map_ref.as_ptr())).add_at(
                                    final_id as usize,
                                    *coker_id,
                                    *coact_val * *val,
                                );
                            }
                        }
                    }
                });
                (*g, g_coaction_transposed.transpose())
            })
            .collect();

        let comodule = kComodule::new(coker_space, GradedLinearMap::from(coaction), tensor);

        (kComoduleMorphism { map: cokernel_map }, comodule)
    }

    fn inject_codomain_to_cofree(
        coalgebra: &kCoalgebra<G, F, M>,
        comodule: &kComodule<G, kCoalgebra<G, F, M>>,
        limit: G,
    ) -> (Self, kCofreeComodule<G, kCoalgebra<G, F, M>>) {
        let mut growing_map: GradedLinearMap<G, F, M> =
            GradedLinearMap::zero_codomain(&comodule.space);
        let mut growing_comodule = kCofreeComodule::zero_comodule();
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
                    .kernel_destroyers(&vec![], &vec![]);

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
                    let (t_gr, t_id) = alg_to_tens.get(&(*alg_gr, a_id as CoalgebraIndexType)).expect("This BasisIndex should exist on the tensor object in the to inject comodule");
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

            let mut f = coalgebra.cofree_comodule(
                iteration,
                pivot_grade,
                fixed_limit,
                <M as Abelian<F>>::Generator::default(),
            );

            if cfg!(debug_assertions) {}

            kCofreeComodule::direct_sum(&mut growing_comodule, &mut f);
            growing_map.vstack(&mut GradedLinearMap::from(cofree_map));

            iteration += 1;
        }

        let m = kComoduleMorphism { map: growing_map };
        assert!(m.verify_dimensions(&comodule, &growing_comodule));
        (m, growing_comodule)
    }

    fn zero_morphism(comodule: &kComodule<G, kCoalgebra<G, F, M>>) -> Self {
        // Verify how we want to handle this zero map
        let mut zero_map = GradedLinearMap::empty();

        for (gr, elements) in comodule.space.0.iter() {
            zero_map.maps.insert(*gr, M::zero(0, elements.len()));
        }

        Self { map: zero_map }
    }

    // domain l == codomain r, l \circ r
    fn compose(l: &Self, r: &Self, _codomain: &kCofreeComodule<G, kCoalgebra<G, F, M>>) -> Self {
        let map = GradedLinearMap::<G, F, M>::compose(&l.map, &r.map);

        Self { map }
    }
}
