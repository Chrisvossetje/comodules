use std::{sync::Arc};

use ahash::HashMap;
use algebra::{abelian::Abelian, field::Field, matrix::Matrix};
use itertools::Itertools;
use rayon::prelude::*;

use crate::{
    basiselement::kBasisElement, graded_space::{BasisIndex, GradedLinearMap}, grading::{Grading, OrderedGrading}, tensor::Tensor
};

use super::{
    kcomodule::kComodule,
    traits::{Comodule, ComoduleMorphism},
};

#[derive(Debug, Clone)]
#[allow(non_camel_case_types)]
pub struct kComoduleMorphism<G: Grading, F: Field, M: Matrix<F>> {
    pub domain: Arc<kComodule<G, F, M>>,
    pub codomain: Arc<kComodule<G, F, M>>,
    pub map: GradedLinearMap<G, F, M>, // Question: Shouldn't this be a module morphism?
}

impl<G: Grading, F: Field, M: Matrix<F>> kComoduleMorphism<G, F, M> {
    fn verify_dimensions(&self) -> bool {
        for k in self.domain.space.0.keys() {
            if !self.map.maps.contains_key(k) {
                return false;
            };
        }

        for k in self.codomain.space.0.keys() {
            if !self.map.maps.contains_key(k) {
                return false;
            };
        }

        for (g, map) in self.map.maps.iter() {
            let dom_dim = self.domain.space.dimension_in_grade(g);
            let codom_dim = self.codomain.space.dimension_in_grade(g);
            if dom_dim != map.domain() {
                return false;
            };
            if codom_dim != map.codomain() {
                return false;
            };
        }

        true
    }

    pub fn new(
        domain: Arc<kComodule<G, F, M>>,
        codomain: Arc<kComodule<G, F, M>>,
        map: GradedLinearMap<G, F, M>,
    ) -> Self {
        let new = Self {
            domain,
            codomain,
            map,
        };
        debug_assert!(new.verify_dimensions());
        new
    }
}

impl<G: Grading, F: Field, M: Abelian<F>> ComoduleMorphism<G, kComodule<G, F, M>>
    for kComoduleMorphism<G, F, M>
{
    fn cokernel(&self) -> Self {
        let cokernel_map = self.map.get_cokernel();

        let coker_space = cokernel_map.codomain_space(kBasisElement::default());

        let coalg = self.codomain.coalgebra.as_ref();
        let tensor = Tensor::generate(&coalg.space, &coker_space);

        let pivots = cokernel_map.pivots();

        // The upcoming should be a "solve commutative square thingy ?"

        let m_lut: HashMap<BasisIndex<G>, Vec<(usize, F)>> = self
            .codomain
            .tensor
            .construct
            .par_iter()
            .map(|((f_gr, f_id), _)| {
                // Transfer a specific codomain grade and id (f_gr, f_id) to a list of elements which it maps to in the cokernel
                let v = (0..cokernel_map
                    .maps
                    .get(f_gr)
                    .map(|map| map.codomain())
                    .unwrap_or(0))
                    .filter_map(|q_index| {
                        let val = cokernel_map.maps.get(f_gr).unwrap().get(*f_id, q_index);
                        match val.is_zero() {
                            true => None,
                            false => Some((q_index, val)),
                        }
                    })
                    .collect();
                ((*f_gr, *f_id), v)
            })
            .collect();

        let coaction: HashMap<G, M> = coker_space
            .0
            .par_iter()
            .map(|(g, v)| {
                let g_tensor_dimen = tensor.get_dimension(g);
                let mut g_coaction = M::zero(v.len(), g_tensor_dimen);

                // TODO:
                // THERE SHOULD? BE A FASTER VERSION OF THIS, BUT THIS IS SIMPLER ?
                // THIS COMMENTS WAS WRITTERN BEFORE M_lUT()

                // (domain, codomain)
                for (codom_id, coker_id) in &pivots[g] { // TODO : THiS COULD BE PARALLEL
                    let coact_size = self.codomain.tensor.dimensions[g];
                    for codom_coact_id in 0..coact_size {
                        let coact_val =
                            self.codomain.coaction.maps[g].get(*codom_id, codom_coact_id);
                        if !coact_val.is_zero() {
                            let ((alg_gr, alg_id), (mod_gr, mod_id)) =
                                self.codomain.tensor.deconstruct[&(*g, codom_coact_id)];

                            for (target_id, val) in m_lut.get(&(mod_gr, mod_id)).unwrap() {
                                let (_, final_id) =
                                    tensor.construct[&(mod_gr, *target_id)][&(alg_gr, alg_id)];
                                g_coaction.add_at(*coker_id, final_id, coact_val * *val);
                            }
                        }
                    }
                }
                (*g, g_coaction)
            })
            .collect();

        let comodule = kComodule::new(
            self.codomain.coalgebra.clone(),
            coker_space,
            GradedLinearMap::from(coaction),
            tensor,
        );

        Self::new(self.codomain.clone(), Arc::new(comodule), cokernel_map)
    }

    fn inject_codomain_to_cofree(&self, limit: G) -> Self 
    where
        G: OrderedGrading,
    {
        let mut growing_map: GradedLinearMap<G, F, M> =
            GradedLinearMap::zero_codomain(&self.codomain.space);
        let mut growing_comodule = kComodule::zero_comodule(self.codomain.coalgebra.clone());
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

            let alg_to_tens = self
                .codomain
                .tensor
                .construct
                .get(&(pivot_grade, pivot))
                .expect("The tensor should exist on the codomain in this grade");

            let coalg_space = &self.codomain.coalgebra.space;

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

                let codomain_len = self.codomain.space.dimension_in_grade(&t_gr);
                let coalg_len = alg_gr_space.len();

                if !alg_to_tens.contains_key(&(*alg_gr, 0)) {
                    let zero_map = M::zero(codomain_len, coalg_len);
                    return Some((t_gr, zero_map));
                };

                let mut map = M::zero(codomain_len, coalg_len);

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

            let mut f = kComodule::cofree_comodule(
                self.codomain.coalgebra.clone(),
                iteration,
                pivot_grade,
                fixed_limit,
                ()
            );

            growing_comodule.direct_sum(&mut f);
            growing_map.vstack(&mut GradedLinearMap::from(cofree_map));

            iteration += 1;
        }

        Self::new(
            self.codomain.clone(),
            Arc::new(growing_comodule),
            growing_map,
        )
    }

    fn zero_morphism(comodule: Arc<kComodule<G, F, M>>) -> Self {
        let codomain = comodule.clone();
        let zero = Arc::new(kComodule::zero_comodule(comodule.coalgebra.clone()));

        // Verify how we want to handle this zero map
        let mut zero_map = GradedLinearMap::empty();

        for (gr, elements) in codomain.space.0.iter() {
            zero_map.maps.insert(*gr, M::zero(0, elements.len()));
        }

        Self {
            domain: zero,
            codomain: codomain,
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

        let map = GradedLinearMap::<G, F, M>::compose(&l.map, &r.map);

        Self::new(domain, codomain, map)
    }

    fn get_structure_lines(&self) -> Vec<((usize, G, usize), (usize, G, usize), F, String)> {
        let mut lines = vec![];

        for (gr, gr_map) in self.map.maps.iter() {
            for el_id in 0..self.domain.space.dimension_in_grade(gr) {
                let els = self.domain.space.0.get(gr).unwrap();
                match els[el_id].primitive {
                    Some(prim_id) => {
                        for t_id in 0..gr_map.codomain() {
                            let t_el = &self.codomain.space.0.get(gr).expect("As codomain of the map is non-zero this vector space should contain an element in this grade.")[t_id];
                            if t_el.generator {
                                if !gr_map.get(el_id, t_id).is_zero() {
                                    lines.push((
                                        // TODO 
                                        (els[el_id].generated_index, G::zero(), 0),
                                        (t_el.generated_index, G::zero(), 0),
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

    fn get_codomain(&self) -> Arc<kComodule<G, F, M>> {
        self.codomain.clone()
    }
}
