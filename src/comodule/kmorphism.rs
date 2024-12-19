use std::{collections::HashMap, sync::Arc};

use itertools::Itertools;

use crate::{
    comodule::{comodule::Tensor, kcomodule::kBasisElement, ktensor::kTensor},
    linalg::{
        field::Field,
        graded::{BasisIndex, GradedLinearMap, Grading},
        matrix::{FieldMatrix, Matrix},
    },
};

use super::{
    comodule::{Comodule, ComoduleMorphism},
    kcomodule::kComodule,
};

#[derive(Debug, Clone)]
#[allow(non_camel_case_types)]
pub struct kComoduleMorphism<G: Grading, F: Field> {
    pub domain: Arc<kComodule<G, F>>,
    pub codomain: Arc<kComodule<G, F>>,

    pub map: GradedLinearMap<G, F, FieldMatrix<F>>, // Question: Shouldn't this be a module morphism?
}

impl<G: Grading, F: Field> kComoduleMorphism<G, F> {
    fn verify_dimensions(&self) {
        for k in self.domain.space.0.keys() {
            debug_assert!(self.map.maps.contains_key(k));
        }

        for k in self.codomain.space.0.keys() {
            debug_assert!(self.map.maps.contains_key(k));
        }

        for (g, map) in self.map.maps.iter() {
            let dom_dim = self.domain.space.dimension_in_grade(g);
            let codom_dim = self.codomain.space.dimension_in_grade(g);
            assert_eq!(dom_dim, map.domain);
            assert_eq!(codom_dim, map.codomain);
        }
    }
}

impl<G: Grading, F: Field> ComoduleMorphism<G, kComodule<G, F>> for kComoduleMorphism<G, F> {
    fn cokernel(&self) -> Self {
        let cokernel_map = self.map.get_cokernel();

        let coker_space = cokernel_map.codomain_space(kBasisElement::default());

        let coalg = self.codomain.coalgebra.as_ref();
        let tensor = kTensor::generate(&coalg.space, &coker_space);

        tensor.is_correct();

        let pivots = cokernel_map.pivots();

        let m_lut: HashMap<BasisIndex<G>, Vec<(usize, F)>> = self
            .codomain
            .tensor
            .construct
            .iter()
            .map(|((f_gr, f_id), _)| {
                // Transfer a specific codomain grade and id (f_gr, f_id) to a list of elements which it maps to in the cokernel
                let v = (0..cokernel_map
                    .maps
                    .get(f_gr)
                    .map(|map| map.codomain)
                    .unwrap_or(0))
                    .filter_map(|q_index| {
                        let val = cokernel_map.maps.get(f_gr).unwrap().data[q_index][*f_id];
                        match val.is_zero() {
                            true => None,
                            false => Some((q_index, val)),
                        }
                    })
                    .collect();
                ((*f_gr, *f_id), v)
            })
            .collect();

        let coaction: HashMap<G, FieldMatrix<F>> = coker_space
            .0
            .iter()
            .map(|(g, v)| {
                let g_tensor_dimen = tensor.get_dimension(g);
                let mut g_coaction = FieldMatrix::<F>::zero(v.len(), g_tensor_dimen);

                // TODO:
                // THERE SHOULD? BE A FASTER VERSION OF THIS, BUT THIS IS SIMPLER ?
                // THIS COMMENTS WAS WRITTERN BEFORE M_lUT()

                // (domain, codomain)
                for (codom_id, coker_id) in &pivots[g] {
                    let coact_size = self.codomain.tensor.dimensions[g];
                    for codom_coact_id in 0..coact_size {
                        let coact_val =
                            self.codomain.coaction.maps[g].data[codom_coact_id][*codom_id];
                        if !coact_val.is_zero() {
                            let ((alg_gr, alg_id), (mod_gr, mod_id)) =
                                self.codomain.tensor.deconstruct[&(*g, codom_coact_id)];

                            for (target_id, val) in m_lut.get(&(mod_gr, mod_id)).unwrap() {
                                let (_, final_id) =
                                    tensor.construct[&(mod_gr, *target_id)][&(alg_gr, alg_id)];
                                g_coaction.data[final_id][*coker_id] += coact_val * *val;
                            }
                        }
                    }
                }
                (*g, g_coaction)
            })
            .collect();

        let comodule = kComodule {
            coalgebra: self.codomain.coalgebra.clone(),
            space: coker_space,
            coaction: GradedLinearMap::from(coaction),
            tensor,
        };
        comodule.verify();
        let coker_morph = Self {
            codomain: Arc::new(comodule),
            domain: self.codomain.clone(),
            map: cokernel_map,
        };
        coker_morph.verify_dimensions();
        coker_morph
    }

    fn inject_codomain_to_cofree(&self, limit: G, fixed_limit: G) -> Self {
        let mut growing_map: GradedLinearMap<G, F, FieldMatrix<F>> =
            GradedLinearMap::zero_codomain(&self.codomain.space);
        let mut growing_comodule = kComodule::zero_comodule(self.codomain.coalgebra.clone());
        let mut iteration = 0;

        let grades: Vec<G> = growing_map.maps.iter().map(|(g, _)| *g).sorted().filter(|&g| g <= limit).collect();
        let mut prev_grade = 0;

        loop {
            // Get lowest graded pivot element
            let mut pivot = None;
            for grade_id in prev_grade..grades.len() {
                let grade = grades[grade_id];
                let kernel = growing_map
                    .maps
                    .get(&grade)
                    .expect("This should exist")
                    .kernel();

                match kernel.first_non_zero_entry() {
                    Some(loc) => {
                        prev_grade = grade_id;
                        pivot = Some((loc, grade));
                        break;
                    }
                    None => {}
                }
            }
            let (pivot, pivot_grade) = match pivot {
                Some(p) => p,
                None => {
                    break;
                }
            };

            // Create a map to a cofree comodule
            let mut map_to_cofree = GradedLinearMap::empty();

            let alg_to_tens = self
                .codomain
                .tensor
                .construct
                .get(&(pivot_grade, pivot.1))
                .expect("The tensor should exist on the codomain in this grade");

            let coalg_space = &self.codomain.coalgebra.space;
            // CHECK
            for alg_gr in coalg_space.0.keys() {
                let t_gr = *alg_gr + pivot_grade;

                if t_gr > fixed_limit {
                    continue;
                }

                let codomain_len = self.codomain.space.dimension_in_grade(&t_gr);
                let coalg_len = coalg_space.dimension_in_grade(alg_gr);

                if !alg_to_tens.contains_key(&(*alg_gr, 0)) {
                    let zero_map = FieldMatrix::zero(codomain_len, coalg_len);
                    map_to_cofree.maps.insert(t_gr, zero_map);
                    continue;
                };

                let mut map = FieldMatrix::zero(codomain_len, coalg_len);

                for a_id in 0..coalg_len {
                    let (t_gr, t_id) = alg_to_tens.get(&(*alg_gr,a_id)).expect("This BasisIndex should exist on the tensor object in the to inject comodule");
                    let slice = self
                        .codomain
                        .coaction
                        .maps
                        .get(t_gr)
                        .expect("This grade should exist on the coaction of the injecting comodule")
                        .data[*t_id]
                        .as_slice();
                    map.data[a_id].clone_from_slice(slice);
                }

                map_to_cofree.maps.insert(t_gr, map);
            }

            let mut f = kComodule::cofree_comodule(
                self.codomain.coalgebra.clone(),
                iteration,
                pivot_grade,
                fixed_limit,
            );
            growing_comodule.direct_sum(&mut f);
            growing_comodule.verify();

            growing_comodule.tensor.is_correct();

            growing_map.vstack(&mut map_to_cofree);
            iteration += 1;
        }
        let final_morph = Self {
            domain: self.codomain.clone(),
            codomain: Arc::new(growing_comodule),
            map: growing_map,
        };
        final_morph.verify_dimensions();
        final_morph
    }

    fn zero_morphism(comodule: Arc<kComodule<G, F>>) -> Self {
        let codomain = comodule.clone();
        let zero = Arc::new(kComodule::zero_comodule(comodule.coalgebra.clone()));

        // Verify how we want to handle this zero map
        let mut zero_map = GradedLinearMap::empty();

        for (gr, elements) in codomain.space.0.iter() {
            zero_map
                .maps
                .insert(*gr, FieldMatrix::zero(0, elements.len()));
        }

        Self {
            domain: zero,
            codomain: codomain,
            map: zero_map,
        }
    }

    fn compose(l: Self, r: Self) -> Self {
        debug_assert!(
            l.domain == r.codomain,
            "l-domain and r-codomain should be equal when composing"
        );

        let codomain = l.codomain;
        let domain = r.domain;

        let map = l.map.compose(r.map);

        let final_morph = Self {
            domain,
            codomain,
            map,
        };
        final_morph.verify_dimensions();
        final_morph
    }

    fn get_structure_lines(&self) -> Vec<(usize, usize, usize, String)> {
        let mut lines = vec![];

        for (gr, gr_map) in self.map.maps.iter() {
            for el_id in 0..self.domain.space.dimension_in_grade(gr) {
                let els = self.domain.space.0.get(gr).unwrap();
                match els[el_id].primitive {
                    Some(prim_id) => {
                        for t_id in 0..gr_map.codomain {
                            let t_el = &self.codomain.space.0.get(gr).expect("As codomain of the map is non-zero this vector space should contain an element in this grade.")[t_id];
                            if t_el.generator {
                                if !gr_map.data[t_id][el_id].is_zero() {
                                    lines.push((
                                        els[el_id].generated_index,
                                        t_el.generated_index,
                                        gr_map.data[t_id][el_id].as_usize(),
                                        prim_id.to_string(),
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

    fn get_codomain(&self) -> Arc<kComodule<G, F>> {
        self.codomain.clone()
    }
}
