use std::collections::HashMap;

use ahash::RandomState;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::linalg::{
    graded::{BasisElement, BasisIndex, GradedVectorSpace},
    grading::Grading,
};

use super::traits::Tensor;

pub type TensorConstruct<G> =
    HashMap<BasisIndex<G>, HashMap<BasisIndex<G>, BasisIndex<G>, RandomState>, RandomState>;
pub type TensorDeconstruct<G> = HashMap<BasisIndex<G>, (BasisIndex<G>, BasisIndex<G>), RandomState>;
pub type TensorDimension<G> = HashMap<G, usize, RandomState>;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[allow(non_camel_case_types)]
pub struct kTensor<G: Grading> {
    // Module BasisIndex -> Coalgebra BasisIndex -> Tensor BasisIndex
    pub construct: TensorConstruct<G>,

    // Tensor BasisIndex -> (Coalgebra BasisIndex, Comodule BasisIndex)
    pub deconstruct: TensorDeconstruct<G>,
    
    // Tensor Grade -> usize (Dimension of tensor in that grade)
    pub dimensions: TensorDimension<G>,
}

impl<G: Grading> Tensor<G> for kTensor<G> {
    fn tensor_to_base() {
        unimplemented!()
    }

    fn base_to_tensor() {
        unimplemented!()
    }

    fn get_dimension(&self, grading: &G) -> usize {
        *self.dimensions.get(grading).unwrap_or(&0)
    }
}

impl<G: Grading> Default for kTensor<G> {
    fn default() -> Self {
        Self {
            construct: Default::default(),
            deconstruct: Default::default(),
            dimensions: Default::default(),
        }
    }
}

impl<G: Grading> kTensor<G> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn generate<B: BasisElement>(
        left: &GradedVectorSpace<G, B>,
        right: &GradedVectorSpace<G, B>,
    ) -> Self {
        let mut construct = HashMap::default();
        let mut deconstruct = HashMap::default();
        let mut dimensions = HashMap::default();

        // This sorted is important for consistency !
        for (r_grade, r_elements) in right.0.iter().sorted_by_key(|(&g, _)| g) {
            for r_id in 0..r_elements.len() {
                let entry = construct
                    .entry((*r_grade, r_id))
                    .or_insert(HashMap::default());
                for (l_grade, l_elements) in left.0.iter() {
                    let t_grade = *l_grade + *r_grade;

                    if !right.0.contains_key(&t_grade) {
                        continue;
                    }
                    let t_id = dimensions.entry(t_grade).or_insert(0);

                    for l_id in 0..l_elements.len() {
                        entry.insert((*l_grade, l_id), (t_grade, *t_id));
                        deconstruct.insert((t_grade, *t_id), ((*l_grade, l_id), (*r_grade, r_id)));
                        *t_id = *t_id + 1;
                    }
                }
            }
        }

        let tensor = Self {
            construct,
            deconstruct,
            dimensions,
        };
        debug_assert!(tensor.is_correct());
        tensor
    }

    pub fn add_and_restrict(&self, add: G, limit: G) -> kTensor<G> {
        let cons = self
            .construct
            .iter()
            .filter_map(|((m_gr, m_id), map)| {
                let m_sum = *m_gr + add;
                if m_sum > limit {
                    return None;
                }

                let alg_map = map
                    .iter()
                    .filter_map(|((a_gr, a_id), (t_gr, t_id))| {
                        let t_sum = *t_gr + add;
                        if t_sum <= limit {
                            Some(((*a_gr, *a_id), (t_sum, *t_id)))
                        } else {
                            None
                        }
                    })
                    .collect();

                Some(((m_sum, *m_id), alg_map))
            })
            .collect();

        let decon = self
            .deconstruct
            .iter()
            .filter_map(|((t_gr, t_id), ((a_gr, a_id), (m_gr, m_id)))| {
                let t_sum = *t_gr + add;
                let m_sum = *m_gr + add;
                if t_sum <= limit {
                    Some(((t_sum, *t_id), ((*a_gr, *a_id), (m_sum, *m_id))))
                } else {
                    None
                }
            })
            .collect();

        let dims = self
            .dimensions
            .iter()
            .filter_map(|(t_gr, size)| {
                let t_sum = *t_gr + add;
                if t_sum <= limit {
                    Some((t_sum, *size))
                } else {
                    None
                }
            })
            .collect();
        let tens = kTensor {
            construct: cons,
            deconstruct: decon,
            dimensions: dims,
        };
        debug_assert!(tens.is_correct());
        tens
    }



    /// Special care should be taken here,
    /// as this tensor object usually relates to certain graded linear maps
    /// We should expect that direct summing the underlying vector spaces creates the correct new tensor object
    ///
    /// After careful thinking, this direct sum is not dependent on the non-determinism of the hashmap
    pub fn direct_sum_add_g(&mut self, other: &Self, self_space_dimensions: &HashMap<G, usize>, grade: G, limit: G) {
        // TODO: Consider using parallel here for bigger coalgebras ?
        self.construct
            .extend(other.construct.iter().filter_map(|(&(m_gr, m_id), maps)| {
                let new_m_gr = m_gr + grade;
                if new_m_gr > limit {
                    return None
                }

                let m_id_new = self_space_dimensions.get(&new_m_gr).unwrap_or(&0) + m_id;

                let new_map = maps
                    .iter()
                    .filter_map(|((a_gr, a_id), (t_gr, t_id))| {
                        let new_t_gr = *t_gr + grade;
                        if new_t_gr > limit {
                            return None
                        }

                        let t_id_new = self.dimensions.get(&new_t_gr).unwrap_or(&0) + t_id;
                        Some(((*a_gr, *a_id), (new_t_gr, t_id_new)))
                    })
                    .collect();

                Some(((new_m_gr, m_id_new), new_map))
            }));

        self.deconstruct.extend(other.deconstruct.iter().filter_map(
            |(&(t_gr, t_id), &((a_gr, a_id), (m_gr, m_id)))| {
                let new_m_gr = m_gr + grade;
                let new_t_gr = t_gr + grade;

                if new_t_gr > limit {
                   return None
                }

                let m_id_new = self_space_dimensions.get(&new_m_gr).unwrap_or(&0) + m_id;
                let t_id_new = self.dimensions.get(&new_t_gr).unwrap_or(&0) + t_id;
                Some(((new_t_gr, t_id_new), ((a_gr, a_id), (new_m_gr, m_id_new))))
            },
        ));

        other.dimensions.iter().for_each(|(gr, other_size)| {
            let new_gr = *gr + grade;
            if new_gr <= limit {
                self.dimensions
                    .entry(new_gr)
                    .and_modify(|size| {
                        *size += *other_size;
                    })
                    .or_insert(*other_size);
            }
        });

        debug_assert!(self.is_correct());
    }


    /// Special care should be taken here,
    /// as this tensor object usually relates to certain graded linear maps
    /// We should expect that direct summing the underlying vector spaces creates the correct new tensor object
    ///
    /// After careful thinking, this direct sum is not dependent on the non-determinism of the hashmap
    pub fn direct_sum(&mut self, other: &Self, self_space_dimensions: &HashMap<G, usize>) {
        // TODO: Consider using parallel here for bigger coalgebras ?
        self.construct
            .extend(other.construct.iter().map(|(&(m_gr, m_id), maps)| {
                let m_id_new = self_space_dimensions.get(&m_gr).unwrap_or(&0) + m_id;

                let new_map = maps
                    .iter()
                    .map(|((a_gr, a_id), (t_gr, t_id))| {
                        let t_id_new = self.dimensions.get(t_gr).unwrap_or(&0) + t_id;
                        ((*a_gr, *a_id), (*t_gr, t_id_new))
                    })
                    .collect();

                ((m_gr, m_id_new), new_map)
            }));

        self.deconstruct.extend(other.deconstruct.iter().map(
            |(&(t_gr, t_id), &((a_gr, a_id), (m_gr, m_id)))| {
                let m_id_new = self_space_dimensions.get(&m_gr).unwrap_or(&0) + m_id;
                let t_id_new = self.dimensions.get(&t_gr).unwrap_or(&0) + t_id;
                ((t_gr, t_id_new), ((a_gr, a_id), (m_gr, m_id_new)))
            },
        ));

        other.dimensions.iter().for_each(|(gr, other_size)| {
            self.dimensions
                .entry(*gr)
                .and_modify(|size| {
                    *size += *other_size;
                })
                .or_insert(*other_size);
        });

        debug_assert!(self.is_correct());
    }

    pub fn is_correct(&self) -> bool {
        let decon_wrong =
            self.deconstruct
                .iter()
                .any(|(t_id, (a_id, m_id))| match self.construct.get(m_id) {
                    Some(map) => match map.get(a_id) {
                        Some(comp_t_id) => t_id != comp_t_id,
                        None => true,
                    },
                    None => true,
                });
        let con_wrong = self.construct.iter().any(|(m_id, map)| {
            map.iter()
                .any(|(a_id, t_id)| match self.deconstruct.get(t_id) {
                    Some((comp_a_id, comp_m_id)) => (comp_a_id != a_id) || (comp_m_id != m_id),
                    None => true,
                })
        });

        let mut dims = HashMap::default();

        let mut found: HashMap<G, Vec<bool>, RandomState> = self
            .dimensions
            .iter()
            .map(|(&gr, &size)| (gr, (0..size).map(|_| false).collect()))
            .collect();
        for ((t_gr, t_id), _) in self.deconstruct.iter() {
            found.get_mut(t_gr).unwrap()[*t_id] = true;
            dims.entry(*t_gr).and_modify(|x| *x += 1).or_insert(1);
        }

        for v in found {
            debug_assert!(v.1.iter().all(|x| *x));
        }

        let dim_wrong = dims != self.dimensions;

        !decon_wrong && !con_wrong && !dim_wrong
    }
}
