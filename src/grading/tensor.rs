use ahash::HashMap;
use deepsize::DeepSizeOf;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use std::fmt::Debug;

use crate::{
    grading::grading::Grading,
    k_comodule::graded_space::GradedVectorSpace,
    types::{CoalgebraIndex, CoalgebraIndexType, ComoduleIndex, ComoduleIndexType},
};

pub type TensorConstruct<G> =
    HashMap<ComoduleIndex<G>, HashMap<CoalgebraIndex<G>, ComoduleIndex<G>>>;
pub type TensorDeconstruct<G> = HashMap<ComoduleIndex<G>, (CoalgebraIndex<G>, ComoduleIndex<G>)>;
pub type TensorDimension<G> = HashMap<G, usize>;

pub trait ObjectGenerator<G: Grading> {
    fn contains_grade(&self, grade: &G) -> bool;
    fn els(&self, grade: &G) -> usize;
    fn sorted_els(&self) -> Vec<(G, usize)>;
}

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Default, DeepSizeOf)]
pub struct TensorMap<G: Grading> {
    // # Module Grade + Index -> Algebra Grading + index -> Tensor Grading + index
    pub construct: TensorConstruct<G>,

    // # Tensor Grade + Index -> (Algebra Grading + index, Module Grading + index)
    pub deconstruct: TensorDeconstruct<G>,
    pub dimensions: TensorDimension<G>,
}

impl<G: Grading, B> ObjectGenerator<G> for GradedVectorSpace<G, B> {
    fn contains_grade(&self, grade: &G) -> bool {
        self.0.contains_key(grade)
    }

    fn els(&self, grade: &G) -> usize {
        self.0.get(grade).map_or(0, |x| x.len())
    }

    fn sorted_els(&self) -> Vec<(G, usize)> {
        self.0
            .iter()
            .sorted_by_key(|(g, _)| **g)
            .map(|(&g, l)| (g, l.len()))
            .collect()
    }
}

impl<G: Grading> TensorMap<G> {
    pub fn default() -> Self {
        Self {
            construct: Default::default(),
            deconstruct: Default::default(),
            dimensions: Default::default(),
        }
    }

    pub fn get_dimension(&self, grading: &G) -> usize {
        *self.dimensions.get(grading).unwrap_or(&0)
    }
}

impl<G: Grading> TensorMap<G> {
    pub fn generate<G1: ObjectGenerator<G>, G2: ObjectGenerator<G>>(left: &G1, right: &G2) -> Self {
        let mut construct = HashMap::default();
        let mut deconstruct = HashMap::default();
        let mut dimensions = HashMap::default();

        // This sorted is important for consistency !
        for (l_grade, l_elements) in left.sorted_els() {
            for l_id in 0..l_elements as CoalgebraIndexType {
                for (r_grade, r_elements) in right.sorted_els() {
                    let t_grade = l_grade + r_grade;

                    if !right.contains_grade(&t_grade) {
                        continue;
                    }

                    for r_id in 0..r_elements as ComoduleIndexType {
                        let t_id = dimensions.entry(t_grade).or_insert(0 as usize);

                        construct
                            .entry((r_grade, r_id))
                            .or_insert(HashMap::default())
                            .insert((l_grade, l_id), (t_grade, *t_id as ComoduleIndexType));

                        deconstruct.insert(
                            (t_grade, *t_id as ComoduleIndexType),
                            ((l_grade, l_id), (r_grade, r_id)),
                        );
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

    pub fn add_and_restrict(&self, add: G, limit: G) -> Self {
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
        let tens = Self {
            construct: cons,
            deconstruct: decon,
            dimensions: dims,
        };
        debug_assert!(tens.is_correct());
        tens
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

        let mut found: HashMap<G, Vec<bool>> = self
            .dimensions
            .iter()
            .map(|(&gr, &size)| (gr, (0..size).map(|_| false).collect()))
            .collect();
        for ((t_gr, t_id), _) in self.deconstruct.iter() {
            found.get_mut(t_gr).unwrap()[*t_id as usize] = true;
            dims.entry(*t_gr).and_modify(|x| *x += 1).or_insert(1);
        }

        for v in found {
            debug_assert!(v.1.iter().all(|x| *x));
        }

        let dim_wrong = dims != self.dimensions;

        !decon_wrong && !con_wrong && !dim_wrong
    }
}
