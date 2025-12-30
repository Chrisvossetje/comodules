use ahash::HashMap;
use deepsize::DeepSizeOf;
use serde::{Deserialize, Serialize};

use crate::grading::{grading::Grading, tensor::TensorMap, unigrading::UniGrading};

// TODO :
// impl<G: Grading> GradedModule<G> {
//     pub fn dimensions(&self) -> HashMap<G, usize> {
//         self.0.iter().map(|(g, v)| (*g, v.len())).collect()
//     }

//     pub fn dimension_in_grade(&self, grade: &G) -> usize {
//         self.0.get(grade).map(|x| x.len()).unwrap_or(0)
//     }

//     pub fn generate_tensor_as_module(
//         &self,
//         coalgebra: &Self,
//         tensor: &TensorMap<G>,
//     ) -> GradedModule<G> {
//         let mut map = HashMap::default();
//         for (&gr, dim) in &tensor.dimensions {
//             // TODO : Unigrading ?
//             map.insert(gr, vec![(UniGrading(0), None); *dim]);
//         }
//         for (&(t_gr, t_id), &((a_gr, a_id), (m_gr, m_id))) in &tensor.deconstruct {
//             let module = &mut map.get_mut(&t_gr).unwrap()[t_id as usize];
//             let a_unigrade = coalgebra.0.get(&a_gr).unwrap()[a_id as usize].0;
//             let (m_unigrade, m_quotient) = self.0.get(&m_gr).unwrap()[m_id as usize];
//             module.0 = a_unigrade + m_unigrade;
//             module.1 = m_quotient;
//         }

//         GradedModule(map)
//     }
// }
