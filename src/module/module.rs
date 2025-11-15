use ahash::HashMap;
use serde::{Deserialize, Serialize};

use crate::{
    basiselement::BasisElement,
    grading::{Grading, UniGrading}, tensor::Tensor
};


#[allow(type_alias_bounds)]
type Module<B: BasisElement> = Vec<(B, UniGrading, Option<usize>)>;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Default)]
pub struct GradedModule<G: Grading, B: BasisElement>(pub HashMap<G, Module<B>>);



impl<G: Grading, B: BasisElement> GradedModule<G, B> {
    pub fn dimensions(&self) -> HashMap<G, usize> {
        self.0.iter().map(|(g, v)| (*g,v.len())).collect()
    }

    pub fn dimension_in_grade(&self, grade: &G) -> usize {
        self.0.get(grade).map(|x| x.len()).unwrap_or(0)
    }

    pub fn generate_tensor_as_module(&self, coalgebra: &Self, tensor: &Tensor<G>) -> GradedModule<G, B> {
        let mut map = HashMap::default();
        for (&gr, dim) in &tensor.dimensions {
            map.insert(gr, vec![(B::default(), UniGrading(0), None); *dim]);
        }
        for (&(t_gr, t_id), &((a_gr, a_id),(m_gr, m_id))) in &tensor.deconstruct {
            let module = &mut map.get_mut(&t_gr).unwrap()[t_id];
            let a_unigrade = coalgebra.0.get(&a_gr).unwrap()[a_id].1;
            let (_, m_unigrade, m_quotient) = self.0.get(&m_gr).unwrap()[m_id];
            module.1 = a_unigrade + m_unigrade;
            module.2 = m_quotient;
        }

        GradedModule(map)
    }
}
