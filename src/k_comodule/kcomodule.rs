use algebra::{abelian::Abelian, field::Field, matrices::flat_matrix::FlatMatrix, matrix::Matrix, ring::CRing};
use deepsize::DeepSizeOf;

use crate::{
    grading::{grading::Grading, tensor::TensorMap},
    k_comodule::graded_space::{GradedLinearMap, GradedVectorSpace},
    traits::{Coalgebra, CofreeComodule, Comodule},
    types::CoalgebraIndex,
};

use super::kcoalgebra::kCoalgebra;

#[derive(Clone, PartialEq, DeepSizeOf)]
#[allow(non_camel_case_types)]
pub struct kComodule<G: Grading, C: Coalgebra<G>> where C::BaseRing: Field {
    pub space: GradedVectorSpace<G, <C::RingMorph as Abelian<C::BaseRing>>::Generator>,
    pub coaction: GradedLinearMap<G, C::BaseRing, C::RingMorph>,
    pub tensor: TensorMap<G>,
}

#[derive(Debug, Clone, PartialEq, DeepSizeOf)]
#[allow(non_camel_case_types)]
pub struct kCofreeComodule<G: Grading, C: Coalgebra<G>> where C::BaseRing: Field  {
    pub space: GradedVectorSpace<G, ((CoalgebraIndex<G>, u16), <C::RingMorph as Abelian<C::BaseRing>>::Generator)>,
}

impl<G: Grading, C: Coalgebra<G>> std::fmt::Debug for kComodule<G, C> where C::BaseRing: Field {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.space.0)
    }
}

impl<G: Grading, C: Coalgebra<G>> kComodule<G, C> where C::BaseRing: Field {
    pub fn verify(&self) -> bool {
        for (&(m_gr, m_id), map) in &self.tensor.construct {
            let &(t_gr, t_id) = map.get(&(G::zero(), 0)).unwrap();
            if t_gr != m_gr {
                return false;
            };

            let val = self
                .coaction
                .maps
                .get(&t_gr)
                .unwrap()
                .get(m_id as usize, t_id as usize);
            if val != C::BaseRing::one() {
                return false;
            };
        }
        true
    }

    pub fn new(
        space: GradedVectorSpace<G, <C::RingMorph as Abelian<C::BaseRing>>::Generator>,
        coaction: GradedLinearMap<G, C::BaseRing, C::RingMorph>,
        tensor: TensorMap<G>,
    ) -> Self {
        let com = Self {
            space,
            coaction,
            tensor,
        };
        debug_assert!(com.verify());
        com
    }
}

impl<G: Grading, C: Coalgebra<G>> CofreeComodule<G, C> 
    for kCofreeComodule<G, C> where C::BaseRing: Field 
{
    fn get_generators(&self) -> Vec<(usize, G, Option<String>)> {
        self.space
            .0
            .iter()
            .flat_map(|(k, v)| {
                v.iter().filter_map(|b| {
                    if b.0.0.0 == G::zero() {
                        Some((b.0.1 as usize, *k, None))
                    } else {
                        None
                    }
                })
            })
            .collect()
    }

    fn direct_sum(&mut self, other: &mut Self) {
        other.space.0.iter_mut().for_each(|(g, other_els)| {
            self.space
                .0
                .entry(*g)
                .and_modify(|self_els| {
                    self_els.extend(other_els.drain(0..));
                })
                .or_insert(other_els.drain(0..).collect());
        });
    }

    fn zero_comodule() -> Self {
        Self {
            space: GradedVectorSpace::new(),
        }
    }
}

impl<G: Grading, C: Coalgebra<G>> Comodule<G, kCoalgebra<G, C::BaseRing, C::RingMorph>> for kComodule<G, C> where C::BaseRing: Field {
    fn zero_comodule() -> Self {
        Self {
            space: GradedVectorSpace::new(),
            coaction: GradedLinearMap::empty(),
            tensor: TensorMap::default(),
        }
    }
}
