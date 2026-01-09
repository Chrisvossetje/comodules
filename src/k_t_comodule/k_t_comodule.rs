use ahash::HashMap;
use algebra::{abelian::Abelian, field::Field, rings::univariate_polynomial_ring::UniPolRing};
use deepsize::DeepSizeOf;

use crate::{
    grading::{grading::Grading, tensor::TensorMap},
    k_comodule::graded_space::GradedVectorSpace,
    k_t_comodule::{graded_module::GradedktFieldMap, k_t_coalgebra::ktCoalgebra},
    traits::{Coalgebra, CofreeComodule, Comodule},
    types::CoalgebraIndex,
};

#[derive(Clone, PartialEq, DeepSizeOf)]
#[allow(non_camel_case_types)]
pub struct ktComodule<G: Grading, F: Field, M: Abelian<UniPolRing<F>>> {
    pub space: GradedVectorSpace<G, M::Generator>,
    pub coaction: GradedktFieldMap<G, F, M>,
    pub tensor: TensorMap<G>,
}

#[derive(Debug, Clone, PartialEq, DeepSizeOf)]
#[allow(non_camel_case_types)]
pub struct ktCofreeComodule<G: Grading, C: Coalgebra<G>> {
    pub space: GradedVectorSpace<
        G,
        (
            (CoalgebraIndex<G>, u16),
            <C::RingMorph as Abelian<C::BaseRing>>::Generator,
        ),
    >,
}

impl<G: Grading, C: Coalgebra<G>> ktCofreeComodule<G, C> {
    pub(crate) fn to_module_generators(
        &self,
    ) -> GradedVectorSpace<G, <C::RingMorph as Abelian<C::BaseRing>>::Generator> {
        let a: HashMap<G, Vec<_>> = self
            .space
            .0
            .iter()
            .map(|(g, v)| (*g, v.iter().map(|x| x.1.clone()).collect()))
            .collect();
        GradedVectorSpace(a)
    }
}

impl<G: Grading, F: Field, M: Abelian<UniPolRing<F>>> Comodule<G, ktCoalgebra<G, F, M>>
    for ktComodule<G, F, M>
{
    fn zero_comodule() -> Self {
        Self {
            space: GradedVectorSpace::new(),
            coaction: GradedktFieldMap {
                maps: HashMap::default(),
                _p: std::marker::PhantomData,
            },
            tensor: TensorMap::default(),
        }
    }
}

impl<G: Grading, C: Coalgebra<G>> CofreeComodule<G, C> for ktCofreeComodule<G, C> {
    fn zero_comodule() -> Self {
        Self {
            space: GradedVectorSpace::new(),
        }
    }

    fn get_generators(&self) -> Vec<(usize, G, Option<String>)> {
        self.space
            .0
            .iter()
            .flat_map(|(k, v)| {
                v.iter().filter_map(|b| {
                    if b.0.0.0 == G::zero() {
                        let s = format!("{:?} | {:?}", b.1, b.0);
                        Some((b.0.1 as usize, *k, Some(s)))
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
}
