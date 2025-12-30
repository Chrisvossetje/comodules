use std::fmt::Debug;

use algebra::{abelian::Abelian, field::Field, ring::CRing};
use deepsize::DeepSizeOf;
use serde::{Deserialize, Serialize};

use crate::{grading::grading::Grading, types::CoalgebraIndex};

pub trait Coalgebra<G: Grading>: DeepSizeOf + Clone + Sync + Send + Debug {
    type BaseField: Field;
    type BaseRing: CRing;
    type RingMorph: Abelian<Self::BaseRing>;

    type Comod: Comodule<G, Self>;
    type CofMod: CofreeComodule<G, Self>;
    type ComodMorph: ComoduleMorphism<G, Self>;

    fn size_in_degree(&self, g: G) -> usize;
    fn coaction(
        &self,
        i: CoalgebraIndex<G>,
    ) -> &[(CoalgebraIndex<G>, CoalgebraIndex<G>, Self::BaseRing)];
    fn basering_comodule(&self, shift: G) -> Self::Comod;
    fn cofree_comodule(&self, index: usize, shift: G, limit: G, generator: <Self::RingMorph as Abelian<Self::BaseRing>>::Generator) -> Self::CofMod;
}

pub trait Comodule<G: Grading, C: Coalgebra<G>>: DeepSizeOf + Clone + Send + Sync {
    fn zero_comodule() -> Self;
}

pub trait CofreeComodule<G: Grading, C: Coalgebra<G>>: DeepSizeOf + Clone {
    // This is not the correct type yet ???
    fn zero_comodule() -> Self;
    fn get_generators(&self) -> Vec<(usize, G, Option<String>)>;
    fn direct_sum(&mut self, other: &mut Self);
}

pub trait ComoduleMorphism<G: Grading, C: Coalgebra<G>>: Sized + DeepSizeOf + Send + Sync {
    fn cokernel(&self, coalgebra: &C, codomain: &C::CofMod) -> (Self, C::Comod);

    fn inject_codomain_to_cofree(coalgebra: &C, comodule: &C::Comod, limit: G)
    -> (Self, C::CofMod);

    fn zero_morphism(comodule: &C::Comod) -> Self;

    // domain l == codomain r, l \circ r
    fn compose(l: &Self, r: &Self, codomain: &C::CofMod) -> Self;
}
