use std::sync::Arc;
use std::fmt::Debug;

use ahash::HashMap;
use algebra::{abelian::Abelian, ring::CRing};
use deepsize::DeepSizeOf;

use crate::{graded_space::BasisIndex, grading::{Grading, OrderedGrading}, tensor::ObjectGenerator};

pub trait Coalgebra<G: Grading>: DeepSizeOf + Clone + Sync + Send + Debug {
    type BaseRing: CRing;
    type RingMorph: Abelian<Self::BaseRing>;

    type Comod: Comodule<G,Self>;
    type CofMod: CofreeComodule<G,Self>;
    type ComodMorph: ComoduleMorphism<G, Self>;

    fn size_in_degree(&self, g: G) -> usize;
    // TODO: to slice
    fn coaction(&self, i: BasisIndex<G>) -> &[(BasisIndex<G>, BasisIndex<G>, Self::BaseRing)];
}

pub trait Comodule<G: Grading, C: Coalgebra<G>>: DeepSizeOf + Clone + Send + Sync {
    fn fp_comodule(coalgebra: &C, degree: G) -> Self;
}

pub trait CofreeComodule<G: Grading, C: Coalgebra<G>>: DeepSizeOf + Clone {
    type Generator: Default;

    // This is not the correct type yet ???
    fn zero_comodule() -> Self;
    fn get_generators(&self) -> Vec<(usize, G, Option<String>)>;
    fn direct_sum(&mut self, other: &mut Self);
    fn cofree_comodule(
        coalgebra: &C,
        index: usize,
        grade: G,
        limit: G,
        generator: Self::Generator,
    ) -> Self;
}

// pub trait CofreeComodule<G: Grading>: Sized {

// }

pub trait ComoduleMorphism<G: Grading, C: Coalgebra<G>>: Sized + DeepSizeOf + Send + Sync {
    type BaseRing: CRing;

    fn cokernel(&self, coalgebra: &C, codomain: &C::CofMod) -> (Self, C::Comod);

    fn inject_codomain_to_cofree(
        coalgebra: &C,
        comodule: &C::Comod,
        limit: G,
    ) -> (Self, C::CofMod)
    where
        G: OrderedGrading;

    fn zero_morphism(comodule: &C::Comod) -> Self;

    // domain l == codomain r, l \circ r
    fn compose(l: &Self, r: &Self, codomain: &C::CofMod) -> Self;

    fn direct_sum(
        a: &mut Self,
        b: &mut Self,
        a_codom: &mut C::CofMod,
        b_codom: &mut C::CofMod,
    );

    // fn get_codomain(&self) -> Arc<M>;

    /// (s, gen_index) uniquely defines a generator of Ext
    /// in a specific morphism we only need to know its gen_index
    /// in the resolution we add the s
    /// (from_dot, to_dot, value, line_type)
    fn get_structure_lines(
        &self,
        coalgebra: &C,
        domain: &C::CofMod,
        codomain: &C::CofMod,
    ) -> Vec<((usize, G, usize), (usize, G, usize), Self::BaseRing, String)>;
}
