use std::sync::Arc;

use crate::linalg::graded::{BasisElement, BasisIndex, Grading};

pub trait Comodule<G: Grading> {
    type Element: BasisElement;
    type Coalgebra;

    fn zero_comodule(coalgebra: Arc<Self::Coalgebra>) -> Self;

    // This is not the correct type yet
    fn get_generators(&self) -> Vec<(usize, G, Option<String>)>;

    fn fp_comodule(coalgebra: Arc<Self::Coalgebra>) -> Self;

    fn direct_sum(&mut self, other: &mut Self);

    fn cofree_comodule(coalgebra: Arc<Self::Coalgebra>, index: usize, grade: G, limit: G) -> Self;
}

pub trait ComoduleMorphism<G: Grading, M: Comodule<G>> {
    fn cokernel(&self) -> Self;
    fn inject_codomain_to_cofree(&self, limit: G) -> Self; // Question: Shouldn't 'codomain' be 'cokernel'/'comodule'?

    fn zero_morphism(comodule: Arc<M>) -> Self;

    // domain l == codomain r, l \circ r
    fn compose(l: Self, r: Self) -> Self;

    fn get_codomain(&self) -> Arc<M>;

    fn get_structure_lines(&self) -> Vec<(BasisIndex<G>, BasisIndex<G>, usize)>; // Question: Should usize here be Field?
}

pub trait Tensor<G: Grading> {
    fn tensor_to_base();
    fn base_to_tensor();

    fn get_dimension(&self, grading: &G) -> usize;
}

// pub struct BasicComodule<M: Module, MM: ModuleMorphism, H: HopfAlgebra> {
//     underlying_module: M,
//     underlying_hopf_algebra: H,

//     coaction: MM,
// }

// impl<M: Module, MM: ModuleMorphism, H: HopfAlgebra> Comodule for BasicComodule<M, MM, H> {

//     fn get_underlying_module(&self) -> &impl Module {
//         &self.underlying_module
//     }

//     fn get_hopf_algebra(&self) -> &impl HopfAlgebra {
//         &self.underlying_hopf_algebra
//     }

//     fn get_cogenerating_module(&self) -> impl ModuleMorphism {
//         todo!()
//     }

//     fn cofree_comodule(hopf: impl HopfAlgebra, module: impl Module) -> impl Comodule {
//         todo!()
//     }

// }

// pub struct BasicComoduleMorphism<M: Module, MM: ModuleMorphism, H: HopfAlgebra> {
//     domain: BasicComodule<M, MM, H>,
//     codomain: BasicComodule<M, MM, H>,

//     morphism: MM
// }

// impl <M: Module, MM: ModuleMorphism, H: HopfAlgebra> ComoduleMorphism for BasicComoduleMorphism<M, MM, H> {
//     fn get_underlying_morphism(&self) -> &impl ModuleMorphism {
//         &self.morphism
//     }

//     fn get_domain(&self) -> &impl Comodule {
//         &self.domain
//     }

//     fn get_codomain(&self) -> &impl Comodule {
//         &self.codomain
//     }

//     fn compute_cokernel(&self) -> impl ComoduleMorphism {
//         let mod_coker = &self.morphism.get_cokernel();

//         //    f      q
//         // C  ->  D  ->  K
//         // |      |      |
//         // v      v      v
//         //CxC -> DxD -> KxK

//         // let k be a basis element in K
//         // q is surjective, thus we get an element b in D which maps to k
//         todo!()

//     }

//     fn get_zero_morphism_to(comod: BasicComodule<M, MM, H>) -> Self {
//         todo!()
//     }
// }
