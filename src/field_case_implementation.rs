// use crate::{coalgebra::Coalgebra, comodule::{Comodule, ComoduleMorphism}, field::Field, graded::{GradedLinearMap, GradedVectorSpace, Grading}, matrix::Matrix, module::{Module, ModuleMorphism}};

// //
// // Module
// //

// #[derive(Clone)]
// pub struct FieldModule<G: Grading, F: Field, M: Matrix<F>> {
//     object: GradedVectorSpace<G>,
//     action: GradedLinearMap<G, F, M>
// }

// impl<G: Grading, F: Field, M: Matrix<F>> Module<G, F, M> for FieldModule<G, F, M>{
//     fn get_vector_space(&self) -> &GradedVectorSpace<G> {
//         &self.object
//     }
    
//     fn get_action(&self) -> &GradedLinearMap<G, F, M> {
//         &self.action
//     }
// }

// #[derive(Clone)]
// pub struct FieldModuleMorphism<G: Grading, F: Field, M: Matrix<F>> {
//     map: GradedLinearMap<G, F, M>,
//     domain: FieldModule<G, F, M>,
//     codomain: FieldModule<G, F, M>
// }

// impl<G: Grading, F: Field, M: Matrix<F>> ModuleMorphism<G, F, M> for FieldModuleMorphism<G, F, M> {
//     fn get_domain(&self) -> &impl Module<G, F, M> {
//         &self.domain
//     }

//     fn get_codomain(&self) -> &impl Module<G, F, M> {
//         &self.codomain
//     }

//     fn get_kernel(&self) -> Self {
//         todo!()
//     }

//     fn get_cokernel(&self) -> Self {
//         todo!()
//     }

//     fn transpose(&self) -> Self {
//         FieldModuleMorphism {
//             map: self.map.transpose(),
//             domain: self.codomain.clone(),
//             codomain: self.domain.clone()
//         }
//     }
// }

// //
// // Coalgebra
// //
// #[derive(Clone)]
// pub struct FieldCoalgebra<G: Grading, F: Field, M: Matrix<F>> {
//     module: FieldModule<G, F, M>,
//     coaction: FieldModuleMorphism<G, F, M>
// }

// impl<G: Grading, F: Field, M: Matrix<F>> Coalgebra<G, F, M> for FieldCoalgebra<G, F, M> {
//     fn get_module(&self) -> &impl Module<G, F, M> {
//         &self.module
//     }

//     fn get_comultiplication(&self) -> &impl ModuleMorphism<G, F, M> {
//         &self.coaction
//     }
// }


// //
// // Comodule
// //
// #[derive(Clone)]
// pub struct FieldComodule<G: Grading, F: Field, M: Matrix<F>> {
//     module: FieldModule<G, F, M>,
//     coalgebra: FieldCoalgebra<G, F, M>,
//     coaction: FieldModuleMorphism<G, F, M>
// }

// // impl<G: Grading, F: Field, M: Matrix<F>> Comodule<G, F, M> for FieldComodule<G, F, M> {
// //     fn get_underlying_module(&self) -> &impl Module<G, F, M> {
// //         &self.module
// //     }
    
// //     fn get_coalgebra(&self) -> &impl Coalgebra<G, F, M> {
// //         &self.coalgebra
// //     }
    
// //     fn get_cogenerating_module(&self) -> impl ModuleMorphism<G, F, M> {
// //         unimplemented!();
// //         self.coaction
// //     }
    
// //     fn cofree_comodule(hopf: impl Coalgebra<G, F, M>, module: impl Module<G, F, M>) -> Self {
// //         todo!()
// //     }
    
// //     fn create_tensor_product(left: &Self, right: &Self) -> Self {
// //         todo!()
// //     }

// // }


// // #[derive(Clone)]
// // pub struct FieldComoduleMorphism<'a, G: Grading, F: Field, M: Matrix<F>> {
// //     map: FieldModuleMorphism<G, F, M>,
// //     domain: &'a FieldComodule<G, F, M>,
// //     codomain: &'a FieldComodule<G, F, M>
// // }

// // impl<'a, G: Grading, F: Field, M: Matrix<F>> ComoduleMorphism<G, F, M> for FieldComoduleMorphism<'a, G, F, M> {
// //     fn get_underlying_morphism(&self) -> &impl ModuleMorphism<G, F, M> {
// //         &self.map
// //     }

// //     fn get_domain(&self) -> &impl Comodule<G, F, M> {
// //         &self.domain
// //     }

// //     fn get_codomain(&self) -> &impl Comodule<G, F, M> {
// //         &self.codomain
// //     }

// //     fn compute_cokernel(&self) -> Self {
// //         todo!()
// //     }

// //     fn get_zero_morphism_to(comod: impl Comodule<G, F, M>) -> Self {
// //         todo!()
// //     }
// // }



