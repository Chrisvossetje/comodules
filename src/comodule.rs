use crate::{coalgebra::Coalgebra, field::Field, graded::Grading, matrix::Matrix, module::{Module, ModuleMorphism, TensorModule}};

pub trait Comodule<G: Grading, F: Field, M: Matrix<F>> {
    // Applies forgetful functor
    // haha, "fancy math boy"
    fn get_underlying_module(&self) -> &impl Module<G, F, M>;
    
    fn get_coalgebra(&self) -> &impl Coalgebra<G, F, M>;

    // returns the smallest submodule of which 'cogenerates' the entire comodule.
    // Let \Delta(x) = \sum_i (a_i \otimes y_i),
    // if the coefficient of y is non-zero, we say that x 'cogenerates' y
    // We define cogeneration to be the smallest transitive relation with the above property.
    fn get_cogenerating_module(&self) -> impl ModuleMorphism<G, F, M>;
    // !!! Note: We want morphisms instead of modules, so we eventually get the injection to the cofree module
    // fn get_cogenerating_module(&self) -> impl ModuleMorphism;

    // Produces the tensor product of given the Hopf algebra A with the given module V.
    // This therefore returns the A-comodule: A \otimes V
    fn cofree_comodule(hopf: impl Coalgebra<G, F, M>, module: impl Module<G, F, M>) -> Self;
    // !!! Note: We want to return the injection from the coker to the cofree comodule, not just the comod itself.

    fn create_tensor_product(left: &Self, right: &Self) -> Self;
}

pub trait ComoduleMorphism<G: Grading, F: Field, M: Matrix<F>> {
    // Applies forgetful functor
    fn get_underlying_morphism(&self) -> &impl ModuleMorphism<G, F, M>;

    fn get_domain(&self) -> &impl Comodule<G, F, M>;
    fn get_codomain(&self) -> &impl Comodule<G, F, M>;

    fn compute_cokernel(&self) -> Self;

    fn get_zero_morphism_to(comod: impl Comodule<G, F, M>) -> Self; 
}


pub trait SimpleComodule {

}

pub trait SimpleComoduleMorphism<M: SimpleComodule> {
    fn cokernel(&self) -> Self;
    fn injection_codomain_to_cofree(&self) -> Self;

    fn zero_morphism(comodule: M) -> Self;

    // codomain r == codomain l, l \circ r
    fn compose(l: Self, r: Self) -> Self;
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
