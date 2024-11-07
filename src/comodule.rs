use crate::{hopfalgebra::HopfAlgebra, module::{Module, ModuleMorphism, TensorModule}};

pub trait Comodule {
    // Applies forgetful functor
    fn get_underlying_module(&self) -> &impl Module;
    
    fn get_hopf_algebra(&self) -> &impl HopfAlgebra;

    // returns the smallest submodule of which 'cogenerates' the entire comodule.
    // Let \Delta(x) = \sum_i (a_i \otimes y_i),
    // if the coefficient of y is non-zero, we say that x 'cogenerates' y
    // We define cogeneration to be the smallest transitive relation with the above property.
    fn get_cogenerating_module(&self) -> impl Module;
    // !!! Note: We want morphisms instead of modules, so we eventually get the injection to the cofree module
    // fn get_cogenerating_module(&self) -> impl ModuleMorphism;

    // Produces the tensor product of given the Hopf algebra A with the given module V.
    // This therefore returns the A-comodule: A \otimes V
    fn cofree_comodule(hopf: impl HopfAlgebra, module: impl Module) -> impl Comodule;
    // !!! Note: We want to return the injection from the coker to the cofree comodule, not just the comod itself.

    fn create_tensor_product(left: &Self, right: &Self) -> Self;
}

pub trait ComoduleMorphism {
    // Applies forgetful functor
    fn get_underlying_morphism(&self) -> &impl ModuleMorphism;

    fn get_domain(&self) -> &impl Comodule;
    fn get_codomain(&self) -> &impl Comodule;

    fn compute_cokernel(&self) -> impl ComoduleMorphism;

    fn get_zero_morphism_to(comod: impl Comodule) -> Self;
}

pub struct BasicComodule<M: Module, MM: ModuleMorphism, H: HopfAlgebra> {
    underlying_module: M,
    underlying_hopf_algebra: H,

    coaction: MM,
}


impl<M: Module, MM: ModuleMorphism, H: HopfAlgebra> Comodule for BasicComodule<M, MM, H> {

    fn get_underlying_module(&self) -> &impl Module {
        &self.underlying_module
    }

    fn get_hopf_algebra(&self) -> &impl HopfAlgebra {
        &self.underlying_hopf_algebra
    }

    fn get_cogenerating_module(&self) -> impl Module {
        todo!()
    }

    fn cofree_comodule(hopf: impl HopfAlgebra, module: impl Module) -> impl Comodule {
        todo!()
    }

    fn create_tensor_product(left: &Self, right: &Self) -> Self {
        todo!()
    }
}


pub struct BasicComoduleMorphism<M: Module, MM: ModuleMorphism, H: HopfAlgebra> {
    domain: BasicComodule<M, MM, H>,
    codomain: BasicComodule<M, MM, H>,

    morphism: MM
}

impl <M: Module, MM: ModuleMorphism, H: HopfAlgebra> ComoduleMorphism for BasicComoduleMorphism<M, MM, H> {
    fn get_underlying_morphism(&self) -> &impl ModuleMorphism {
        &self.morphism
    }

    fn get_domain(&self) -> &impl Comodule {
        &self.domain
    }

    fn get_codomain(&self) -> &impl Comodule {
        &self.codomain
    }

    fn compute_cokernel(&self) -> impl ComoduleMorphism {
        todo!()
    }
    
    fn get_zero_morphism_to(comod: BasicComodule<M, MM, H>) -> Self {
        todo!()
    }
}
