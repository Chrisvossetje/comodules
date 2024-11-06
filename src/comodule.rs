use crate::{hopfalgebra::{HopfAlgebra}, module::{Module, ModuleMorphism}};

pub trait Comodule {
    // Applies forgetful functor
    fn get_underlying_module(&self) -> impl Module;
    
    fn get_hopf_algebra(&self) -> impl HopfAlgebra;

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

    fn get_zero_morphism_to(&self) -> impl ComoduleMorphism;
}

pub trait ComoduleMorphism {
    // Applies forgetful functor
    fn get_underlying_morphism(&self) -> impl ModuleMorphism;

    fn get_domain(&self) -> impl Comodule;
    fn get_codomain(&self) -> impl Comodule;

    fn get_cokernel(&self) -> impl ComoduleMorphism;

    // I couldn't put this function here because it started a fight with the compiler.
    // fn get_zero_morphism_to(comod: impl Comodule) -> impl ComoduleMorphism;
}