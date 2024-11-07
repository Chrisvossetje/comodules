use crate::linearalgebra::{GradedBasisElement, GradedVectorSpace};

pub trait Module {
    fn get_vector_space(&self) -> &impl GradedVectorSpace;
    fn get_basis(&self) -> &Vec<impl GradedBasisElement> { self.get_vector_space().get_basis() }
}

pub trait TensorModule: Module {
    fn get_left_module(&self) -> &impl Module;
    fn get_right_module(&self) -> &impl Module;

    fn compose_basis_element(&self, left: &impl GradedBasisElement, right: &impl GradedBasisElement) -> &impl GradedBasisElement;
    fn decompose_basis_element(&self, composed: &impl GradedBasisElement) -> (&impl GradedBasisElement, &impl GradedBasisElement);
}

pub trait ModuleMorphism {
    fn get_domain(&self) -> impl Module;
    fn get_codomain(&self) -> impl Module;

    fn get_cokernel(&self) -> impl ModuleMorphism;
    fn transpose(&self) -> Self;
}