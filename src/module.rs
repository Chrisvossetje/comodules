use crate::{fp::Field, graded::{GradedBasisElement, GradedVectorSpace, Grading}};


// EXAMPLE OF A kt MODULE 

#[derive(Debug, Clone)]
pub struct ktBasisElement {

}

pub struct ktModule {
    space: GradedVectorSpace<(i32, i32), ktBasisElement>,
    t: GradedLinearMap<UniGrading, F2, FieldMatrix<F2>>
}




pub trait Module<G: Grading, F: Field> {
    fn get_vector_space(&self) -> &impl GradedVectorSpace<G, F>;
    fn get_basis(&self) -> &Vec<&impl GradedBasisElement<G>> { self.get_vector_space().get_basis() }
}

pub trait TensorModule<G: Grading, F: Field> : Module<G, F> {
    fn get_left_module(&self) -> &impl Module<G, F>;
    fn get_right_module(&self) -> &impl Module<G, F>;

    fn compose_basis_element(&self, left: &impl GradedBasisElement<G>, right: &impl GradedBasisElement<G>) -> &impl GradedBasisElement<G>;
    fn decompose_basis_element(&self, composed: &impl GradedBasisElement<G>) -> (&impl GradedBasisElement<G>, &impl GradedBasisElement<G>);
}

pub trait ModuleMorphism<G: Grading, F: Field> {
    fn get_domain(&self) -> impl Module<G, F>;
    fn get_codomain(&self) -> impl Module<G, F>;

    fn get_cokernel(&self) -> impl ModuleMorphism<G, F>;
    fn transpose(&self) -> Self;
}