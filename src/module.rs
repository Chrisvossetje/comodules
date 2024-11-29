use crate::{field::Field, graded::{GradedLinearMap, GradedVectorSpace, Grading}, matrix::{FieldMatrix, Matrix}};


// EXAMPLE OF A kt MODULE 

// #[derive(Debug, Clone)]
// pub struct ktBasisElement {

// }



// pub struct ktModule {
//     space: GradedVectorSpace<(i32, i32), ktBasisElement>,
//     t: GradedLinearMap<UniGrading, F2, FieldMatrix<F2>>
// }




pub trait Module<G: Grading, F: Field, M: Matrix<F>> : Clone {
    fn get_vector_space(&self) -> &GradedVectorSpace<G>;

    fn get_action(&self) -> &GradedLinearMap<G, F, M>;
}


pub trait TensorModule<G: Grading, F: Field, M: Matrix<F>> : Module<G, F, M> {
    fn get_left_module(&self) -> &impl Module<G, F, M>;
    fn get_right_module(&self) -> &impl Module<G, F, M>;

    // fn compose_basis_element(&self, left: &impl GradedBasisElement<G>, right: &impl GradedBasisElement<G>) -> &impl GradedBasisElement<G>;
    // fn decompose_basis_element(&self, composed: &impl GradedBasisElement<G>) -> (&impl GradedBasisElement<G>, &impl GradedBasisElement<G>);
}

pub trait ModuleMorphism<G: Grading, F: Field, M: Matrix<F>> {
    fn get_domain(&self) -> &impl Module<G, F, M>;
    fn get_codomain(&self) -> &impl Module<G, F, M>;

    fn get_kernel(&self) -> Self;
    fn get_cokernel(&self) -> Self;
    fn transpose(&self) -> Self;
}








pub struct kModule<G: Grading, F: Field> {
    space: GradedVectorSpace<G, G>,
    action: GradedLinearMap<G, F, FieldMatrix<F>>
}

pub struct kModuleMorphism<'a, G: Grading, F: Field> {
    domain: &'a kModule<G, F>,
    codomain: &'a kModule<G, F>,

    map: GradedLinearMap<G, F, FieldMatrix<F>>
}