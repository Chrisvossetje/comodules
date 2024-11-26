use crate::{field::Field, graded::Grading, matrix::Matrix, module::{Module, ModuleMorphism}};

pub trait Coalgebra<G: Grading, F: Field, M: Matrix<F>> : Clone {
    fn get_module(&self) -> &impl Module<G, F, M>;
    fn get_comultiplication(&self) -> &impl ModuleMorphism<G, F, M>;
}