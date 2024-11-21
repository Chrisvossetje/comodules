use crate::fp::PrimeField;

pub trait Matrix<F: PrimeField> {
    fn get_matrix(&self) -> Vec<Vec<F>>;
}


pub trait GradedBasisElement {

}

pub trait GradedLinearTransformation<F: PrimeField> {
    fn get_domain(&self) -> &impl GradedVectorSpace<F>;
    fn get_codomain(&self) -> &impl GradedVectorSpace<F>;

    fn get_matrix(&self) -> &impl Matrix<F>;
}

pub trait GradedVectorSpace<F: PrimeField> {
    fn get_basis(&self) -> &Vec<impl GradedBasisElement>;
}

