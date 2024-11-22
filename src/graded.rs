pub trait Grading : Clone + Add<Output=Self> + Sub<Output=Self> + PartialEq + Eq + AddAssign + SubAssign + PartialOrd + Display {}

pub trait GradedBasisElement<G: Grading> {

}

pub trait GradedLinearMap<F: PrimeField, > {
    fn get_domain(&self) -> &impl GradedVectorSpace<F>;
    fn get_codomain(&self) -> &impl GradedVectorSpace<F>;
}

pub trait GradedVectorSpace<F: PrimeField> {
    fn get_basis(&self) -> &Vec<impl GradedBasisElement>;
}



impl Grading for usize {}