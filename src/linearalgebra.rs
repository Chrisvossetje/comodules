pub trait GradedBasisElement {

}

pub trait GradedLinearTransformation {

}

pub trait GradedVectorSpace {
    fn get_basis(&self) -> &Vec<impl GradedBasisElement>;
}