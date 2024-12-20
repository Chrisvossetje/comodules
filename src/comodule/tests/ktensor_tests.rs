use crate::linalg::graded::{BasisElement, GradedVectorSpace};
use crate::linalg::grading::Grading;

use std::collections::HashMap;

#[derive(Debug, Clone)]
struct MockBasisElement;

impl BasisElement for MockBasisElement {}

#[cfg(test)]
mod tests {
    use crate::comodule::{kcoalgebra::A0_coalgebra, ktensor::kTensor};

    use super::*;

    fn create_mock_vector_space<G: Grading, B: BasisElement>(
        elements: Vec<(G, Vec<B>)>,
    ) -> GradedVectorSpace<G, B> {
        GradedVectorSpace(elements.into_iter().collect())
    }

    #[test]
    fn test_new_ktensor() {
        let tensor: kTensor<i32> = kTensor::new();
        assert!(tensor.construct.is_empty());
        assert!(tensor.deconstruct.is_empty());
        assert!(tensor.dimensions.is_empty());
    }

    #[test]
    fn test_generate_ktensor() {
        let left_elements = vec![(0, vec![MockBasisElement]), (1, vec![MockBasisElement])];
        let right_elements = vec![(0, vec![MockBasisElement]), (1, vec![MockBasisElement])];

        let left_space = create_mock_vector_space(left_elements);
        let right_space = create_mock_vector_space(right_elements);

        let tensor = kTensor::generate(&left_space, &right_space);

        assert_eq!(tensor.dimensions.get(&0), Some(&1));
        assert_eq!(tensor.dimensions.get(&1), Some(&2));

        let mut map = HashMap::default();
        map.insert((0, 0), ((0, 0), (0, 0)));
        map.insert((1, 0), ((0, 0), (1, 0)));
        map.insert((1, 1), ((1, 0), (0, 0)));

        assert_eq!(tensor.deconstruct, map);
    }

    #[test]
    fn test_add() {
        let left_elements = vec![(0, vec![MockBasisElement]), (1, vec![MockBasisElement])];
        let right_elements = vec![(0, vec![MockBasisElement]), (1, vec![MockBasisElement])];

        let left_space = create_mock_vector_space(left_elements);
        let right_space = create_mock_vector_space(right_elements);

        let tensor = kTensor::<i32>::generate(&left_space, &right_space);
        assert_eq!(tensor.dimensions.get(&0), Some(&1));
        assert_eq!(tensor.dimensions.get(&1), Some(&2));

        let new_tensor = tensor.add_and_restrict(1, 2);

        assert_eq!(new_tensor.dimensions.get(&1), Some(&1));
        assert_eq!(new_tensor.dimensions.get(&2), Some(&2));
        assert!(new_tensor.dimensions.get(&0).is_none());
    }

    #[test]
    fn test_restrict() {
        let left_elements = vec![(0, vec![MockBasisElement]), (1, vec![MockBasisElement])];
        let right_elements = vec![(0, vec![MockBasisElement]), (1, vec![MockBasisElement])];

        let left_space = create_mock_vector_space(left_elements);
        let right_space = create_mock_vector_space(right_elements);

        let tensor = kTensor::<i32>::generate(&left_space, &right_space);
        assert_eq!(tensor.dimensions.get(&0), Some(&1));
        assert_eq!(tensor.dimensions.get(&1), Some(&2));

        let new_tensor = tensor.add_and_restrict(1, 1);

        assert_eq!(new_tensor.dimensions.get(&1), Some(&1));
        assert_eq!(new_tensor.dimensions.get(&2), None);
        assert!(new_tensor.dimensions.get(&0).is_none());
    }

    #[test]
    fn test_direct_sum() {
        let left_elements = vec![(0, vec![MockBasisElement]), (1, vec![MockBasisElement])];
        let right_elements = vec![(0, vec![MockBasisElement]), (1, vec![MockBasisElement])];

        let left_space = create_mock_vector_space(left_elements.clone());
        let right_space = create_mock_vector_space(right_elements.clone());

        let mut tensor1 = kTensor::<i32>::generate(&left_space, &right_space);
        let mut tensor2 = kTensor::<i32>::generate(
            &create_mock_vector_space(left_elements),
            &create_mock_vector_space(right_elements),
        );

        let self_space_dimensions = HashMap::from([(0, 1), (1, 1)]);
        tensor1.direct_sum(&mut tensor2, &self_space_dimensions);

        assert_eq!(tensor1.dimensions.get(&0), Some(&2));
        assert_eq!(tensor1.dimensions.get(&1), Some(&4));
        assert_eq!(tensor1.dimensions.get(&2), None);
        tensor1.is_correct();
    }

    #[test]
    fn test_correctness_a0() {
        let a0 = A0_coalgebra();
        a0.tensor.is_correct();
        assert_eq!(a0.tensor.dimensions[&1], 2);
    }
}
