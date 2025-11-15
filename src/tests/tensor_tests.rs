use crate::basiselement::BasisElement;
use crate::linalg::graded::GradedVectorSpace;
use crate::grading::OrderedGrading;

#[derive(Debug, Clone, Default)]
struct MockBasisElement;

impl BasisElement for MockBasisElement {}

#[cfg(test)]
mod tests {
    use ahash::HashMap;

    use crate::{comodule::kcoalgebra::A0_coalgebra, grading::UniGrading, tensor::Tensor};

    use super::*;

    fn create_mock_vector_space<G: OrderedGrading, B: BasisElement>(
        elements: Vec<(G, Vec<B>)>,
    ) -> GradedVectorSpace<G, B> {
        GradedVectorSpace(elements.into_iter().collect())
    }

    #[test]
    fn test_new_ktensor() {
        let tensor: Tensor<UniGrading> = Tensor::default();
        assert!(tensor.construct.is_empty());
        assert!(tensor.deconstruct.is_empty());
        assert!(tensor.dimensions.is_empty());
    }

    #[test]
    fn test_generate_ktensor() {
        let left_elements = vec![(UniGrading(0), vec![MockBasisElement]), (UniGrading(1), vec![MockBasisElement])];
        let right_elements = vec![(UniGrading(0), vec![MockBasisElement]), (UniGrading(1), vec![MockBasisElement])];

        let left_space = create_mock_vector_space(left_elements);
        let right_space = create_mock_vector_space(right_elements);

        let tensor = Tensor::generate(&left_space, &right_space);

        assert_eq!(tensor.dimensions.get(&UniGrading(0)), Some(&1));
        assert_eq!(tensor.dimensions.get(&UniGrading(1)), Some(&2));

        let mut map = HashMap::default();
        map.insert((UniGrading(0), 0), ((UniGrading(0), 0), (UniGrading(0), 0)));
        map.insert((UniGrading(1), 0), ((UniGrading(0), 0), (UniGrading(1), 0)));
        map.insert((UniGrading(1), 1), ((UniGrading(1), 0), (UniGrading(0), 0)));

        assert_eq!(tensor.deconstruct, map);
    }

    #[test]
    fn test_add() {
        let left_elements = vec![(UniGrading(0), vec![MockBasisElement]), (UniGrading(1), vec![MockBasisElement])];
        let right_elements = vec![(UniGrading(0), vec![MockBasisElement]), (UniGrading(1), vec![MockBasisElement])];

        let left_space = create_mock_vector_space(left_elements);
        let right_space = create_mock_vector_space(right_elements);

        let tensor = Tensor::<UniGrading>::generate(&left_space, &right_space);
        assert_eq!(tensor.dimensions.get(&UniGrading(0)), Some(&1));
        assert_eq!(tensor.dimensions.get(&UniGrading(1)), Some(&2));

        let new_tensor = tensor.add_and_restrict(UniGrading(1), UniGrading(2));

        assert_eq!(new_tensor.dimensions.get(&UniGrading(1)), Some(&1));
        assert_eq!(new_tensor.dimensions.get(&UniGrading(2)), Some(&2));
        assert!(new_tensor.dimensions.get(&UniGrading(0)).is_none());
    }

    #[test]
    fn test_restrict() {
        let left_elements = vec![(UniGrading(0), vec![MockBasisElement]), (UniGrading(1), vec![MockBasisElement])];
        let right_elements = vec![(UniGrading(0), vec![MockBasisElement]), (UniGrading(1), vec![MockBasisElement])];

        let left_space = create_mock_vector_space(left_elements);
        let right_space = create_mock_vector_space(right_elements);

        let tensor = Tensor::<UniGrading>::generate(&left_space, &right_space);
        assert_eq!(tensor.dimensions.get(&UniGrading(0)), Some(&1));
        assert_eq!(tensor.dimensions.get(&UniGrading(1)), Some(&2));

        let new_tensor = tensor.add_and_restrict(UniGrading(1), UniGrading(1));

        assert_eq!(new_tensor.dimensions.get(&UniGrading(1)), Some(&1));
        assert_eq!(new_tensor.dimensions.get(&UniGrading(2)), None);
        assert!(new_tensor.dimensions.get(&UniGrading(0)).is_none());
    }

    #[test]
    fn test_direct_sum() {
        let left_elements = vec![(UniGrading(0), vec![MockBasisElement]), (UniGrading(1), vec![MockBasisElement])];
        let right_elements = vec![(UniGrading(0), vec![MockBasisElement]), (UniGrading(1), vec![MockBasisElement])];

        let left_space = create_mock_vector_space(left_elements.clone());
        let right_space = create_mock_vector_space(right_elements.clone());

        let mut tensor1 = Tensor::<UniGrading>::generate(&left_space, &right_space);
        let mut tensor2 = Tensor::<UniGrading>::generate(
            &create_mock_vector_space(left_elements),
            &create_mock_vector_space(right_elements),
        );

        let mut self_space_dimensions = HashMap::default();
        self_space_dimensions.insert(UniGrading(0), 1);
        self_space_dimensions.insert(UniGrading(1), 1);
        tensor1.direct_sum(&mut tensor2, &self_space_dimensions);

        assert_eq!(tensor1.dimensions.get(&UniGrading(0)), Some(&2));
        assert_eq!(tensor1.dimensions.get(&UniGrading(1)), Some(&4));
        assert_eq!(tensor1.dimensions.get(&UniGrading(2)), None);
        tensor1.is_correct();
    }

    #[test]
    fn test_correctness_a0() {
        let a0 = A0_coalgebra();
        a0.tensor.is_correct();
        assert_eq!(a0.tensor.dimensions[&UniGrading(1)], 2);
    }
}
