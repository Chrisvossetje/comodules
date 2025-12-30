#[cfg(test)]
mod tests {
    use ahash::HashMap;
    use deepsize::DeepSizeOf;

    #[derive(Debug, Clone, Default, DeepSizeOf)]
    struct MockBasisElement;

    use crate::{
        grading::{grading::Grading, tensor::TensorMap, unigrading::UniGrading},
        k_comodule::{graded_space::GradedVectorSpace, kcoalgebra::A0_coalgebra},
    };

    fn create_mock_vector_space<G: Grading, B>(
        elements: Vec<(G, Vec<B>)>,
    ) -> GradedVectorSpace<G, B> {
        GradedVectorSpace(elements.into_iter().collect())
    }

    #[test]
    fn test_new_ktensor() {
        let tensor: TensorMap<UniGrading> = TensorMap::default();
        assert!(tensor.construct.is_empty());
        assert!(tensor.deconstruct.is_empty());
        assert!(tensor.dimensions.is_empty());
    }

    #[test]
    fn test_generate_ktensor() {
        let left_elements = vec![
            (UniGrading(0), vec![MockBasisElement]),
            (UniGrading(1), vec![MockBasisElement]),
        ];
        let right_elements = vec![
            (UniGrading(0), vec![MockBasisElement]),
            (UniGrading(1), vec![MockBasisElement]),
        ];

        let left_space = create_mock_vector_space(left_elements);
        let right_space = create_mock_vector_space(right_elements);

        let tensor = TensorMap::generate(&left_space, &right_space);

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
        let left_elements = vec![
            (UniGrading(0), vec![MockBasisElement]),
            (UniGrading(1), vec![MockBasisElement]),
        ];
        let right_elements = vec![
            (UniGrading(0), vec![MockBasisElement]),
            (UniGrading(1), vec![MockBasisElement]),
        ];

        let left_space = create_mock_vector_space(left_elements);
        let right_space = create_mock_vector_space(right_elements);

        let tensor = TensorMap::<UniGrading>::generate(&left_space, &right_space);
        assert_eq!(tensor.dimensions.get(&UniGrading(0)), Some(&1));
        assert_eq!(tensor.dimensions.get(&UniGrading(1)), Some(&2));

        let new_tensor = tensor.add_and_restrict(UniGrading(1), UniGrading(2));

        assert_eq!(new_tensor.dimensions.get(&UniGrading(1)), Some(&1));
        assert_eq!(new_tensor.dimensions.get(&UniGrading(2)), Some(&2));
        assert!(new_tensor.dimensions.get(&UniGrading(0)).is_none());
    }

    #[test]
    fn test_restrict() {
        let left_elements = vec![
            (UniGrading(0), vec![MockBasisElement]),
            (UniGrading(1), vec![MockBasisElement]),
        ];
        let right_elements = vec![
            (UniGrading(0), vec![MockBasisElement]),
            (UniGrading(1), vec![MockBasisElement]),
        ];

        let left_space = create_mock_vector_space(left_elements);
        let right_space = create_mock_vector_space(right_elements);

        let tensor = TensorMap::<UniGrading>::generate(&left_space, &right_space);
        assert_eq!(tensor.dimensions.get(&UniGrading(0)), Some(&1));
        assert_eq!(tensor.dimensions.get(&UniGrading(1)), Some(&2));

        let new_tensor = tensor.add_and_restrict(UniGrading(1), UniGrading(1));

        assert_eq!(new_tensor.dimensions.get(&UniGrading(1)), Some(&1));
        assert_eq!(new_tensor.dimensions.get(&UniGrading(2)), None);
        assert!(new_tensor.dimensions.get(&UniGrading(0)).is_none());
    }
}
