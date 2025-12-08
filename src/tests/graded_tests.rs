#[cfg(test)]
mod tests {

    use ahash::HashMap;

    use crate::linalg::field::F2;
    use crate::basiselement::BasisElement;
    use crate::linalg::linalg::FlatMatrix;
    use crate::linalg::graded::{GradedLinearMap, GradedVectorSpace};
    use crate::grading::UniGrading;
    use crate::linalg::matrix::Matrix;
    
    // Assuming F2 is implemented as a Field type.

    type G = UniGrading; // OrderedGrading type
    type F = F2; // Field type
    type M = FlatMatrix<F>; // Matrix type

    impl BasisElement for usize {}

    #[test]
    fn test_graded_vector_space_new() {
        // Test creation of a new empty GradedVectorSpace
        let graded_space: GradedVectorSpace<G, usize> = GradedVectorSpace::new();
        assert!(graded_space.0.is_empty());
    }

    #[test]
    fn test_graded_vector_space_from() {
        // Test creating a GradedVectorSpace from a HashMap
        let mut space = HashMap::default();
        space.insert(UniGrading(0), vec![1, 2, 3]);
        space.insert(UniGrading(1), vec![4, 5]);

        let graded_space: GradedVectorSpace<G, usize> = GradedVectorSpace::from(space.clone());
        assert_eq!(graded_space.0, space);
    }

    #[test]
    fn test_graded_linear_map_from() {
        // Test creating a GradedLinearMap from a HashMap
        let mut map = HashMap::default();
        map.insert(UniGrading(0), FlatMatrix::zero(2, 2));
        map.insert(UniGrading(1), FlatMatrix::zero(3, 3));

        let linear_map: GradedLinearMap<G, F, M> = GradedLinearMap::from(map.clone());
        assert_eq!(linear_map.maps, map);
    }

    #[test]
    fn test_graded_linear_map_get_cokernel() {
        // Test get_cokernel() produces a correct GradedLinearMap
        let mut map = HashMap::default();
        map.insert(UniGrading(0), FlatMatrix::zero(2, 3));
        map.insert(UniGrading(1), FlatMatrix::zero(3, 4));

        let linear_map: GradedLinearMap<G, F, M> = GradedLinearMap::from(map);
        let cokernel = linear_map.get_cokernel();

        // The dimensions of the cokernel matrices depend on the implementation of FlatMatrix::cokernel()
        assert!(cokernel.maps.get(&UniGrading(0)).is_some());
        assert!(cokernel.maps.get(&UniGrading(1)).is_some());
    }

    #[test]
    fn test_graded_linear_map_get_kernel() {
        // Test get_kernel() produces a correct GradedLinearMap
        let mut map = HashMap::default();
        map.insert(UniGrading(0), FlatMatrix::zero(3, 2));
        map.insert(UniGrading(1), FlatMatrix::zero(4, 3));

        let linear_map: GradedLinearMap<G, F, M> = GradedLinearMap::from(map);
        let kernel = linear_map.get_kernel();

        // The dimensions of the kernel matrices depend on the implementation of FlatMatrix::kernel()
        assert!(kernel.maps.get(&UniGrading(0)).is_some());
        assert!(kernel.maps.get(&UniGrading(1)).is_some());
    }

    #[test]
    fn test_graded_linear_map_vstack() {
        // Test vstack combines matrices vertically for each grade
        let mut map1 = HashMap::default();
        map1.insert(UniGrading(0), FlatMatrix::zero(2, 3));
        map1.insert(UniGrading(1), FlatMatrix::zero(2, 4));

        let mut map2 = HashMap::default();
        map2.insert(UniGrading(0), FlatMatrix::zero(2, 6));
        map2.insert(UniGrading(1), FlatMatrix::zero(2, 2));

        let mut linear_map1: GradedLinearMap<G, F, M> = GradedLinearMap::from(map1);
        let mut linear_map2: GradedLinearMap<G, F, M> = GradedLinearMap::from(map2);

        linear_map1.vstack(&mut linear_map2);

        assert_eq!(linear_map1.maps[&UniGrading(0)].domain, 2); // 2 + 3 rows
        assert_eq!(linear_map1.maps[&UniGrading(1)].domain, 2); // 2 + 3 rows
        assert_eq!(linear_map1.maps[&UniGrading(0)].codomain, 9); // 2 + 3 rows
        assert_eq!(linear_map1.maps[&UniGrading(1)].codomain, 6); // 1 + 2 rows
    }

    #[test]
    fn test_graded_linear_map_block_sum() {
        // Test block_sum combines matrices block-wise
        let mut map1 = HashMap::default();
        map1.insert(UniGrading(0), FlatMatrix::zero(2, 3));
        map1.insert(UniGrading(1), FlatMatrix::zero(1, 2));

        let mut map2 = HashMap::default();
        map2.insert(UniGrading(0), FlatMatrix::zero(2, 4));
        map2.insert(UniGrading(1), FlatMatrix::zero(1, 3));

        let mut linear_map1: GradedLinearMap<G, F, M> = GradedLinearMap::from(map1);
        let mut linear_map2: GradedLinearMap<G, F, M> = GradedLinearMap::from(map2);

        linear_map1.block_sum(&mut linear_map2);

        assert_eq!(linear_map1.maps[&UniGrading(0)].domain, 4);
        assert_eq!(linear_map1.maps[&UniGrading(0)].codomain, 7);
        assert_eq!(linear_map1.maps[&UniGrading(1)].domain, 2);
        assert_eq!(linear_map1.maps[&UniGrading(1)].codomain, 5);
    }

    #[test]
    fn test_graded_linear_map_compose() {
        // Test composing two graded linear maps
        let mut map1 = HashMap::default();
        map1.insert(UniGrading(0), FlatMatrix::zero(2, 3));
        map1.insert(UniGrading(1), FlatMatrix::zero(3, 2));

        let mut map2 = HashMap::default();
        map2.insert(UniGrading(0), FlatMatrix::zero(3, 4));
        map2.insert(UniGrading(1), FlatMatrix::zero(2, 5));

        let linear_map1: GradedLinearMap<G, F, M> = GradedLinearMap::from(map1);
        let linear_map2: GradedLinearMap<G, F, M> = GradedLinearMap::from(map2);

        let composed = linear_map2.compose(&linear_map1);

        assert!(composed.maps.get(&UniGrading(0)).is_some());
        assert!(composed.maps.get(&UniGrading(1)).is_some());
        assert_eq!(composed.maps[&UniGrading(0)].domain, 2);
        assert_eq!(composed.maps[&UniGrading(0)].codomain, 4);
        assert_eq!(composed.maps[&UniGrading(1)].domain, 3);
        assert_eq!(composed.maps[&UniGrading(1)].codomain, 5);
    }

    #[test]
    fn test_graded_linear_map_pivots() {
        // Test pivots retrieves pivot positions
        let mut map = HashMap::default();
        map.insert(UniGrading(0), FlatMatrix::identity(3));
        map.insert(UniGrading(1), FlatMatrix::zero(4, 4));

        let linear_map: GradedLinearMap<G, F, M> = GradedLinearMap::from(map);
        let pivots = linear_map.pivots();

        assert!(pivots.get(&UniGrading(0)).is_some());
        assert!(pivots.get(&UniGrading(1)).is_some());
        assert_eq!(pivots[&UniGrading(0)], vec![(0, 0), (1, 1), (2, 2)]);
        assert_eq!(pivots[&UniGrading(1)], vec![]);
    }

    #[test]
    fn test_graded_linear_map_zero_codomain() {
        // Test zero_codomain produces a map with empty domain
        let mut codomain_space = HashMap::default();
        codomain_space.insert(UniGrading(0), vec![0; 3]);
        codomain_space.insert(UniGrading(1), vec![0; 2]);

        let codomain = GradedVectorSpace::from(codomain_space);
        let zero_map: GradedLinearMap<G, F, M> = GradedLinearMap::zero_codomain(&codomain);

        assert_eq!(zero_map.maps[&UniGrading(0)].codomain, 0);
        assert_eq!(zero_map.maps[&UniGrading(0)].domain, 3);
        assert_eq!(zero_map.maps[&UniGrading(1)].codomain, 0);
        assert_eq!(zero_map.maps[&UniGrading(1)].domain, 2);
    }
}
