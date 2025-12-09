#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use ahash::HashMap;
    use algebra::{field::Field, matrices::flat_matrix::FlatMatrix, matrix::Matrix, rings::{finite_fields::{F2, Fp}, univariate_polynomial_ring::UniPolRing}};

    use crate::{
        basiselement::kBasisElement, comodule::{
            kcoalgebra::kCoalgebra, rcoalgebra::{A0_C, A1_C, tensor_k_coalgebra}, rcomodule::{RCoalgebra, RComodule}, traits::Comodule
        }, graded_module::GradedModule, graded_module_morphism::GradedModuleMap, grading::{BiGrading, Grading, UniGrading}, tensor::TensorMap
    };

    /// Helper function to create a simple test RCoalgebra
    fn create_test_coalgebra<G: Grading, F: Field>() -> RCoalgebra<G, F> {
        let zero = G::zero();
        
        // Create basis element
        let el = kBasisElement {
            name: "1".to_string(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let space_map = [(zero, vec![(el, UniGrading(0), None)])].into_iter().collect();
        let space = GradedModule(space_map);

        // Create coaction map
        let coact_map: HashMap<G, FlatMatrix<UniPolRing<F>>> =
            [(zero, FlatMatrix::identity(1))].into_iter().collect();
        let mut domain_explain = HashMap::default();
        domain_explain.insert(zero, vec![((zero, 0), 0)]);
        let coaction = GradedModuleMap {
            maps: coact_map,
        };

        // Create tensor structure
        let mut dimensions = HashMap::default();
        dimensions.insert(zero, 1);
        let mut construct = HashMap::default();
        let mut first_entry = HashMap::default();
        first_entry.insert((zero, 0), (zero, 0));
        construct.insert((zero, 0), first_entry);
        let mut deconstruct = HashMap::default();
        deconstruct.insert((zero, 0), ((zero, 0), (zero, 0)));

        let tensor = TensorMap {
            construct,
            deconstruct,
            dimensions,
        };

        RCoalgebra {
            space,
            coaction,
            tensor,
        }
    }

    type BiG = crate::grading::BiGrading;

    // Test for RComodule::zero_comodule
    #[test]
    fn test_zero_comodule() {
        let coalgebra = Arc::new(create_test_coalgebra::<BiG, F2>());
        let comodule = RComodule::zero_comodule(coalgebra.clone());

        assert_eq!(
            comodule.coalgebra.as_ref() as *const _,
            coalgebra.as_ref() as *const _
        );
        assert!(comodule.space.0.is_empty());
        assert!(comodule.coaction.maps.is_empty());
        assert!(comodule.tensor.dimensions.is_empty());
    }

    // Test for RComodule::fp_comodule
    #[test]
    fn test_fp_comodule() {
        let coalgebra = Arc::new(create_test_coalgebra::<BiG, F2>());
        let comodule = RComodule::fp_comodule(coalgebra.clone(), BiGrading::zero());

        let zero = BiG::zero();
        
        // Check space structure
        assert_eq!(comodule.space.0.len(), 1);
        assert!(comodule.space.0.contains_key(&zero));
        
        let elements = comodule.space.0.get(&zero).unwrap();
        assert_eq!(elements.len(), 1);
        assert_eq!(elements[0].0.name, "fp");
        assert_eq!(elements[0].0.generator, false);
        assert_eq!(elements[0].0.generated_index, 0);

        // Check coaction structure
        assert!(comodule.coaction.maps.contains_key(&zero));
        assert_eq!(comodule.coaction.maps.len(), 1);
        assert_eq!(comodule.coaction.maps[&zero], FlatMatrix::identity(1));

        // Check tensor structure
        assert_eq!(comodule.tensor.dimensions.get(&zero), Some(&1));
        assert!(comodule.tensor.construct.contains_key(&(zero, 0)));
        assert!(comodule.tensor.deconstruct.contains_key(&(zero, 0)));
    }

    // Test for RComodule::get_generators
    #[test]
    fn test_get_generators_empty() {
        let coalgebra = Arc::new(create_test_coalgebra::<UniGrading, F2>());
        let comodule = RComodule::zero_comodule(coalgebra);
        
        let generators = comodule.get_generators();
        assert!(generators.is_empty());
    }

    #[test]
    fn test_get_generators_with_generators() {
        let coalgebra = Arc::new(create_test_coalgebra::<UniGrading, F2>());
        let zero = UniGrading(0);
        let one = UniGrading(1);
        
        let mut space_map = HashMap::default();
        space_map.insert(
            zero,
            vec![
                (kBasisElement {
                    name: "gen1".to_string(),
                    generator: true,
                    primitive: None,
                    generated_index: 0,
                }, UniGrading(0), None),
                (kBasisElement {
                    name: "non_gen".to_string(),
                    generator: false,
                    primitive: None,
                    generated_index: 1,
                }, UniGrading(0), None),
            ],
        );
        space_map.insert(
            one,
            vec![
                (kBasisElement {
                    name: "gen2".to_string(),
                    generator: true,
                    primitive: None,
                    generated_index: 2,
                }, UniGrading(0), None),
            ],
        );

        let comodule = RComodule {
            coalgebra,
            space: GradedModule(space_map),
            coaction: GradedModuleMap::default(),
            tensor: TensorMap::default(),
        };

        let generators = comodule.get_generators();
        
        // Should have 2 generators
        assert_eq!(generators.len(), 2);
        
        // Check that we have generators with correct indices and grades
        let gen_indices: Vec<_> = generators.iter().map(|(idx, _g, _name)| *idx).collect();
        assert!(gen_indices.contains(&0));
        assert!(gen_indices.contains(&2));
        
        let gen_grades: Vec<_> = generators.iter().map(|(_idx, g, _name)| *g).collect();
        assert!(gen_grades.contains(&zero));
        assert!(gen_grades.contains(&one));
    }

    // Test for RComodule::direct_sum
    #[test]
    fn test_direct_sum_with_zero() {
        let coalgebra = Arc::new(create_test_coalgebra::<BiG, F2>());
        let mut comodule1 = RComodule::fp_comodule(coalgebra.clone(), BiGrading::zero());
        let mut comodule2 = RComodule::zero_comodule(coalgebra);

        let original_dims = comodule1.space.0.len();
        
        comodule1.direct_sum(&mut comodule2);

        // Should still have the same structure as before
        assert_eq!(comodule1.space.0.len(), original_dims);
        assert!(comodule2.space.0.is_empty());
    }

    #[test]
    fn test_direct_sum_two_fp_comodules() {
        let coalgebra = Arc::new(create_test_coalgebra::<BiG, F2>());
        let mut comodule1 = RComodule::fp_comodule(coalgebra.clone(), BiGrading::zero());
        let mut comodule2 = RComodule::fp_comodule(coalgebra, BiGrading::zero());

        let zero = BiG::zero();
        
        comodule1.direct_sum(&mut comodule2);

        // Check that we now have 2 elements at grade 0
        assert_eq!(comodule1.space.0.get(&zero).unwrap().len(), 2);
        
        // Check that both elements have the correct name
        let elements = comodule1.space.0.get(&zero).unwrap();
        assert_eq!(elements[0].0.name, "fp");
        assert_eq!(elements[1].0.name, "fp");
        
        // Check tensor dimensions updated correctly
        assert_eq!(comodule1.tensor.dimensions.get(&zero), Some(&2));
        
        // The second comodule should have empty vectors after draining
        // (but the HashMap keys might still exist)
        if let Some(elements) = comodule2.space.0.get(&zero) {
            assert!(elements.is_empty());
        }
    }

    // Test for RComodule::cofree_comodule
    #[test]
    fn test_cofree_comodule_basic() {
        let coalgebra = Arc::new(create_test_coalgebra::<BiG, F2>());
        let zero = BiG::zero();
        let one = BiGrading(1, 0);
        
        let cofree = RComodule::cofree_comodule(coalgebra.clone(), 5, zero, one, (UniGrading(0), None));

        // Should have elements at grade 0
        assert!(cofree.space.0.contains_key(&zero));
        
        let elements = cofree.space.0.get(&zero).unwrap();
        assert!(!elements.is_empty());
        
        // All elements should have the specified generated_index
        for (element, _, _) in elements {
            assert_eq!(element.generated_index, 5);
        }

        // Should have coaction maps
        assert!(cofree.coaction.maps.contains_key(&zero));
    }

    #[test]
    fn test_cofree_shifted_is_correct() {
        let coalgebra = Arc::new(A0_C());
        let shift = RComodule::cofree_comodule(coalgebra, 0, UniGrading(2), UniGrading(50), (UniGrading(0), None));
        shift.verify().unwrap();
    }

    #[test]
    fn test_cofree_shifted_is_correct_difficult() {
        let input = include_str!("../../../examples/direct/A(1).txt");
        let coalgebra = kCoalgebra::parse(input, UniGrading::infty()).unwrap().0;
        let tensor_coalgebra = Arc::new(tensor_k_coalgebra(coalgebra));

        let shift = RComodule::cofree_comodule(tensor_coalgebra, 0, UniGrading(0), UniGrading(50), (UniGrading(0), None));
        shift.verify().unwrap();
    }

    #[test]
    fn test_cofree_comodule_with_limit() {
        let coalgebra = Arc::new(create_test_coalgebra::<BiG, F2>());
        let zero = BiG::zero();
        let two = BiGrading(2, 0);
        
        let cofree = RComodule::cofree_comodule(coalgebra.clone(), 10, zero, two, (UniGrading(0), None));

        // Check that we respect the grading limit
        for (grade, _) in &cofree.space.0 {
            assert!(grade <= &two);
        }
        
        // Should have tensor structure within the limit
        for (grade, _) in &cofree.tensor.dimensions {
            assert!(grade <= &two);
        }
    }

    // Test type compatibility with different fields
    #[test] 
    fn test_different_fields() {
        // Test with F2
        let coalgebra_f2 = Arc::new(create_test_coalgebra::<BiG, F2>());
        let _comodule_f2 = RComodule::fp_comodule(coalgebra_f2, BiGrading::zero());

        // Test with F3
        let coalgebra_f3 = Arc::new(create_test_coalgebra::<BiG, Fp<3>>());
        let _comodule_f3 = RComodule::fp_comodule(coalgebra_f3, BiGrading::zero());

        // Test with F5  
        let coalgebra_f5 = Arc::new(create_test_coalgebra::<BiG, Fp<5>>());
        let _comodule_f5 = RComodule::fp_comodule(coalgebra_f5, BiGrading::zero());
    }

    // Test that RCoalgebra can be cloned
    #[test]
    fn test_coalgebra_clone() {
        let coalgebra1 = create_test_coalgebra::<BiG, F2>();
        let coalgebra2 = coalgebra1.clone();
        
        assert_eq!(coalgebra1, coalgebra2);
    }

    // Test clone and equality for RComodule
    #[test]
    fn test_comodule_clone_and_equality() {
        let coalgebra = Arc::new(create_test_coalgebra::<BiG, F2>());
        let comodule1 = RComodule::fp_comodule(coalgebra.clone(), BiGrading::zero());
        let comodule2 = comodule1.clone();

        assert_eq!(comodule1, comodule2);
        
        // Verify they share the same coalgebra reference
        assert_eq!(
            Arc::as_ptr(&comodule1.coalgebra),
            Arc::as_ptr(&comodule2.coalgebra)
        );
    }

    // Test edge cases
    #[test]
    fn test_cofree_comodule_zero_grade() {
        let coalgebra = Arc::new(create_test_coalgebra::<BiG, F2>());
        let zero = BiG::zero();
        
        let cofree = RComodule::cofree_comodule(coalgebra, 0, zero, zero, (UniGrading(0), None));
        
        // Should have elements only at zero grade
        assert_eq!(cofree.space.0.len(), 1);
        assert!(cofree.space.0.contains_key(&zero));
    }

    #[test] 
    fn test_multiple_direct_sums() {
        let coalgebra = Arc::new(create_test_coalgebra::<BiG, F2>());
        let mut base = RComodule::zero_comodule(coalgebra.clone());
        
        // Add multiple fp_comodules via direct sum
        for _ in 0..5 {
            let mut fp = RComodule::fp_comodule(coalgebra.clone(), BiGrading::zero());
            base.direct_sum(&mut fp);
        }
        
        let zero = BiG::zero();
        assert_eq!(base.space.0.get(&zero).unwrap().len(), 5);
        assert_eq!(base.tensor.dimensions.get(&zero), Some(&5));
    }


    #[test] 
    fn test_multiple_direct_sums_in_different_grade_1() {
        let coalgebra = Arc::new(A1_C());
        let mut base = RComodule::zero_comodule(coalgebra.clone());
        
        
        let mut s = RComodule::cofree_comodule(coalgebra.clone(), 0, UniGrading(0), UniGrading(70), (UniGrading(0), None));
        base.direct_sum(&mut s);
        
        
        let mut s = RComodule::cofree_comodule(coalgebra.clone(), 1, UniGrading(20), UniGrading(70), (UniGrading(0), None));
        base.direct_sum(&mut s);
    }

    #[test] 
    fn test_multiple_direct_sums_in_different_grade_2() {
        let coalgebra = Arc::new(A1_C());
        let mut base = RComodule::zero_comodule(coalgebra.clone());
        
        
        let mut s = RComodule::cofree_comodule(coalgebra.clone(), 0, UniGrading(0), UniGrading(70), (UniGrading(0), None));
        base.direct_sum(&mut s);
        
        let mut s = RComodule::cofree_comodule(coalgebra.clone(), 1, UniGrading(1), UniGrading(70), (UniGrading(0), None));
        base.direct_sum(&mut s);
        
        
    }
}