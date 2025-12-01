#[cfg(test)]
mod tests {
    use ahash::{HashMap};
    use crate::{
        basiselement::BasisElement,
        grading::{BiGrading, UniGrading},
        linalg::{
            field::F2,
            flat_matrix::FlatMatrix,
            matrix::RModMorphism,
            ring::{CRing, UniPolRing}
        },
        module::{
            module::GradedModule,
            morphism::GradedModuleMap
        }
    };

    #[derive(Debug, Clone, PartialEq, Default)]
    struct TestBasis;
    
    impl BasisElement for TestBasis {}



    #[test]
    fn test_k_t_to_k_t_morphism_cokernel() {
        use crate::linalg::ring::UniPolRing;
        
        // Create k[t] modules - domain and codomain are both k[t]
        let mut codomain_elements = HashMap::default(); 
        codomain_elements.insert(BiGrading(0, 0), vec![(TestBasis, UniGrading(0), None)]); // k[t] generator (no quotient)
        let codomain = GradedModule(codomain_elements);
        
        // Create morphism that maps 1 -> t^2
        let mut maps = HashMap::default();
        let mut matrix = FlatMatrix::zero(1, 1);
        matrix.set(0, 0, UniPolRing(F2::one(), 2)); // coefficient 1, degree 2 (represents t^2)
        maps.insert(BiGrading(0, 0), matrix);
        
        let mut explained = HashMap::default();
        explained.insert(BiGrading(0, 0), vec![((BiGrading(0, 0), 0), 2)]); // maps to grade (0,0), element 0, with t^2
        
        let morphism = GradedModuleMap {
            maps
        };
        
        // Compute cokernel
        let (coker_map, _coker_inv, cokernel) = morphism.cokernel(&codomain);
        
        // The cokernel should be k[t] / (t^2) which has one generator with quotient t^2
        assert_eq!(cokernel.0.len(), 1, "Cokernel should have one grade");
        
        let coker_elements = cokernel.0.get(&BiGrading(0, 0)).expect("Should have elements in grade (0,0)");
        assert_eq!(coker_elements.len(), 1, "Cokernel should have exactly one element");
        
        let (_, _, quotient) = &coker_elements[0];
        assert_eq!(*quotient, Some(2), "Element should be quotient by t^2");
        
        // Verify the cokernel map dimensions
        let coker_map_matrix = coker_map.maps.get(&BiGrading(0, 0)).expect("Cokernel map should exist in grade (0,0)");
        assert_eq!(coker_map_matrix.domain, 1, "Cokernel map should have domain dimension 1");
        assert_eq!(coker_map_matrix.codomain, 1, "Cokernel map should have codomain dimension 1");
    }

    #[test]
    fn test_tau_grading_cokernel_1() {
        let mut codomain_elements = HashMap::default(); 
        codomain_elements.insert(UniGrading(0), vec![(TestBasis, UniGrading(5), None)]); 
        let codomain = GradedModule(codomain_elements);

        let mut maps = HashMap::default();
        let mut matrix = FlatMatrix::zero(1, 1);
        matrix.set(0, 0, UniPolRing(F2::one(), 0)); // coefficient 1, degree 2 (represents t^2)
        maps.insert(UniGrading(0), matrix);

               
        let morphism = GradedModuleMap { maps };
        
        let (_, _, coker) = morphism.cokernel(&codomain);
        assert_eq!(coker.dimension_in_grade(&UniGrading(0)), 0); 
    }

    #[test]
    fn test_tau_grading_cokernel_2() {
        use crate::linalg::ring::UniPolRing;
        
        // Create k[t] modules - domain and codomain are both k[t]

        let mut codomain_elements = HashMap::default(); 
        codomain_elements.insert(UniGrading(0), vec![(TestBasis, UniGrading(0), None)]); // k[t] generator (no quotient)
        let codomain = GradedModule(codomain_elements);
        
        
        let mut maps = HashMap::default();
        let mut matrix = FlatMatrix::zero(1, 1);
        matrix.set(0, 0, UniPolRing(F2::one(), 2)); // coefficient 1, degree 2 (represents t^2)
        maps.insert(UniGrading(0), matrix);
        
        let morphism = GradedModuleMap { maps };
        
        // Compute cokernel
        let (coker_map, _coker_inv, cokernel) = morphism.cokernel(&codomain);
        
        // The cokernel should be k[t] / (t^2) which has one generator with quotient t^2
        assert_eq!(cokernel.0.len(), 1, "Cokernel should have one grade");
        
        let coker_elements_2 = cokernel.0.get(&UniGrading(0)).expect("Should have elements in grade (2,0)");
        assert_eq!(coker_elements_2.len(), 1, "Cokernel should have exactly one element");
        
        let (_, _, quotient) = &coker_elements_2[0];
        assert_eq!(*quotient, Some(2), "Element should be quotient by t^2");
        
        // Verify the cokernel map dimensions
        let coker_map_matrix = coker_map.maps.get(&UniGrading(0)).expect("Cokernel map should exist in grade (2,0)");
        assert_eq!(coker_map_matrix.domain, 1, "Cokernel map should have domain dimension 1 in (2,0)");
        assert_eq!(coker_map_matrix.codomain, 1, "Cokernel map should have codomain dimension 1 in (2,0)");
    }


    #[test]
    fn test_tau_grading_cokernel_3() {
        let mut map = FlatMatrix::zero(1, 3);
        map.set(0, 0, UniPolRing::<F2>::zero());
        map.set(0, 1, UniPolRing::one());
        map.set(0, 2, UniPolRing::one());

        // 1 -> (0,1,1)

        let mut graded_map = HashMap::default();
        graded_map.insert(UniGrading(7), map);

        let graded = GradedModuleMap {
            maps: graded_map,
        };

        let mut hashmap: HashMap<UniGrading, Vec<(TestBasis, UniGrading, Option<u16>)>> = HashMap::default();
        hashmap.insert(UniGrading(7), vec![(TestBasis, UniGrading(1), None), (TestBasis, UniGrading(3), Some(1)), (TestBasis, UniGrading(3), None)]);
        let codomain = GradedModule(hashmap);
        
        let mut hashmap: HashMap<UniGrading, Vec<(TestBasis, UniGrading, Option<u16>)>> = HashMap::default();
        hashmap.insert(UniGrading(7), vec![(TestBasis, UniGrading(3), None)]);
        let domain = GradedModule(hashmap);
        
        graded.verify(&domain, &codomain).unwrap();

        let (_, _, coker) = graded.cokernel(&codomain);

        let el_1 = &coker.0.get(&UniGrading(7)).expect("Expected coker to contain an element.")[0];
        let el_2 = &coker.0.get(&UniGrading(7)).expect("Expected coker to contain an element.")[1];
        assert_eq!(el_1.1, UniGrading(3));
        assert_eq!(el_1.2, Some(1));
        assert_eq!(el_2.1, UniGrading(1));
        assert_eq!(el_2.2, None);
    }

    #[test]
    fn test_correct_coker_module_1() {
        let mut map = FlatMatrix::zero(1, 2);
        map.set(0, 0, UniPolRing::<F2>::one());
        map.set(0, 1, UniPolRing(F2::one(), 1));

        // 1 -> (1,t)

        let mut graded_map = HashMap::default();
        graded_map.insert(UniGrading(5), map);

        let graded = GradedModuleMap {
            maps: graded_map,
        };

        let mut hashmap: HashMap<UniGrading, Vec<(TestBasis, UniGrading, Option<u16>)>> = HashMap::default();
        hashmap.insert(UniGrading(5), vec![(TestBasis, UniGrading(0), None), (TestBasis, UniGrading(1), None)]);
        let codomain = GradedModule(hashmap);
        
        let mut hashmap: HashMap<UniGrading, Vec<(TestBasis, UniGrading, Option<u16>)>> = HashMap::default();
        hashmap.insert(UniGrading(5), vec![(TestBasis, UniGrading(0), None)]);
        let domain = GradedModule(hashmap);
        
        graded.verify(&domain, &codomain).unwrap();

        let (_, _, coker) = graded.cokernel(&codomain);

        let el = &coker.0.get(&UniGrading(5)).expect("Expected coker to contain an element.")[0];
        assert_eq!(el.1, UniGrading(1));
        assert_eq!(el.2, None);
    }  
    
    #[test]
    fn test_correct_coker_module_2() {
        let mut map = FlatMatrix::zero(1, 1);
        map.set(0, 0, UniPolRing::<F2>(F2::one(), 2));
        
        // k[t]/t^7 -> k[t]/t^4
        // 1 -> t^2

        let mut graded_map = HashMap::default();
        graded_map.insert(UniGrading(0), map);

        let graded = GradedModuleMap {
            maps: graded_map,
        };

        let mut hashmap: HashMap<UniGrading, Vec<(TestBasis, UniGrading, Option<u16>)>> = HashMap::default();
        hashmap.insert(UniGrading(0), vec![(TestBasis, UniGrading(1), Some(7))]);
        let domain = GradedModule(hashmap);
        
        let mut hashmap: HashMap<UniGrading, Vec<(TestBasis, UniGrading, Option<u16>)>> = HashMap::default();
        hashmap.insert(UniGrading(0), vec![(TestBasis, UniGrading(3), Some(4))]);
        let codomain = GradedModule(hashmap);
        
        graded.verify(&domain, &codomain).unwrap();

        let (_, _, coker) = graded.cokernel(&codomain);

        let el = &coker.0.get(&UniGrading(0)).expect("Expected coker to contain an element.")[0];
        assert_eq!(el.1, UniGrading(3));
        assert_eq!(el.2, Some(2));
    }     

    #[test]
    fn test_kernel() {
        let mut map = FlatMatrix::zero(1, 1);
        map.set(0, 0, UniPolRing::<F2>::one());
        
        // k[t]/t -> k[t]/t
        // 1 -> 1

        let mut graded_map = HashMap::default();
        graded_map.insert(UniGrading(0), map);

        let graded = GradedModuleMap {
            maps: graded_map,
        };

        let mut hashmap: HashMap<UniGrading, Vec<(TestBasis, UniGrading, Option<u16>)>> = HashMap::default();
        hashmap.insert(UniGrading(0), vec![(TestBasis, UniGrading(0), Some(1))]);
        let domain = GradedModule(hashmap);
        
        let mut hashmap: HashMap<UniGrading, Vec<(TestBasis, UniGrading, Option<u16>)>> = HashMap::default();
        hashmap.insert(UniGrading(0), vec![(TestBasis, UniGrading(0), Some(1))]);
        let codomain = GradedModule(hashmap);
        
        graded.verify(&domain, &codomain).unwrap();

        let pivot = graded.basis_element_kernel_pivot_in_grade(&domain, &codomain, UniGrading(0));

        assert!(pivot.is_none())
    }      

    #[test]
    fn test_coker_module_structure() {
        let z = UniPolRing::<F2>::zero();
        let o = UniPolRing::<F2>::one();
        let t = UniPolRing(F2::one(), 1);

        let mut mat = FlatMatrix::zero(1, 2);
        mat.set(0, 0, t);
        mat.set(0, 1, o);


        let codom = vec![
            (TestBasis, UniGrading(1), Some(3)),
            (TestBasis, UniGrading(0), Some(1)),
        ];


        let mut a = HashMap::default();
        a.insert(UniGrading(0), codom);
        let codomain = GradedModule(a);
        
        let mut a = HashMap::default();
        a.insert(UniGrading(0), mat);        
        let maps = GradedModuleMap {
            maps: a,
        };

        let (_,_,coker) = maps.cokernel(&codomain);
        

        assert_eq!(coker.0.get(&UniGrading(0)).unwrap()[0].2, Some(2));
    }
}