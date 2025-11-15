#[cfg(test)]
mod tests {
    use std::{sync::Arc};

    use ahash::HashMap;

    use crate::{
        comodule::{
            kcoalgebra::A0_coalgebra,
            kcomodule::kComodule,
            kmorphism::kComoduleMorphism,
            traits::{Comodule, ComoduleMorphism},
        },
        grading::{Grading, UniGrading},
        linalg::{
            field::F2, flat_matrix::FlatMatrix, graded::{GradedLinearMap, GradedVectorSpace}, matrix::RModMorphism, ring::CRing
        },
    };

    #[test]
    fn test_inject_codomain_to_cofree() {
        let coalgebra = Arc::new(A0_coalgebra());
        let comodule = Arc::new(kComodule::fp_comodule(coalgebra.clone(), UniGrading::zero()));

        let morphism = kComoduleMorphism::zero_morphism(comodule);

        let cofree_morphism = morphism.inject_codomain_to_cofree(UniGrading(5));

        let comp = kComodule::cofree_comodule(coalgebra, 0, UniGrading(0), UniGrading(5), ());

        // Assertions
        assert_eq!(cofree_morphism.codomain.space, comp.space);
        assert_eq!(cofree_morphism.codomain.tensor, comp.tensor);
        assert_eq!(cofree_morphism.codomain.coaction, comp.coaction);

        // This should not fail if the above succeeds
        assert_eq!(cofree_morphism.codomain.as_ref(), &comp);
    }

    #[test]
    fn test_cokernel() {
        let coalgebra = Arc::new(A0_coalgebra());

        let domain = Arc::new(kComodule::fp_comodule(coalgebra.clone(), UniGrading::zero()));

        let codomain = Arc::new(kComodule::cofree_comodule(coalgebra, 0, UniGrading(0), UniGrading(5), ()));

        let mut map: GradedLinearMap<UniGrading, F2, FlatMatrix<F2>> =
            GradedLinearMap::zero(&domain.space, &codomain.space);
        map.maps.get_mut(&UniGrading(0)).unwrap().data[0] = F2::one();

        let morphism = kComoduleMorphism {
            domain: domain.clone(),
            codomain: codomain.clone(),
            map,
        };

        let cokernel_morphism = morphism.cokernel();

        // Assertions
        assert_eq!(cokernel_morphism.domain, morphism.codomain);

        let mut map = HashMap::default();
        map.insert(UniGrading(0), FlatMatrix::zero(1, 0));
        map.insert(UniGrading(1), FlatMatrix::identity(1));
        let expected_map: GradedLinearMap<UniGrading, F2, FlatMatrix<F2>> = GradedLinearMap::from(map);
        assert_eq!(cokernel_morphism.map, expected_map);
    }

    #[test]
    fn test_structure_lines() {
        let coalgebra = Arc::new(A0_coalgebra());

        let domain = Arc::new(kComodule::cofree_comodule(coalgebra.clone(), 0, UniGrading(0), UniGrading(5), ()));

        let codomain = Arc::new(kComodule::cofree_comodule(coalgebra, 0, UniGrading(1), UniGrading(5), ()));

        let mut map: GradedLinearMap<UniGrading, F2, FlatMatrix<F2>> =
            GradedLinearMap::zero(&domain.space, &codomain.space);
        map.maps.get_mut(&UniGrading(1)).unwrap().data[0] = F2::one();

        let morphism = kComoduleMorphism {
            domain: domain.clone(),
            codomain: codomain.clone(),
            map,
        };

        assert_eq!(
            morphism.get_structure_lines(),
            [(0, 0, 1, "h_0".to_string())]
        );
    }

    #[test]
    fn test_zero_morphism() {
        let coalgebra = Arc::new(A0_coalgebra());
        let comodule = Arc::new(kComodule::fp_comodule(coalgebra, UniGrading::zero()));

        let zero_morphism = kComoduleMorphism::zero_morphism(comodule.clone());

        // Assertions
        assert_eq!(zero_morphism.map.maps.len(), 1); // Assuming `is_zero` verifies all entries are zero
        assert_eq!(zero_morphism.domain.space, GradedVectorSpace::new());
    }
}
