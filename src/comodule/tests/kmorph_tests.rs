#[cfg(test)]
mod tests {
    use std::{collections::HashMap, sync::Arc};

    use crate::{
        comodule::{
            comodule::{Comodule, ComoduleMorphism},
            kcomodule::{kComodule, A0_coalgebra},
            kmorphism::kComoduleMorphism,
        },
        linalg::{field::{Field, F2}, graded::{GradedLinearMap, GradedVectorSpace}, matrix::{FieldMatrix, Matrix}},
    };

    #[test]
    fn test_inject_codomain_to_cofree() {
        let coalgebra = Arc::new(A0_coalgebra());
        let comodule = Arc::new(kComodule::fp_comodule(coalgebra.clone()));

        let morphism = kComoduleMorphism::zero_morphism(comodule);

        let cofree_morphism = morphism.inject_codomain_to_cofree(5);

        let comp = kComodule::cofree_comodule(coalgebra, 0, 0, 5);

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
        
        let domain = Arc::new(kComodule::fp_comodule(coalgebra.clone()));

        let codomain = Arc::new(kComodule::cofree_comodule(coalgebra, 0, 0, 5));

        let mut map: GradedLinearMap<i32, F2, FieldMatrix<F2>> = GradedLinearMap::zero(&domain.space, &codomain.space);
        map.maps.get_mut(&0).unwrap().data[0][0] = F2::one();

        let morphism = kComoduleMorphism {
            domain: domain.clone(),
            codomain: codomain.clone(),
            map,
        };

        let cokernel_morphism = morphism.cokernel();

        // Assertions
        assert_eq!(cokernel_morphism.domain, morphism.codomain);
        
        let expected_map: GradedLinearMap<i32, F2, FieldMatrix<F2>> = GradedLinearMap::from(HashMap::from([
            (0, FieldMatrix::zero(1, 0)),
            (1, FieldMatrix::identity(1))
        ]));
        assert_eq!(cokernel_morphism.map, expected_map); 
    }

    #[test]
    fn test_zero_morphism() {
        let coalgebra = Arc::new(A0_coalgebra());
        let comodule = Arc::new(kComodule::fp_comodule(coalgebra));

        let zero_morphism = kComoduleMorphism::zero_morphism(comodule.clone());

        // Assertions
        assert_eq!(zero_morphism.map.maps.len(), 1); // Assuming `is_zero` verifies all entries are zero
        assert_eq!(zero_morphism.domain.space, GradedVectorSpace::new());
    }
}
