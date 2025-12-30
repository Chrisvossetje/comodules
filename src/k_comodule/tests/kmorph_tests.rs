#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use ahash::HashMap;
    use algebra::{
        matrices::flat_matrix::FlatMatrix, matrix::Matrix, ring::CRing, rings::finite_fields::F2,
    };

    use crate::{grading::{grading::Grading, unigrading::UniGrading}, k_comodule::{graded_space::GradedLinearMap, kcoalgebra::{A0_coalgebra, kCoalgebra}, kcomodule::{kCofreeComodule, kComodule}, kmorphism::kComoduleMorphism}, traits::{Coalgebra, CofreeComodule}};
    use crate::traits::ComoduleMorphism;

    #[test]
    fn test_inject_codomain_to_cofree() {
        let coalgebra = Arc::new(A0_coalgebra());
        let comodule = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());

        let (_cofree_morphism, cofree_comod) =
            kComoduleMorphism::inject_codomain_to_cofree(&coalgebra, &comodule, UniGrading(5));

        let comp = coalgebra.cofree_comodule(0, UniGrading(0), UniGrading(5), ());

        // Assertions
        assert_eq!(cofree_comod, comp);
    }

    #[test]
    fn test_cokernel() {
        let coalgebra = Arc::new(A0_coalgebra());

        let domain = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());

        let codomain =
            coalgebra.cofree_comodule(0, UniGrading(0), UniGrading(5), ());

        // Manually create a zero map
        let mut maps = HashMap::default();
        maps.insert(
            UniGrading(0),
            FlatMatrix::zero(
                domain.space.dimension_in_grade(&UniGrading(0)),
                codomain.space.dimension_in_grade(&UniGrading(0)),
            ),
        );
        let mut map: GradedLinearMap<UniGrading, F2, FlatMatrix<F2>> = GradedLinearMap::from(maps);

        // Set one element to F2::one()
        map.maps.get_mut(&UniGrading(0)).unwrap().set(0, 0, F2::one());

        let morphism = kComoduleMorphism { map };

        let (cokernel_morphism, _cokernel_comodule) = morphism.cokernel(&coalgebra, &codomain);

        // Assertions - test that cokernel morphism is created
        assert!(cokernel_morphism.map.maps.len() > 0);
    }

    #[test]
    fn test_zero_morphism() {
        let coalgebra = Arc::new(A0_coalgebra());
        let comodule = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());

        let zero_morphism = kComoduleMorphism::zero_morphism(&comodule);

        // Assertions
        assert_eq!(zero_morphism.map.maps.len(), 1); // Should have one grade with zero dimensions
    }
}
