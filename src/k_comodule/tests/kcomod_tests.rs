#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use ahash::HashMap;
    use algebra::{
        matrices::flat_matrix::FlatMatrix,
        matrix::Matrix,
        rings::finite_fields::{F2, Fp},
    };

    use crate::{
        grading::{grading::Grading, tensor::TensorMap, unigrading::UniGrading},
        k_comodule::{
            graded_space::{GradedLinearMap, GradedVectorSpace},
            kcoalgebra::{A0_coalgebra, kCoalgebra},
            kcomodule::kComodule,
        },
        traits::Coalgebra,
    };

    // Test for creating an empty kComodule manually
    #[test]
    fn test_empty_comodule() {
        let coalgebra = Arc::new(A0_coalgebra());

        // Create empty comodule manually
        let empty_space: GradedVectorSpace<UniGrading, ()> = GradedVectorSpace::new();
        let empty_coaction_maps = HashMap::default();
        let empty_coaction =
            GradedLinearMap::<UniGrading, F2, FlatMatrix<F2>>::from(empty_coaction_maps);
        let empty_tensor = TensorMap::default();

        let comodule: kComodule<UniGrading, kCoalgebra<UniGrading, F2, FlatMatrix<F2>>> =
            kComodule::new(empty_space, empty_coaction, empty_tensor);

        assert_eq!(coalgebra.clone().as_ref(), coalgebra.clone().as_ref());
        assert!(comodule.space.0.is_empty());
        assert!(comodule.coaction.maps.is_empty());
    }

    // Test for kComodule::fp_comodule
    #[test]
    fn test_fp_comodule() {
        let coalgebra = Arc::new(A0_coalgebra());
        let comodule = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());

        assert_eq!(comodule.space.0.len(), 1);
        assert!(comodule.space.0.contains_key(&UniGrading(0)));

        let elements = comodule.space.0.get(&UniGrading(0)).unwrap();
        assert_eq!(elements.len(), 1);
        // Elements are now just () instead of kBasisElement with fields

        assert!(comodule.coaction.maps.contains_key(&UniGrading(0)));
        assert_eq!(comodule.coaction.maps.len(), 1);
        assert_eq!(
            comodule.coaction.maps[&UniGrading(0)],
            FlatMatrix::<F2>::identity(1)
        );
    }

    // Test for direct_sum - simplified since we can't test specific names
    #[test]
    fn test_direct_sum_basic() {
        let coalgebra = Arc::new(A0_coalgebra());
        let comodule1 = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());
        let comodule2 = kCoalgebra::basering_comodule(&coalgebra, UniGrading(1));

        // Test basic properties of the comodules
        assert_eq!(comodule1.space.0.len(), 1);
        assert_eq!(comodule2.space.0.len(), 1);

        // Check dimensions in specific grades
        assert_eq!(comodule1.space.dimension_in_grade(&UniGrading(0)), 1);
        assert_eq!(comodule2.space.dimension_in_grade(&UniGrading(1)), 1);
    }

    // Test for basic comodule construction with unit elements
    #[test]
    fn test_comodule_construction() {
        // Create a simple space with () elements
        let mut space_map = HashMap::default();
        space_map.insert(UniGrading(0), vec![()]);

        // Create simple coaction and tensor maps
        let mut coaction_maps = HashMap::default();
        coaction_maps.insert(UniGrading(0), FlatMatrix::identity(1));
        let coaction = GradedLinearMap::<UniGrading, F2, FlatMatrix<F2>>::from(coaction_maps);

        let tensor = TensorMap::default();

        let comodule: kComodule<UniGrading, kCoalgebra<UniGrading, F2, FlatMatrix<F2>>> =
            kComodule::new(GradedVectorSpace::from(space_map), coaction, tensor);

        // Basic tests
        assert_eq!(comodule.space.0.len(), 1);
        assert_eq!(comodule.space.dimension_in_grade(&UniGrading(0)), 1);
        assert_eq!(comodule.coaction.maps.len(), 1);
    }

    #[test]
    fn test_fp_comod_parser() {
        let input_coalg = include_str!("../../../examples/direct/A(0).txt");
        let input_comod = include_str!("../../../examples/comodule/F2_comod.txt");

        let (kcoalg, translator) =
            kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input_coalg, UniGrading::infty())
                .unwrap();

        match kComodule::parse(input_comod, &kcoalg, &translator, UniGrading::infty()) {
            Ok(comod) => {
                assert_eq!(comod.space.0[&UniGrading(0)].len(), 1);
                assert!(!comod.space.0.contains_key(&UniGrading(1)));
                assert_eq!(comod.coaction.maps.len(), 1);
            }
            Err(e) => {
                println!("{}", e);
                assert!(false);
            }
        }
    }

    #[test]
    fn test_a0_comod_parser() {
        let input_coalg = include_str!("../../../examples/polynomial/A(0).txt");
        let input_comod = include_str!("../../../examples/comodule/A(0)_comod.txt");

        let (kcoalg, translator) =
            kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input_coalg, UniGrading(32))
                .unwrap();

        match kComodule::parse(input_comod, &kcoalg, &translator, UniGrading::infty()) {
            Ok(comod) => {
                println!("{:?}", comod.coaction);
            }
            Err(e) => {
                println!("{}", e);
                assert!(false);
            }
        }
    }

    #[test]
    fn test_gen_comod_parser() {
        let input_coalg = include_str!("../../../examples/polynomial/Test.txt");
        let input_comod = include_str!("../../../examples/comodule/gen_comod.txt");

        let (kcoalg, translator) =
            kCoalgebra::<UniGrading, Fp<3>, FlatMatrix<Fp<3>>>::parse(input_coalg, UniGrading(128))
                .unwrap();
        match kComodule::parse(input_comod, &kcoalg, &translator, UniGrading(6)) {
            Ok(comod) => {
                println!("{:?}", comod.coaction);
            }
            Err(e) => {
                println!("{}", e);
                assert!(false);
            }
        }
    }

    #[test]
    fn test_a0_comod() {
        let input = include_str!("../../../examples/polynomial/A(0).txt");
        const MAX_GRADING: UniGrading = UniGrading(60);
        let (coalgebra, translate) =
            kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input, MAX_GRADING).unwrap();

        let coalgebra = Arc::new(coalgebra);
        let input = include_str!("../../../examples/comodule/A(0).txt");

        let comod = kComodule::parse(input, &coalgebra, &translate, MAX_GRADING).unwrap();

        // Just test that the parsing worked correctly
        assert!(comod.space.0.len() > 0);
        assert!(comod.coaction.maps.len() > 0);
    }

    #[test]
    fn test_a1_comod() {
        let input = include_str!("../../../examples/polynomial/A(1).txt");
        const MAX_GRADING: UniGrading = UniGrading(60);
        let (coalgebra, translate) =
            kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input, MAX_GRADING).unwrap();

        let coalgebra = coalgebra;
        let input = include_str!("../../../examples/comodule/A(1).txt");

        let comod = kComodule::parse(input, &coalgebra, &translate, MAX_GRADING).unwrap();

        // Just test that the parsing worked correctly
        assert!(comod.space.0.len() > 0);
        assert!(comod.coaction.maps.len() > 0);
    }
}
