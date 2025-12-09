#[cfg(test)]
mod tests {
    use std::{sync::Arc};

    use ahash::HashMap;
    use algebra::{matrices::flat_matrix::FlatMatrix, matrix::Matrix, rings::finite_fields::{F2, Fp}};

    use crate::{
        basiselement::kBasisElement, comodule::{
            kcoalgebra::{A0_coalgebra, kCoalgebra},
            kcomodule::kComodule,
            traits::Comodule,
        }, graded_space::{GradedLinearMap, GradedVectorSpace}, grading::{Grading, UniGrading}, resolution::Resolution, tensor::TensorMap
    };

    // Test for kComodule::zero_comodule
    #[test]
    fn test_zero_comodule() {
        let coalgebra = Arc::new(A0_coalgebra());
        let comodule = kComodule::zero_comodule(coalgebra.clone());

        assert_eq!(
            comodule.coalgebra.clone().as_ref(),
            coalgebra.clone().as_ref()
        );
        assert!(comodule.space.0.is_empty());
        assert!(comodule.coaction.maps.is_empty());
    }

    // Test for kComodule::fp_comodule
    #[test]
    fn test_fp_comodule() {
        let coalgebra = Arc::new(A0_coalgebra());
        let comodule = kComodule::fp_comodule(coalgebra.clone(), UniGrading::zero());

        assert_eq!(comodule.space.0.len(), 1);
        assert!(comodule.space.0.contains_key(&UniGrading(0)));

        let elements = comodule.space.0.get(&UniGrading(0)).unwrap();
        assert_eq!(elements.len(), 1);
        assert_eq!(elements[0].name, "fp");
        assert_eq!(elements[0].generator, false);

        assert!(comodule.coaction.maps.contains_key(&UniGrading(0)));
        assert_eq!(comodule.coaction.maps.len(), 1);
        assert_eq!(comodule.coaction.maps[&UniGrading(0)], FlatMatrix::<F2>::identity(1));
    }

    // Test for kComodule::direct_sum
    #[test]
    fn test_direct_sum() {
        let coalgebra = Arc::new(A0_coalgebra());
        let mut comodule1 = kComodule::fp_comodule(coalgebra.clone(), UniGrading::zero());
        let mut comodule2 = kComodule::cofree_comodule(coalgebra.clone(), 0, UniGrading(0), UniGrading(4), ());

        assert_eq!(comodule2.tensor.dimensions.get(&UniGrading(1)).unwrap(), &2);

        comodule1.direct_sum(&mut comodule2);

        assert_eq!(comodule1.space.0.get(&UniGrading(0)).unwrap().len(), 2);
        assert_eq!(comodule1.space.0.get(&UniGrading(1)).unwrap().len(), 1);

        let elements = &comodule1.space.0[&UniGrading(0)];
        assert_eq!(elements[0].name, "fp");
        assert_eq!(elements[1].name, "1");

        let elements = &comodule1.space.0[&UniGrading(1)];
        assert_eq!(elements[0].name, "xi1");

        assert_eq!(comodule2.tensor.dimensions[&UniGrading(1)], 2);

        let dims = &comodule1.tensor.dimensions;
        assert_eq!(dims.get(&UniGrading(0)), Some(&2));
        assert_eq!(dims.get(&UniGrading(1)), Some(&2));
        comodule1.tensor.is_correct();
    }

    // Test for kComodule::get_generators
    #[test]
    fn test_get_generators() {
        let mut space_map = HashMap::default();
        space_map.insert(
            UniGrading(0),
            vec![
                kBasisElement {
                    name: "gen1".to_string(),
                    generator: true,
                    primitive: None,
                    generated_index: 0,
                },
                kBasisElement {
                    name: "non_gen".to_string(),
                    generator: false,
                    primitive: None,
                    generated_index: 1,
                },
            ],
        );

        let coalgebra = Arc::new(A0_coalgebra());

        let comodule = kComodule {
            coalgebra,
            space: GradedVectorSpace::from(space_map),
            coaction: GradedLinearMap::empty(),
            tensor: TensorMap::default(),
        };

        let generators = comodule.get_generators();

        assert_eq!(generators.len(), 1);
        assert_eq!(generators[0].1, UniGrading(0));
        assert_eq!(generators[0].2, None);
    }

    #[test]
    fn test_fp_comod_parser() {
        let input_coalg = include_str!("../../../examples/direct/A(0).txt");
        let input_comod = include_str!("../../../examples/comodule/F2_comod.txt");

        let (kcoalg, translator) =
            kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input_coalg, UniGrading::infty()).unwrap();

        match kComodule::<UniGrading, F2, FlatMatrix<F2>>::parse(
            input_comod,
            Arc::new(kcoalg),
            &translator,
            UniGrading::infty(),
        ) {
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
            kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input_coalg, UniGrading(32)).unwrap();

        match kComodule::<UniGrading, F2, FlatMatrix<F2>>::parse(
            input_comod,
            Arc::new(kcoalg),
            &translator,
            UniGrading::infty(),
        ) {
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
            kCoalgebra::<UniGrading, Fp<3>, FlatMatrix<Fp<3>>>::parse(input_coalg, UniGrading(128)).unwrap();
        match kComodule::<UniGrading, Fp<3>, FlatMatrix<Fp<3>>>::parse(
            input_comod,
            Arc::new(kcoalg),
            &translator,
            UniGrading(6),
        ) {
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

        let comod = kComodule::parse(input, coalgebra, &translate, MAX_GRADING).unwrap();

        let mut res: Resolution<UniGrading, kComodule<UniGrading, F2, FlatMatrix<F2>>> = Resolution::new(comod);
        res.resolve_to_s(10, UniGrading(1234));
        let sseq = res.generate_sseq("A(1)-comod");

        assert_eq!(sseq.pages[0].generators.len(), 1);
    }

    #[test]
    fn test_a1_comod() {
        let input = include_str!("../../../examples/polynomial/A(1).txt");
        const MAX_GRADING: UniGrading = UniGrading(60);
        let (coalgebra, translate) =
            kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input, MAX_GRADING).unwrap();

        let coalclcone = coalgebra.clone();
        let coalgebra = Arc::new(coalgebra);
        let input = include_str!("../../../examples/comodule/A(1).txt");

        let comod = kComodule::parse(input, coalgebra, &translate, MAX_GRADING).unwrap();
        for n in 0..=6 {
            let n = UniGrading(n);
            println!("{:?}", coalclcone.coaction.maps[&n]);
            println!("{:?}", comod.coaction.maps[&n]);
            println!("");
        }

        let mut res: Resolution<UniGrading, kComodule<UniGrading, F2, FlatMatrix<F2>>> = Resolution::new(comod);
        res.resolve_to_s(10, UniGrading(1234));
        let sseq = res.generate_sseq("A(1)-comod");

        assert_eq!(sseq.pages[0].generators.len(), 1);
    }
}
