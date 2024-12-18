#[cfg(test)]
mod tests {
    use std::{collections::HashMap, sync::Arc};

    use crate::{
        comodule::{
            comodule::{Comodule, ComoduleMorphism},
            kcoalgebra::A0_coalgebra,
            kcomodule::{kBasisElement, kComodule},
            kmorphism::kComoduleMorphism,
            ktensor::kTensor,
        },
        linalg::{
            field::F2,
            graded::{GradedLinearMap, GradedVectorSpace},
            matrix::{FieldMatrix, Matrix},
        },
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
        let comodule = kComodule::fp_comodule(coalgebra.clone());

        assert_eq!(comodule.space.0.len(), 1);
        assert!(comodule.space.0.contains_key(&0));

        let elements = comodule.space.0.get(&0).unwrap();
        assert_eq!(elements.len(), 1);
        assert_eq!(elements[0].name, "fp");
        assert_eq!(elements[0].generator, false);

        assert!(comodule.coaction.maps.contains_key(&0));
        assert_eq!(comodule.coaction.maps.len(), 1);
        assert_eq!(comodule.coaction.maps[&0], FieldMatrix::<F2>::identity(1));
    }

    // Test for kComodule::direct_sum
    #[test]
    fn test_direct_sum() {
        let coalgebra = Arc::new(A0_coalgebra());
        let mut comodule1 = kComodule::fp_comodule(coalgebra.clone());
        let mut comodule2 = kComodule::cofree_comodule(coalgebra.clone(), 0, 0, 4);

        assert_eq!(comodule2.tensor.dimensions.get(&1).unwrap(), &2);

        comodule1.direct_sum(&mut comodule2);

        assert_eq!(comodule1.space.0.get(&0).unwrap().len(), 2);
        assert_eq!(comodule1.space.0.get(&1).unwrap().len(), 1);

        let elements = &comodule1.space.0[&0];
        assert_eq!(elements[0].name, "fp");
        assert_eq!(elements[1].name, "1");

        let elements = &comodule1.space.0[&1];
        assert_eq!(elements[0].name, "xi1");

        assert_eq!(comodule2.tensor.dimensions[&1], 2);

        let dims = &comodule1.tensor.dimensions;
        assert_eq!(dims.get(&0), Some(&2));
        assert_eq!(dims.get(&1), Some(&2));
        comodule1.tensor.is_correct();
    }

    // Test for kComodule::get_generators
    #[test]
    fn test_get_generators() {
        let mut space_map = HashMap::new();
        space_map.insert(
            0,
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
            tensor: kTensor::new(),
        };

        let generators = comodule.get_generators();

        assert_eq!(generators.len(), 1);
        assert_eq!(generators[0].1, 0);
        assert_eq!(generators[0].2, Some("gen1".to_string()));
    }

    #[test]
    fn test_kcomod_cokernel() {
        let coalgebra = Arc::new(A0_coalgebra());

        let fp = Arc::new(kComodule::fp_comodule(coalgebra));

        let initial = kComoduleMorphism::zero_morphism(fp);
    }
}
