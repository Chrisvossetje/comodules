#[cfg(test)]
mod tests {
    use std::{collections::HashMap, sync::Arc};

    use crate::{comodule::{comodule::{Comodule, ComoduleMorphism}, kcomodule::{kBasisElement, kCoalgebra, kComodule, kComoduleMorphism, kTensor}}, linalg::{field::{Field, F2}, graded::{GradedLinearMap, GradedVectorSpace, UniGrading}, matrix::FieldMatrix}};


    fn simple_kcomod() {

    }


    fn simple_kcomorph() {

    }

    fn A0_coalgebra() -> kCoalgebra<UniGrading, F2> {
        let mut space = GradedVectorSpace::new();
        space.0.insert(0, vec![kBasisElement { 
            name: "0".to_owned(), 
            generator: true, 
            primitive: None, 
            generated_index: 0, 
        }]);

        space.0.insert(1, vec![kBasisElement { 
            name: "xi1".to_owned(), 
            generator: false, 
            primitive: Some(0), 
            generated_index: 0, 
        }]);


        let mut dimensions = HashMap::new();
        dimensions.insert(0, 1);
        dimensions.insert(1, 1);
        
        let mut construct = HashMap::new();
        let mut first_entry = HashMap::new();
        first_entry.insert((0,0), (0, 0));
        first_entry.insert((1,0), (1, 0));

        let mut second_entry = HashMap::new();
        second_entry.insert((0,0), (1, 1));

        construct.insert((0, 0), first_entry);
        construct.insert((1, 0), second_entry);

        let mut deconstruct = HashMap::new();
        deconstruct.insert((0,0), ((0,0), (0,0)));
        deconstruct.insert((1,0), ((1,0), (0,0)));
        deconstruct.insert((1,1), ((0,0), (1,0)));

        let tensor = kTensor { construct, deconstruct, dimensions };

        let mut coaction = GradedLinearMap::empty();
        coaction.maps.insert(0, FieldMatrix {
            data: vec![vec![F2::one()]],
            domain: 1,
            codomain: 1,  
        });
        coaction.maps.insert(1, FieldMatrix {
            data: vec![vec![F2::one()],vec![F2::one()]],
            domain: 1,
            codomain: 2,
        });

        kCoalgebra {
            space,
            coaction,
            tensor,
        }

    }

    #[test]
    fn test_kcomod_cokernel() {
        let coalgebra = Arc::new(A0_coalgebra());

        let fp = Arc::new(kComodule::fp_comodule(coalgebra));

        let initial = kComoduleMorphism::zero_morphism(fp);
        
        dbg!(initial.cokernel());
    }
}