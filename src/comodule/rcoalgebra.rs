use crate::{
    basiselement::kBasisElement,
    comodule::{kcoalgebra::kCoalgebra, rcomodule::RCoalgebra},
    grading::{Grading, UniGrading},
    linalg::{
        field::{F2, Field}, flat_matrix::FlatMatrix, matrix::RModMorphism, ring::{CRing, UniPolRing}
    },
    module::{module::GradedModule, morphism::GradedModuleMap},
    tensor::Tensor,
};


impl<G: Grading, F: Field> RCoalgebra<G, F> {
    pub fn set_generator(&mut self) -> Result<(), &str> {
        let grade_zero = self.space.0.get_mut(&G::zero());
        if let Some(basis) = grade_zero {
            if basis.len() == 1 {
                if let Some(basis_element) = basis.first_mut() {
                    basis_element.0.generator = true;
                }
            } else {
                Err("Coalgebra is not connected, no unique generator found in grade 0")?;
            }
        } else {
            Err("Grade 0 not found in space")?;
        }
        Ok(())
    }
}



#[allow(non_snake_case)]
pub fn A1_C() -> RCoalgebra<UniGrading, F2> {
    let input = include_str!("../../examples/direct/A(1)_C.txt");
    let coalg = RCoalgebra::parse(input, UniGrading::infty()).unwrap().0;

    coalg
}

#[allow(non_snake_case)]
pub fn A0_C() -> RCoalgebra<UniGrading, F2> {
    let mut space = GradedModule::default();
    space.0.insert(
        UniGrading(0),
        vec![(
            kBasisElement {
                name: "1".to_owned(),
                generator: true,
                primitive: None,
                generated_index: 0,
            },
            UniGrading(0),
            None,
        )],
    );

    space.0.insert(
        UniGrading(1),
        vec![(
            kBasisElement {
                name: "xi1".to_owned(),
                generator: false,
                primitive: Some(0),
                generated_index: 0,
            },
            UniGrading(0),
            None,
        )],
    );

    let tensor = Tensor::generate(&space, &space);

    let mut coaction = GradedModuleMap::default();
    coaction.maps.insert(
        UniGrading(0),
        FlatMatrix {
            data: vec![UniPolRing(F2::one(), 0)],
            domain: 1,
            codomain: 1,
        },
    );

    coaction.maps.insert(
        UniGrading(1),
        FlatMatrix {
            data: vec![UniPolRing(F2::one(), 0), UniPolRing(F2::one(), 0)],
            domain: 1,
            codomain: 2,
        },
    );

    RCoalgebra {
        space,
        coaction,
        tensor,
    }
}

pub fn tensor_k_coalgebra(
    coalgebra: kCoalgebra<UniGrading, F2, FlatMatrix<F2>>,
) -> RCoalgebra<UniGrading, F2> {

    let (space, coaction, tensor) = (coalgebra.space, coalgebra.coaction, coalgebra.tensor);

    let space = space
        .0
        .into_iter()
        .map(|x| {
            let gr = UniGrading(x.0 .0);
            let module: Vec<(kBasisElement, UniGrading, Option<usize>)> =
                x.1.into_iter().map(|y| (y, UniGrading(0), None)).collect();
            (gr, module)
        })
        .collect();

    let coaction_maps = coaction.maps.into_iter().map(|(gr, m)| {        
        let mut new_m = FlatMatrix::zero(m.domain, m.codomain);
        for x in 0..m.domain {
            for y in 0..m.codomain {
                let el = m.get(x, y);
                new_m.set(x, y, UniPolRing(el, 0));
            }
        }
        (gr, new_m)
    }).collect();

    if cfg!(debug_assertions) {
        if !tensor.is_correct() {
            panic!("Tensor is not correct");
        }
    }

    RCoalgebra {
        space: GradedModule(space),
        coaction: GradedModuleMap {
            maps: coaction_maps,
        },
        tensor 
    }
}
