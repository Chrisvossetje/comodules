use std::marker::PhantomData;

use ahash::HashMap;
use algebra::{
    abelian::Abelian, field::Field, matrices::flat_matrix::FlatMatrix, matrix::Matrix, ring::CRing, rings::finite_fields::F2
};
use deepsize::DeepSizeOf;
use itertools::Itertools;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::{
    basiselement::kBasisElement, comodule::{kcomodule::{kCofreeComodule, kComodule}, kmorphism::kComoduleMorphism, traits::Coalgebra}, graded_space::{BasisIndex, GradedLinearMap, GradedVectorSpace}, grading::{Grading, UniGrading}, tensor::TensorMap
};

use serde::{Deserialize, Serialize};

// #[derive(Debug, Clone, PartialEq, Deserialize, Serialize, DeepSizeOf)]
// #[allow(non_camel_case_types)]
// pub struct kCoalgebra<G: Grading, F: Field, M: Matrix<F>> {
//     pub space: GradedVectorSpace<G, kBasisElement>,
//     pub coaction: GradedLinearMap<G, F, M>,
//     // pub coaction: HashMap<BasisIndex<G>, Vec<(BasisIndex<G>,BasisIndex<G>, F)>>,
//     pub tensor: TensorMap<G>,
// }

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, DeepSizeOf)]
#[allow(non_camel_case_types)]
pub struct kCoalgebra<G: Grading, F: Field, M: Matrix<F>> {
    pub space: GradedVectorSpace<G, kBasisElement>,
    pub coaction: HashMap<BasisIndex<G>, Vec<(BasisIndex<G>,BasisIndex<G>, F)>>,
    pub(crate) _p: PhantomData<M>,
}

impl<G: Grading, F: Field, M: Abelian<F>> Coalgebra<G> for kCoalgebra<G,F,M> {
    type BaseRing = F;
    type RingMorph = M;
    type Comod = kComodule<G,F,M>; 
    type CofMod = kCofreeComodule<G,F,M>;
    type ComodMorph = kComoduleMorphism<G,F,M>;
    
    fn size_in_degree(&self, g: G) -> usize {
        self.space.dimension_in_grade(&g)
    }
    
    fn coaction(&self, i: BasisIndex<G>) -> &[(BasisIndex<G>, BasisIndex<G>, Self::BaseRing)] {
        self.coaction.get(&i).unwrap().as_slice()
    }
}

impl<G: Grading, F: Field, M: Matrix<F>> kCoalgebra<G, F, M> {
    pub fn set_primitives(&mut self) {
        let mut primitive_index = 0;

        for (id, c) in &self.coaction {
            if c.len() == 2 {
                let el = self.space.0.get_mut(&id.0).unwrap().get_mut(id.1 as usize).unwrap();
                el.primitive = Some(primitive_index);
                primitive_index += 1;
            }
        }
    }

    pub fn set_generator(&mut self) -> Result<(), &str> {
        let grade_zero = self.space.0.get_mut(&G::zero());
        if let Some(basis) = grade_zero {
            if basis.len() == 1 {
                if let Some(basis_element) = basis.first_mut() {
                    basis_element.generator = true;
                }
            } else {
                Err("Coalgebra is not connected, no unique generator found in grade 0")?;
            }
        } else {
            Err("Grade 0 not found in space")?;
        }
        Ok(())
    }

    pub fn reduce(&mut self) {
        // reduce_helper(&mut self.coaction, &mut self.tensor);
    }
}

// TODO :

// #[allow(non_snake_case)]
// pub fn A0_coalgebra() -> kCoalgebra<UniGrading, F2, FlatMatrix<F2>> {
//     let mut space = GradedVectorSpace::new();
//     space.0.insert(
//         UniGrading(0),
//         vec![kBasisElement {
//             name: "1".to_owned(),
//             generator: true,
//             primitive: None,
//             generated_index: 0,
//         }],
//     );

//     space.0.insert(
//         UniGrading(1),
//         vec![kBasisElement {
//             name: "xi1".to_owned(),
//             generator: false,
//             primitive: Some(0),
//             generated_index: 0,
//         }],
//     );

//     let mut dimensions = HashMap::default();
//     dimensions.insert(UniGrading(0), 1);
//     dimensions.insert(UniGrading(1), 2);

//     let mut construct = HashMap::default();
//     let mut first_entry = HashMap::default();
//     first_entry.insert((UniGrading(0), 0), (UniGrading(0), 0));
//     first_entry.insert((UniGrading(1), 0), (UniGrading(1), 1));

//     let mut second_entry = HashMap::default();
//     second_entry.insert((UniGrading(0), 0), (UniGrading(1), 0));

//     construct.insert((UniGrading(0), 0), first_entry);
//     construct.insert((UniGrading(1), 0), second_entry);

//     let mut deconstruct = HashMap::default();
//     deconstruct.insert((UniGrading(0), 0), ((UniGrading(0), 0), (UniGrading(0), 0)));
//     deconstruct.insert((UniGrading(1), 1), ((UniGrading(1), 0), (UniGrading(0), 0)));
//     deconstruct.insert((UniGrading(1), 0), ((UniGrading(0), 0), (UniGrading(1), 0)));

//     let tensor = TensorMap {
//         construct,
//         deconstruct,
//         dimensions,
//     };

//     let mut first_mat = FlatMatrix::zero(1, 1);
//     first_mat.set(0, 0, F2::one());

//     let mut second_mat = FlatMatrix::zero(1, 2);
//     second_mat.set(0, 0, F2::one());
//     second_mat.set(0, 1, F2::one());

//     let mut coaction = GradedLinearMap::empty();
//     coaction.maps.insert(UniGrading(0), first_mat);
//     coaction.maps.insert(UniGrading(1), second_mat);

//     kCoalgebra {
//         space,
//         coaction,
//         tensor,
//     }
// }

// pub fn reduce_helper<G: Grading, F: Field, M: Matrix<F>>(
//     coaction: &mut GradedLinearMap<G, F, M>,
//     tensor: &mut TensorMap<G>,
// ) {
//     // Left vec goes from old t_id -> new t_id, usize is dimension
//     let mapping: HashMap<G, (Vec<Option<usize>>, usize)> = tensor
//         .dimensions
//         .par_iter()
//         .map(|(grade, _)| {
//             let dimension = tensor.get_dimension(grade);
//             let mut new_dimension: usize = 0;
//             let mut to_new = vec![];
//             for codomain in 0..dimension {
//                 if coaction.maps.get(grade).unwrap().is_row_non_zero(codomain) {
//                     to_new.push(Some(new_dimension));
//                     new_dimension += 1;
//                 } else {
//                     to_new.push(None);
//                 }
//             }
//             (*grade, (to_new, new_dimension))
//         })
//         .collect();

//     // This should happen before the others !
//     tensor.dimensions = mapping
//         .iter()
//         .map(|(g, (_, new_dim))| (*g, *new_dim))
//         .collect();

//     let new_coaction = coaction
//         .maps
//         .par_iter()
//         .map(|(grade, matrix)| {
//             let new_codom = tensor.dimensions.get(grade).unwrap();
//             let mut new_matrix = M::zero(matrix.domain(), *new_codom);
//             for (old_row, new_maybe) in mapping.get(grade).unwrap().0.iter().enumerate() {
//                 match new_maybe {
//                     Some(new_row) => {
//                         new_matrix.set_row(*new_row, matrix.get_row(old_row));
//                     }
//                     None => {}
//                 }
//             }
//             (*grade, new_matrix)
//         })
//         .collect();
//     coaction.maps = new_coaction;

//     let new_deconstruct = tensor
//         .deconstruct
//         .iter()
//         .filter_map(
//             |((grade, old_tid), &target)| match mapping.get(grade).unwrap().0[*old_tid as usize] {
//                 Some(new_id) => Some(((*grade, new_id as u32), target)),
//                 None => None,
//             },
//         )
//         .collect();
//     tensor.deconstruct = new_deconstruct;

//     let new_construct = tensor
//         .construct
//         .par_iter()
//         .map(|((m_gr, m_id), map)| {
//             let new_map = map
//                 .iter()
//                 .filter_map(|((a_gr, a_id), (t_gr, old_t_id))| {
//                     match mapping.get(t_gr).unwrap().0[*old_t_id as usize] {
//                         Some(new_tid) => Some(((*a_gr, *a_id), (*t_gr, new_tid as u32))),
//                         None => None,
//                     }
//                 })
//                 .collect();
//             ((*m_gr, *m_id), new_map)
//         })
//         .collect();
//     tensor.construct = new_construct;

//     debug_assert!(
//         tensor.is_correct(),
//         "Tensor is not correct after reducing :("
//     );
// }
