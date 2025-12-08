// #[cfg(test)]
// mod tests {
//     use std::{i32, sync::Arc};

//     use algebra::{matrices::flat_matrix::FlatMatrix, rings::finite_fields::F2};
//     use comodules::{
//         comodule::{
//             kcoalgebra::{A0_coalgebra, kCoalgebra}, kcomodule::kComodule, traits::Comodule
//         }, export::SSeq, grading::{Grading, UniGrading}, resolution::Resolution
//     };
//     use itertools::Itertools;

//     #[test]
//     fn test_z2_coalgebra_generation() {
//         let coalg = Arc::new(Z2::generate_coalgebra::<F2>().unwrap());
        
//         let fp = kComodule::fp_comodule(coalg, Z2::zero());

//         let mut res = Resolution::new(fp);

//         res.resolve_to_s(20, Z2::infty());

//         let _p = res.generate_sseq("Z2");
//     }
// }
