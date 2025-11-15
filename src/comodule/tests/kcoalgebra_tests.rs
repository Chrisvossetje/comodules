#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::{
        comodule::kcoalgebra::{A0_coalgebra, kCoalgebra},
        grading::{Grading, UniGrading},
        linalg::{
            field::{F2, Fp}, flat_matrix::FlatMatrix, ring::CRing
        },
    };

    #[test]
    fn test_a0() {
        let input = include_str!("../../../examples/direct/A(0).txt");

        let (kcoalg, _) = kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input, UniGrading::infty()).unwrap();

        assert_eq!(kcoalg.coaction, A0_coalgebra().coaction);

        // HASHMAPS are NOT deterministic SO DON'T COMPARE
        // as construct and deconstruct are dependent on insertion order
        assert_eq!(kcoalg.tensor.dimensions, A0_coalgebra().tensor.dimensions);
        assert_eq!(kcoalg.space, A0_coalgebra().space);
    }

    #[test]
    fn test_a2_consistency_direct() {
        let mut comps = vec![];
        for _ in 0..10 {
            let input = include_str!("../../../examples/direct/A(2).txt");

            let (kcoalg, _) = kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input, UniGrading::infty()).unwrap();
            comps.push(kcoalg);
        }
        assert!(comps.iter().all_equal())
    }

    #[test]
    fn test_a2_consistency_poly() {
        let mut comps = vec![];
        for _ in 0..10 {
            let input = include_str!("../../../examples/polynomial/A(2).txt");

            let (kcoalg, _) =
                kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input, UniGrading::infty() - UniGrading(10)).unwrap();
            comps.push(kcoalg);
        }
        assert!(comps.iter().all_equal())
    }

    #[test]
    fn test_p3_imports() {
        let input = include_str!("../../../examples/polynomial/P(3).txt");
        let res = kCoalgebra::<UniGrading, Fp<3>, FlatMatrix<Fp<3>>>::parse(input, UniGrading(129));

        assert!(res.is_ok());
        let (_, trans) = res.unwrap();
        assert!(trans.len() > 3)
    }

    #[test]
    fn test_poly_vs_direct_tensor() {
        let input_direct = include_str!("../../../examples/direct/A(2).txt");
        let input_poly = include_str!("../../../examples/polynomial/A(2).txt");

        let (kcoalg_direct, _) =
            kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input_direct, UniGrading::infty()).unwrap();
        let (kcoalg_poly, _) =
            kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input_poly, UniGrading::infty() - UniGrading(10)).unwrap();

        for grade in kcoalg_direct.tensor.dimensions.keys() {
            assert_eq!(
                (
                    kcoalg_direct.coaction.maps[grade].domain,
                    kcoalg_direct.coaction.maps[grade].codomain
                ),
                (
                    kcoalg_poly.coaction.maps[grade].domain,
                    kcoalg_poly.coaction.maps[grade].codomain
                ),
                "Mismatched dimensions, grade: {}",
                grade
            );

            let num_non_zeros_direct = kcoalg_direct.coaction.maps[grade]
                .data
                .iter()
                .filter(|x| !(*x != &F2::zero())).count();
            let num_non_zeros_poly = kcoalg_poly.coaction.maps[grade]
                .data
                .iter()
                .filter(|x| !(*x != &F2::zero())).count();
            assert_eq!(
                num_non_zeros_direct, num_non_zeros_poly,
                "Mismatched number of non-zero entries, grade: {}",
                grade
            );
        }
    }
}
