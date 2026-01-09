#[cfg(test)]
mod tests {
    use algebra::{
        matrices::flat_matrix::FlatMatrix,
        rings::finite_fields::{F2, Fp},
    };
    use itertools::Itertools;

    use crate::{
        grading::{grading::Grading, unigrading::UniGrading},
        k_comodule::kcoalgebra::{A0_coalgebra, kCoalgebra},
    };

    #[test]
    fn test_a0() {
        let input = include_str!("../../../examples/direct/A(0).txt");

        let (kcoalg, _) =
            kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input, UniGrading::infty())
                .unwrap();

        assert_eq!(kcoalg.coaction, A0_coalgebra().coaction);
        assert_eq!(kcoalg.space, A0_coalgebra().space);
    }

    #[test]
    fn test_a2_consistency_direct() {
        let mut comps = vec![];
        for _ in 0..10 {
            let input = include_str!("../../../examples/direct/A(2).txt");

            let (kcoalg, _) =
                kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input, UniGrading::infty())
                    .unwrap();
            comps.push(kcoalg);
        }
        assert!(comps.iter().all_equal())
    }

    #[test]
    fn test_a2_consistency_poly() {
        let mut comps = vec![];
        for _ in 0..10 {
            let input = include_str!("../../../examples/polynomial/A(2).txt");

            let (kcoalg, _) = kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(
                input,
                UniGrading::infty() - UniGrading(10),
            )
            .unwrap();
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
}
