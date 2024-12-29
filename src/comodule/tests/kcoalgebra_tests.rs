#[cfg(test)]
mod tests {
    use std::i32;

    use itertools::Itertools;

    use crate::{
        comodule::kcoalgebra::{kCoalgebra, A0_coalgebra},
        linalg::{field::{Field, F2}, row_matrix::RowMatrix},
    };

    #[test]
    fn test_a0() {
        let input = include_str!("../../../examples/kcoalgebras/A(0).txt");

        let (kcoalg, _) = kCoalgebra::<i32, F2, RowMatrix<F2>>::parse_direct(input).unwrap();

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
            let input = include_str!("../../../examples/kcoalgebras/A(2).txt");

            let (kcoalg, _) = kCoalgebra::<i32, F2, RowMatrix<F2>>::parse_direct(input).unwrap();
            comps.push(kcoalg);
        }
        assert!(comps.iter().all_equal())
    }

    #[test]
    fn test_a2_consistency_poly() {
        let mut comps = vec![];
        for _ in 0..10 {
            let input = include_str!("../../../examples/kcoalgebras/A(2)_gen.txt");

            let (kcoalg, _) = kCoalgebra::<i32, F2, RowMatrix<F2>>::parse_polynomial_hopf_algebra(input, i32::MAX).unwrap();
            comps.push(kcoalg);
        }
        assert!(comps.iter().all_equal())
    }

    #[test]
    fn test_poly_vs_direct_tensor() {
        let input_direct = include_str!("../../../examples/kcoalgebras/A(2).txt");
        let input_poly = include_str!("../../../examples/kcoalgebras/A(2)_gen.txt");

        let (kcoalg_direct, _) = kCoalgebra::<i32, F2, RowMatrix<F2>>::parse_direct(input_direct).unwrap();
        let (kcoalg_poly, _) = kCoalgebra::<i32, F2, RowMatrix<F2>>::parse_polynomial_hopf_algebra(input_poly, i32::MAX).unwrap();

        for grade in kcoalg_direct.tensor.dimensions.keys() {
            assert_eq!(
                  (kcoalg_direct.coaction.maps[grade].domain, kcoalg_direct.coaction.maps[grade].codomain)
                , (kcoalg_poly.coaction.maps[grade].domain, kcoalg_poly.coaction.maps[grade].codomain)
                , "Mismatched dimensions, grade: {}", grade);
            
                let num_non_zeros_direct = kcoalg_direct.coaction.maps[grade].data.iter().map(
                    |v| v.iter().filter(|x| !(*x != &F2::zero())).count()
                ).sum::<usize>();
                let num_non_zeros_poly = kcoalg_poly.coaction.maps[grade].data.iter().map(
                    |v| v.iter().filter(|x| !(*x != &F2::zero())).count()
                ).sum::<usize>();
                assert_eq!(num_non_zeros_direct, num_non_zeros_poly, "Mismatched number of non-zero entries, grade: {}", grade);
        }
    }
}
