#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use comodules::{
        comodule::{
            traits::Comodule,
            kcoalgebra::{kCoalgebra, A0_coalgebra},
            kcomodule::kComodule,
            kmorphism::kComoduleMorphism,
        },
        linalg::{field::F2, graded::UniGrading, matrix::FieldMatrix},
        resolution::Resolution,
    };
    use itertools::Itertools;

    #[test]
    fn test_a0_resolution() {
        let coalgebra = Arc::new(A0_coalgebra());

        let fp = kComodule::fp_comodule(coalgebra);

        let mut res: Resolution<
            UniGrading,
            kComodule<UniGrading, F2, FieldMatrix<F2>>,
            kComoduleMorphism<UniGrading, F2, FieldMatrix<F2>>,
        > = Resolution::new(fp);

        res.resolve_to_s(4, 10);

        let page = res.generate_page();

        let sorted_gens: Vec<(usize, usize, Vec<i32>)> = page
            .generators
            .iter()
            .map(|x| (x.0, x.1, x.2.clone()))
            .sorted_by_key(|x| x.0 * 1000 + x.2[0] as usize)
            .collect();
        assert_eq!(
            sorted_gens,
            vec![
                (0, 0, vec![0]),
                (1, 0, vec![1]),
                (2, 0, vec![2]),
                (3, 0, vec![3]),
                (4, 0, vec![4]),
            ]
        );

        let sorted_lines: Vec<((usize, usize), (usize, usize), usize, String)> = page
            .structure_lines
            .iter()
            .map(|x| (x.clone()))
            .sorted_by_key(|x| x.0 .0 * 1000 + x.0 .1)
            .collect();
        assert_eq!(
            sorted_lines,
            vec![
                ((0, 0), (1, 0), 1, "0".to_string()),
                ((1, 0), (2, 0), 1, "0".to_string()),
                ((2, 0), (3, 0), 1, "0".to_string()),
                ((3, 0), (4, 0), 1, "0".to_string()),
            ]
        );
    }

    #[test]
    fn test_a1_resolution() {
        let input = include_str!("../examples/kcoalgebras/A(1).txt");
        let coalgebra = Arc::new(kCoalgebra::parse(input).unwrap().0);

        let fp = kComodule::fp_comodule(coalgebra);

        let mut res: Resolution<
            UniGrading,
            kComodule<UniGrading, F2, FieldMatrix<F2>>,
            kComoduleMorphism<UniGrading, F2, FieldMatrix<F2>>,
        > = Resolution::new(fp);

        res.resolve_to_s(10, 10);

        let _ = res.generate_page();
    }
}
