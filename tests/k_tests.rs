// TODO : Only generators get checked here

#[cfg(test)]
mod tests {
    use std::{i32, sync::Arc};

    use algebra::{matrices::flat_matrix::FlatMatrix, rings::finite_fields::F2};
    use comodules::{export::SSeq, grading::{grading::Grading, unigrading::UniGrading}, k_comodule::{kcoalgebra::{A0_coalgebra, kCoalgebra}, kcomodule::kComodule}, resolution::resolve_by_s::Resolution, traits::Coalgebra};
    use itertools::Itertools;

    #[allow(dead_code)]
    fn generate_sseq_jsons() {
        {
            let coalgebra = A0_coalgebra();

            let fp = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());

            let mut res: Resolution<UniGrading, kCoalgebra<UniGrading, F2, FlatMatrix<F2>>> =
                Resolution::new(coalgebra, fp);

            res.resolve_to_s(20, UniGrading(20));

            let a0_sseq = res.generate_sseq("A0");
            a0_sseq.save_to_json("./A(0)sseq.json").unwrap();
        }

        {
            let input = include_str!("../examples/direct/A(1).txt");
            let coalgebra = kCoalgebra::parse(input, UniGrading::infty()).unwrap().0;

            let fp = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());

            let mut res: Resolution<UniGrading, kCoalgebra<UniGrading, F2, FlatMatrix<F2>>> =
                Resolution::new(coalgebra, fp);

            res.resolve_to_s(20, UniGrading(20));

            let p = res.generate_sseq("A(1)");
            p.save_to_json("./A(1)sseq.json").unwrap();
        }

        {
            let input = include_str!("../examples/polynomial/A(2).txt");

            let coalgebra =
                kCoalgebra::parse(input, UniGrading::infty() - UniGrading(10))
                    .unwrap()
                    .0;

            let fp = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());

            let mut res: Resolution<UniGrading, kCoalgebra<UniGrading, F2, FlatMatrix<F2>>> =
                Resolution::new(coalgebra, fp);

            res.resolve_to_s(20, UniGrading(20));
            dbg!(&res);

            let p = res.generate_sseq("A(2)");
            p.save_to_json("./A(2)sseq.json").unwrap();
        }

        {
            let input = include_str!("../examples/polynomial/A.txt");

            let coalgebra = kCoalgebra::parse(input, UniGrading(40)).unwrap().0;

            let fp = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());

            let mut res: Resolution<UniGrading, kCoalgebra<UniGrading, F2, FlatMatrix<F2>>> =
                Resolution::new(coalgebra, fp);

            res.resolve_to_s(25, UniGrading(40));
            dbg!(&res);

            let p = res.generate_sseq("A");
            p.save_to_json("./Asseq.json").unwrap();
        }
    }

    #[test]
    fn test_a0_resolution() {
        let coalgebra = A0_coalgebra();

        let fp = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());

        let mut res: Resolution<UniGrading, kCoalgebra<UniGrading, F2, FlatMatrix<F2>>> =
            Resolution::new(coalgebra, fp);

        res.resolve_to_s(4, UniGrading(10));

        let sseq = res.generate_sseq("A0");
        let page = sseq.pages[0].clone();

        let sorted_gens: Vec<(usize, usize, Vec<i64>)> = page
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
    }

    #[test]
    fn test_a0_consistent_page() {
        let coalgebra = A0_coalgebra();

        let fp = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());

        let mut res: Resolution<UniGrading, kCoalgebra<UniGrading, F2, FlatMatrix<F2>>> =
            Resolution::new(coalgebra, fp);

        res.resolve_to_s(20, UniGrading(20));

        let sseq = res.generate_sseq("A0");
        let comp_sseq: SSeq = serde_json::from_str(include_str!("./A(0).json")).unwrap();

        assert_eq!(sseq.pages[0].generators, comp_sseq.pages[0].generators);
    }

    #[test]
    fn test_a1_resolution() {
        let input = include_str!("../examples/direct/A(1).txt");
        let coalgebra = kCoalgebra::parse(input, UniGrading::infty()).unwrap().0;

        let fp = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());

        let mut res: Resolution<UniGrading, kCoalgebra<UniGrading, F2, FlatMatrix<F2>>> =
            Resolution::new(coalgebra, fp);

        res.resolve_to_s(20, UniGrading(20));

        let p = res.generate_sseq("A(1)");
        let comp_p: SSeq = serde_json::from_str(include_str!("./A(1).json")).unwrap();
        assert_eq!(p.pages[0].generators, comp_p.pages[0].generators);
    }

    #[test]
    fn test_a2_resolution_direct() {
        let input = include_str!("../examples/direct/A(2).txt");
        let coalgebra = kCoalgebra::parse(input, UniGrading::infty()).unwrap().0;

        let fp = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());

        let mut res: Resolution<UniGrading, kCoalgebra<UniGrading, F2, FlatMatrix<F2>>> =
            Resolution::new(coalgebra, fp);

        res.resolve_to_s(20, UniGrading(20));

        let p = res.generate_sseq("A(2)");
        let comp_p: SSeq = serde_json::from_str(include_str!("./A(2).json")).unwrap();
        assert_eq!(p.pages[0].generators, comp_p.pages[0].generators);
    }

    #[test]
    fn test_a2_resolution_poly() {
        let input = include_str!("../examples/polynomial/A(2).txt");

        let coalgebra =
            kCoalgebra::parse(input, UniGrading::infty() - UniGrading(10))
                .unwrap()
                .0;

        let fp = kCoalgebra::basering_comodule(&coalgebra, UniGrading::zero());

        let mut res: Resolution<UniGrading, kCoalgebra<UniGrading, F2, FlatMatrix<F2>>> =
            Resolution::new(coalgebra, fp);

        res.resolve_to_s(20, UniGrading(20));

        let p = res.generate_sseq("A(2)");
        let comp_p: SSeq = serde_json::from_str(include_str!("./A(2).json")).unwrap();
        assert_eq!(p.pages[0].generators, comp_p.pages[0].generators);
    }
}
