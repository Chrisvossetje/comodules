#[cfg(test)]
mod tests {
    use std::{i32, sync::Arc};

    use comodules::{
        comodule::{
            kcoalgebra::{A0_coalgebra, kCoalgebra}, kcomodule::kComodule, traits::Comodule
        }, export::SSeq, grading::{Grading, UniGrading}, groups::{Group, Z2}, linalg::{field::F2, flat_matrix::FlatMatrix}, resolution::Resolution
    };
    use itertools::Itertools;

    // #[test]
    #[allow(dead_code)]
    fn generate_sseq_jsons() {
        {
            let coalgebra = Arc::new(A0_coalgebra());
            let fp = kComodule::fp_comodule(coalgebra, UniGrading::zero());
            let mut res: Resolution<UniGrading, kComodule<UniGrading, F2, FlatMatrix<F2>>> =
                Resolution::new(fp);

            res.resolve_to_s(20, UniGrading(20));

            let a0_sseq = res.generate_sseq("A0");
            a0_sseq.save_to_json("./A(0)sseq.json").unwrap();
        }

        {
            let input = include_str!("../examples/direct/A(1).txt");
            let coalgebra = Arc::new(kCoalgebra::parse(input, UniGrading::infty()).unwrap().0);

            let fp = kComodule::fp_comodule(coalgebra, UniGrading::zero());

            let mut res: Resolution<UniGrading, kComodule<UniGrading, F2, FlatMatrix<F2>>> =
                Resolution::new(fp);

            res.resolve_to_s(20, UniGrading(20));

            let p = res.generate_sseq("A(1)");
            p.save_to_json("./A(1)sseq.json").unwrap();
        }

        {
            let input = include_str!("../examples/polynomial/A(2).txt");

            let coalgebra = Arc::new(kCoalgebra::parse(input, UniGrading::infty() - UniGrading(10)).unwrap().0);

            let fp = kComodule::fp_comodule(coalgebra, UniGrading::zero());

            let mut res: Resolution<UniGrading, kComodule<UniGrading, F2, FlatMatrix<F2>>> =
                Resolution::new(fp);

            res.resolve_to_s(20, UniGrading(20));
            dbg!(&res);

            let p = res.generate_sseq("A(2)");
            p.save_to_json("./A(2)sseq.json").unwrap();
        }
    }

    #[test]
    fn test_a0_resolution() {
        let coalgebra = Arc::new(A0_coalgebra());

        let fp = kComodule::fp_comodule(coalgebra, UniGrading::zero());

        let mut res: Resolution<UniGrading, kComodule<UniGrading, F2, FlatMatrix<F2>>> =
            Resolution::new(fp);

        res.resolve_to_s(4, UniGrading(10));

        let sseq = res.generate_sseq("A0");
        let page = sseq.pages[0].clone();

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
            .map(|x| x.clone())
            .sorted_by_key(|x| x.0 .0 * 1000 + x.0 .1)
            .collect();
        assert_eq!(
            sorted_lines,
            vec![
                ((0, 0), (1, 0), 1, "h_0".to_string()),
                ((1, 0), (2, 0), 1, "h_0".to_string()),
                ((2, 0), (3, 0), 1, "h_0".to_string()),
                ((3, 0), (4, 0), 1, "h_0".to_string()),
            ]
        );
    }

    #[test]
    fn test_a0_consistent_page() {
        let coalgebra = Arc::new(A0_coalgebra());

        let fp = kComodule::fp_comodule(coalgebra, UniGrading::zero());

        let mut res: Resolution<UniGrading, kComodule<UniGrading, F2, FlatMatrix<F2>>> =
            Resolution::new(fp);

        res.resolve_to_s(20, UniGrading(20));

        let sseq = res.generate_sseq("A0");
        let comp_sseq: SSeq = serde_json::from_str(include_str!("./A(0).json")).unwrap();
        assert_eq!(sseq, comp_sseq);
    }

    #[test]
    fn test_a1_resolution() {
        let input = include_str!("../examples/direct/A(1).txt");
        let coalgebra = Arc::new(kCoalgebra::parse(input, UniGrading::infty()).unwrap().0);

        let fp = kComodule::fp_comodule(coalgebra, UniGrading::zero());

        let mut res: Resolution<UniGrading, kComodule<UniGrading, F2, FlatMatrix<F2>>> =
            Resolution::new(fp);

        res.resolve_to_s(20, UniGrading(20));

        let p = res.generate_sseq("A(1)");
        let comp_p: SSeq = serde_json::from_str(include_str!("./A(1).json")).unwrap();
        assert_eq!(p, comp_p);
    }

    #[test]
    fn test_a2_resolution_direct() {
        let input = include_str!("../examples/direct/A(2).txt");
        let coalgebra = Arc::new(kCoalgebra::parse(input, UniGrading::infty()).unwrap().0);

        let fp = kComodule::fp_comodule(coalgebra, UniGrading::zero());

        let mut res: Resolution<UniGrading, kComodule<UniGrading, F2, FlatMatrix<F2>>> =
            Resolution::new(fp);

        res.resolve_to_s(20, UniGrading(20));

        let p = res.generate_sseq("A(2)");
        let comp_p: SSeq = serde_json::from_str(include_str!("./A(2).json")).unwrap();
        assert_eq!(p, comp_p);
    }

    #[test]
    fn test_a2_resolution_poly() {
        let input = include_str!("../examples/polynomial/A(2).txt");

        let coalgebra = Arc::new(kCoalgebra::parse(input, UniGrading::infty() - UniGrading(10)).unwrap().0);

        let fp = kComodule::fp_comodule(coalgebra, UniGrading::zero());

        let mut res: Resolution<UniGrading, kComodule<UniGrading, F2, FlatMatrix<F2>>> =
            Resolution::new(fp);

        res.resolve_to_s(20, UniGrading(20));
        dbg!(&res);

        let p = res.generate_sseq("A(2)");
        let comp_p: SSeq = serde_json::from_str(include_str!("./A(2).json")).unwrap();
        assert_eq!(p, comp_p);
    }

    #[test]
    fn test_z2_coalgebra_generation() {
        let coalg = Arc::new(Z2::generate_coalgebra::<F2>().unwrap());
        
        let fp = kComodule::fp_comodule(coalg, Z2::zero());

        let mut res = Resolution::new(fp);

        res.resolve_to_s(20, Z2::infty());

        let _p = res.generate_sseq("Z2");
    }
}
