#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use comodules::{
        comodule::{
            kcoalgebra::{A0_coalgebra, kCoalgebra}, rcoalgebra::{A0_C, tensor_k_coalgebra}, rcomodule::{RCoalgebra, RComodule}, traits::Comodule
        }, grading::{Grading, UniGrading}, linalg::field::F2, resolution::Resolution
    };

    #[test]
    fn test_a0_c_resolution() {
        let coalgebra = Arc::new(A0_C());

        let kt = RComodule::fp_comodule(coalgebra, UniGrading(0));

        let mut res = Resolution::new(kt);

        res.resolve_to_s(20, UniGrading(40));

        let sseq = res.generate_sseq("");

        sseq.save_to_json("LOL.json").unwrap();
    }

    #[test]
    fn test_a0_tensor_resolution() {
        let coalgebra = A0_coalgebra();

        let tensor_coalgebra = Arc::new(tensor_k_coalgebra(coalgebra));

        let kt = RComodule::fp_comodule(tensor_coalgebra, UniGrading(0));

        let mut res = Resolution::new(kt);

        res.resolve_to_s(20, UniGrading(40));

        let sseq = res.generate_sseq("");

        sseq.save_to_json("LOL.json").unwrap();
    }


    #[test]
    fn test_a1_tensor_resolution() {
        let input = include_str!("../examples/direct/A(1).txt");
        let coalgebra = kCoalgebra::parse(input, UniGrading::infty()).unwrap().0;

        let tensor_coalgebra = Arc::new(tensor_k_coalgebra(coalgebra));

        for (gr, m) in &tensor_coalgebra.coaction.maps {
            println!("{gr}");
            println!("{:?}", m);
        }

        let kt = RComodule::fp_comodule(tensor_coalgebra, UniGrading(0));

        let mut res = Resolution::new(kt);

        res.resolve_to_s(20, UniGrading(40));

        let sseq = res.generate_sseq("");

        sseq.save_to_json("LOL.json").unwrap();
    }

    #[test]
    fn test_a1_tensor_parser_resolution() {
        let input = include_str!("../examples/direct/A(1)_dual_graded.txt");
        let coalgebra = RCoalgebra::<UniGrading, F2>::parse(input, UniGrading::infty()).unwrap().0;
        let coalgebra = Arc::new(coalgebra);

        let kt = RComodule::fp_comodule(coalgebra, UniGrading(0));

        let mut res = Resolution::new(kt);

        res.resolve_to_s(20, UniGrading(40));

        let sseq = res.generate_sseq("");

        sseq.save_to_json("LOL.json").unwrap();
    }
    

    #[test]
    fn test_a_tensor_resolution() {
        let input = include_str!("../examples/polynomial/A.txt");
        let coalgebra = kCoalgebra::parse(input, UniGrading(15)).unwrap().0;

        let tensor_coalgebra = Arc::new(tensor_k_coalgebra(coalgebra));
                
        let kt = RComodule::fp_comodule(tensor_coalgebra, UniGrading(0));

        let mut res = Resolution::new(kt);

        res.resolve_to_s(20, UniGrading(15));

        let sseq = res.generate_sseq("");

        sseq.save_to_json("LOL.json").unwrap();
    }


    #[test]
    fn test_a1_c_resolution() {
        let input = include_str!("../examples/direct/A(1)_C.txt");
        let coalgebra = RCoalgebra::<UniGrading, F2>::parse(input, UniGrading::infty()).unwrap().0;
        let coalgebra = Arc::new(coalgebra);

        let kt = RComodule::fp_comodule(coalgebra, UniGrading(0));

        let mut res = Resolution::new(kt);

        res.resolve_to_s(8, UniGrading(40));

        let sseq = res.generate_sseq("");

        sseq.save_to_json("LOL.json").unwrap();
    }
}
