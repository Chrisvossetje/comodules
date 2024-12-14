#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use comodules::{
        comodule::{
            comodule::Comodule,
            kcomodule::{kComodule, A0_coalgebra},
            kmorphism::kComoduleMorphism,
        },
        linalg::{field::F2, graded::UniGrading},
        resolution::Resolution,
    };

    #[test]
    fn test_a0_resolution() {
        let coalgebra = Arc::new(A0_coalgebra());

        let fp = kComodule::fp_comodule(coalgebra);

        let mut res: Resolution<
            UniGrading,
            kComodule<UniGrading, F2>,
            kComoduleMorphism<UniGrading, F2>,
        > = Resolution::new(fp);

        res.resolve_to_s(4, 10);

        dbg!(res.generate_page());
    }
}
