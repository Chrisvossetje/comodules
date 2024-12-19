use std::sync::Arc;

use comodules::{
    comodule::{
        comodule::Comodule, kcoalgebra::kCoalgebra, kcomodule::kComodule,
        kmorphism::kComoduleMorphism,
    },
    linalg::{
        field::F2,
        graded::UniGrading,
    },
    resolution::Resolution,
};

fn main() {
    let input = include_str!("../../examples/kcoalgebras/A(2).txt");
    let coalgebra = Arc::new(kCoalgebra::parse(input).unwrap().0);

    let fp = kComodule::fp_comodule(coalgebra);

    let mut res: Resolution<
        UniGrading,
        kComodule<UniGrading, F2>,
        kComoduleMorphism<UniGrading, F2>,
    > = Resolution::new(fp);

    res.resolve_to_s(8, 60);

    let page = res.generate_page();

    let _ = page.save_to_json("page.json".to_owned());
}
