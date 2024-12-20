use std::sync::Arc;

use std::time::Instant;

use comodules::{
    comodule::{
        kcoalgebra::kCoalgebra, kcomodule::kComodule, kmorphism::kComoduleMorphism,
        traits::Comodule,
    },
    linalg::{field::F2, flat_matrix::FlatMatrix, grading::UniGrading, row_matrix::RowMatrix},
    resolution::Resolution,
};

fn main() {
    let start = Instant::now();

    let input = include_str!("../../examples/kcoalgebras/A(2).txt");
    let coalgebra = Arc::new(kCoalgebra::parse(input).unwrap().0);

    let fp = kComodule::fp_comodule(coalgebra);

    let mut res: Resolution<
        UniGrading,
        kComodule<UniGrading, F2, FlatMatrix<F2>>,
        kComoduleMorphism<UniGrading, F2, FlatMatrix<F2>>,
    > = Resolution::new(fp);

    res.resolve_to_s(80, 140);

    let page = res.generate_page();

    let _ = page.save_to_json("page.json".to_owned());

    println!("\nProgram took: {:.2?}", start.elapsed());
}
