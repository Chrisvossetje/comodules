use comodules::linalg::field::{Field, Fp};

fn main() {
    let f = Fp::<5>(4);
    dbg!(f);
    dbg!(f.inv().unwrap());
}