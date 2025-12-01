#[cfg(test)]
mod tests {
    use crate::{basiselement::kBasisElement, grading::UniGrading, linalg::{field::F2, flat_matrix::FlatMatrix, matrix::RModMorphism, ring::{CRing, UniPolRing}}, module::module::{Module, cohomology}};

    #[test]
    fn test_cohom_1() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(9)));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), Some(2)));



        f.set(0, 0, UniPolRing(F2::one(), 6));
        g.set(0, 0, UniPolRing::<F2>::one());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0].2, Some(4));
    }

    #[test]
    fn test_cohom_2() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(9)));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), Some(2)));



        f.set(0, 0, UniPolRing(F2::one(), 6));
        g.set(0, 0, UniPolRing::<F2>::one());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0].2, Some(4));
    }

    #[test]
    fn test_cohom_3() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(1, 2);
        let mut g = FlatMatrix::zero(2, 1);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        n.push((fake_el.clone(), UniGrading(0), Some(2)));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), Some(1)));



        f.set(0, 0, UniPolRing(F2::one(), 0));
        f.set(0, 1, UniPolRing(F2::one(), 0));
        g.set(0, 0, UniPolRing::<F2>::one());
        g.set(1, 0, UniPolRing::<F2>::one());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }
    
    #[test]
    fn test_cohom_4() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(2, 3);
        let mut g = FlatMatrix::zero(3, 2);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        n.push((fake_el.clone(), UniGrading(0), Some(2)));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), Some(1)));
        q.push((fake_el.clone(), UniGrading(0), Some(1)));



        f.set(0, 1, UniPolRing(F2::one(), 0));
        f.set(1, 0, UniPolRing(F2::one(), 0));
        f.set(1, 2, UniPolRing(F2::one(), 0));

        g.set(0, 0, UniPolRing::<F2>::one());
        g.set(2, 0, UniPolRing::<F2>::one());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_5() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(2, 3);
        let mut g = FlatMatrix::zero(3, 2);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        n.push((fake_el.clone(), UniGrading(0), Some(2)));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), Some(1)));
        q.push((fake_el.clone(), UniGrading(0), Some(1)));



        f.set(0, 1, UniPolRing(F2::one(), 0));
        f.set(1, 0, UniPolRing(F2::one(), 0));
        f.set(1, 2, UniPolRing(F2::one(), 0));

        g.set(0, 0, UniPolRing::<F2>::one());
        g.set(2, 0, UniPolRing::<F2>::one());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_6() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(2, 2);
        let mut g = FlatMatrix::zero(2, 1);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        n.push((fake_el.clone(), UniGrading(0), Some(2)));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), Some(1)));



        f.set(0, 0, UniPolRing(F2::one(), 0));
        f.set(1, 0, UniPolRing(F2::one(), 0));
        f.set(0, 1, UniPolRing(F2::one(), 1));
        f.set(1, 1, UniPolRing(F2::one(), 1));

        g.set(1, 0, UniPolRing::<F2>::one());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0].2, Some(1));
    }

    #[test]
    fn test_cohom_7() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(2, 3);
        let mut g = FlatMatrix::zero(3, 3);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        n.push((fake_el.clone(), UniGrading(0), Some(2)));
        n.push((fake_el.clone(), UniGrading(0), Some(3)));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), Some(1)));
        q.push((fake_el.clone(), UniGrading(0), Some(2)));
        q.push((fake_el.clone(), UniGrading(0), Some(3)));

        f.set(0, 0, UniPolRing(F2::zero(), 0));
        f.set(0, 1, UniPolRing(F2::one(), 1));
        f.set(0, 2, UniPolRing(F2::one(), 1));
        f.set(1, 0, UniPolRing(F2::one(), 0));
        f.set(1, 1, UniPolRing(F2::one(), 0));
        f.set(1, 2, UniPolRing(F2::one(), 0));
        
        g.set(0, 0, UniPolRing::<F2>::one());
        g.set(2, 0, UniPolRing::<F2>::one());
        g.set(1, 1, UniPolRing::<F2>::one());
        g.set(2, 1, UniPolRing::<F2>::one());
        g.set(1, 2, UniPolRing(F2::one(), 1));
        g.set(2, 2, UniPolRing(F2::one(), 1));


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_8() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(3, 2);
        let mut g = FlatMatrix::zero(2, 3);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(2)));
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), None));
        q.push((fake_el.clone(), UniGrading(0), Some(2)));
        q.push((fake_el.clone(), UniGrading(0), Some(1)));

        f.set(0, 0, UniPolRing(F2::zero(), 0));
        f.set(1, 0, UniPolRing(F2::one(), 1));
        f.set(2, 0, UniPolRing(F2::one(), 0));
        f.set(2, 1, UniPolRing(F2::one(), 0));
        
        g.set(0, 1, UniPolRing(F2::one(), 1));
        g.set(1, 1, UniPolRing(F2::one(), 1));
        g.set(0, 2, UniPolRing::<F2>::one());
        g.set(1, 2, UniPolRing::<F2>::one());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_9() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(2, 4);
        let mut g = FlatMatrix::zero(4, 7);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        n.push((fake_el.clone(), UniGrading(0), Some(2)));
        n.push((fake_el.clone(), UniGrading(0), Some(3)));
        n.push((fake_el.clone(), UniGrading(0), None));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), Some(1)));
        q.push((fake_el.clone(), UniGrading(0), Some(1)));
        q.push((fake_el.clone(), UniGrading(0), Some(2)));
        q.push((fake_el.clone(), UniGrading(0), Some(2)));
        q.push((fake_el.clone(), UniGrading(0), Some(3)));
        q.push((fake_el.clone(), UniGrading(0), Some(3)));
        q.push((fake_el.clone(), UniGrading(0), Some(4)));

        f.set(0, 1, UniPolRing(F2::one(), 1));
        f.set(0, 2, UniPolRing(F2::one(), 1));
        f.set(1, 0, UniPolRing(F2::one(), 0));
        f.set(1, 1, UniPolRing(F2::one(), 0));
        f.set(1, 2, UniPolRing(F2::one(), 0));
        
        g.set(0, 0, UniPolRing(F2::one(), 0));
        g.set(2, 0, UniPolRing(F2::one(), 0));
        g.set(3, 1, UniPolRing(F2::one(), 0));
        g.set(1, 2, UniPolRing(F2::one(), 0));
        g.set(2, 2, UniPolRing(F2::one(), 0));
        g.set(3, 3, UniPolRing(F2::one(), 0));
        g.set(3, 5, UniPolRing(F2::one(), 0));
        g.set(3, 6, UniPolRing(F2::one(), 0));

        g.set(1, 4, UniPolRing(F2::one(), 1));
        g.set(2, 4, UniPolRing(F2::one(), 1));


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0].2, None);
    }
}