#[cfg(test)]
mod tests {
    use crate::{basiselement::kBasisElement, grading::UniGrading, linalg::{field::F2, flat_matrix::FlatMatrix, matrix::RModMorphism, ring::{CRing, UniPolRing}}, module::module::{Module, cohomology}};

    #[test]
    fn test_cohom_simple_1() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), None));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), None));



        f.set(0, 0, UniPolRing::<F2>::zero());
        g.set(0, 0, UniPolRing::<F2>::zero());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0].2, None);
    }

    #[test]
    fn test_cohom_simple_2() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), None));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), None));



        f.set(0, 0, UniPolRing::<F2>::zero());
        g.set(0, 0, UniPolRing::<F2>::one());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_simple_3() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), None));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), None));



        f.set(0, 0, UniPolRing::<F2>::one());
        g.set(0, 0, UniPolRing::<F2>::zero());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }


    #[test]
    fn test_cohom_simple_unfree_module_1() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), Some(1)));


        f.set(0, 0, UniPolRing::<F2>::zero());
        g.set(0, 0, UniPolRing::<F2>::zero());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0].2, Some(1));
    }

    #[test]
    fn test_cohom_simple_unfree_module_3() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), Some(1)));



        f.set(0, 0, UniPolRing::<F2>::zero());
        g.set(0, 0, UniPolRing::<F2>::one());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_simple_unfree_module_4() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), Some(1)));



        f.set(0, 0, UniPolRing::<F2>::one());
        g.set(0, 0, UniPolRing::<F2>::zero());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }


    #[test]
    fn test_cohom_larger_n() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::zero(1, 2);
        let mut g = FlatMatrix::zero(2, 1);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), None));
        n.push((fake_el.clone(), UniGrading(0), None));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), None));



        f.set(0, 0, UniPolRing::<F2>::zero());
        g.set(0, 0, UniPolRing::<F2>::zero());


        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 2);
        assert_eq!(cohom[0].2, None);
        assert_eq!(cohom[1].2, None);
    }

    #[test]
    fn test_cohom_zero_n() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::<UniPolRing<F2>>::zero(1, 0);
        let mut g = FlatMatrix::<UniPolRing<F2>>::zero(0, 1);

        let mut n = Module::new();
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), None));

        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_zero_m() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::<UniPolRing<F2>>::zero(0, 1);
        let mut g = FlatMatrix::<UniPolRing<F2>>::zero(1, 1);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), None));
        
        let mut q = Module::new();
        q.push((fake_el.clone(), UniGrading(0), None));

        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0].2, None);
    }

    #[test]
    fn test_cohom_zero_q_1() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::<UniPolRing<F2>>::zero(1, 1);
        let mut g = FlatMatrix::<UniPolRing<F2>>::zero(1, 0);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), None));
        
        let mut q = Module::new();

        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0].2, None);
    }

    #[test]
    fn test_cohom_zero_q_2() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::<UniPolRing<F2>>::zero(1, 1);
        let mut g = FlatMatrix::<UniPolRing<F2>>::zero(1, 0);

        f.set(0, 0, UniPolRing::<F2>::one());


        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), None));
        
        let mut q = Module::new();

        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_zero_q_3() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::<UniPolRing<F2>>::zero(2, 1);
        let mut g = FlatMatrix::<UniPolRing<F2>>::zero(1, 0);

        f.set(1, 0, UniPolRing::<F2>::one());


        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        
        let mut q = Module::new();

        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_zero_q_m() {
        let fake_el = kBasisElement {
            name: "".to_owned(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let mut f = FlatMatrix::<UniPolRing<F2>>::zero(0, 1);
        let mut g = FlatMatrix::<UniPolRing<F2>>::zero(1, 0);

        let mut n = Module::new();
        n.push((fake_el.clone(), UniGrading(0), Some(1)));
        
        let mut q = Module::new();

        let (cohom, _) = cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0].2, Some(1));
    }
}