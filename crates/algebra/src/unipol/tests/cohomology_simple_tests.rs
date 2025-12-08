#[cfg(test)]
mod tests {
    use crate::{matrices::flat_matrix::FlatMatrix, rings::{finite_fields::F2, univariate_polynomial_ring::UniPolRing}, unipol::{UniPolModule, cohomology::internal_cohomology}, ring::CRing, matrix::Matrix};

    #[test]
    fn test_cohom_simple_1() {
        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = UniPolModule::new();
        n.push(None);
        
        let mut q = UniPolModule::new();
        q.push(None);

        f.set(0, 0, UniPolRing::<F2>::zero());
        g.set(0, 0, UniPolRing::<F2>::zero());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0], None);
    }

    #[test]
    fn test_cohom_simple_2() {
        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = UniPolModule::new();
        n.push(None);
        
        let mut q = UniPolModule::new();
        q.push(None);

        f.set(0, 0, UniPolRing::<F2>::zero());
        g.set(0, 0, UniPolRing::<F2>::one());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_simple_3() {
        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = UniPolModule::new();
        n.push(None);
        
        let mut q = UniPolModule::new();
        q.push(None);

        f.set(0, 0, UniPolRing::<F2>::one());
        g.set(0, 0, UniPolRing::<F2>::zero());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_simple_unfree_module_1() {
        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = UniPolModule::new();
        n.push(Some(1));
        
        let mut q = UniPolModule::new();
        q.push(Some(1));

        f.set(0, 0, UniPolRing::<F2>::zero());
        g.set(0, 0, UniPolRing::<F2>::zero());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0], Some(1));
    }

    #[test]
    fn test_cohom_simple_unfree_module_3() {
        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = UniPolModule::new();
        n.push(Some(1));
        
        let mut q = UniPolModule::new();
        q.push(Some(1));

        f.set(0, 0, UniPolRing::<F2>::zero());
        g.set(0, 0, UniPolRing::<F2>::one());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_simple_unfree_module_4() {
        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = UniPolModule::new();
        n.push(Some(1));
        
        let mut q = UniPolModule::new();
        q.push(Some(1));

        f.set(0, 0, UniPolRing::<F2>::one());
        g.set(0, 0, UniPolRing::<F2>::zero());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_larger_n() {
        let mut f = FlatMatrix::zero(1, 2);
        let mut g = FlatMatrix::zero(2, 1);

        let mut n = UniPolModule::new();
        n.push(None);
        n.push(None);
        
        let mut q = UniPolModule::new();
        q.push(None);

        f.set(0, 0, UniPolRing::<F2>::zero());
        g.set(0, 0, UniPolRing::<F2>::zero());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 2);
        assert_eq!(cohom[0], None);
        assert_eq!(cohom[1], None);
    }

    #[test]
    fn test_cohom_zero_n() {
        let f = FlatMatrix::<UniPolRing<F2>>::zero(1, 0);
        let g = FlatMatrix::<UniPolRing<F2>>::zero(0, 1);

        let n = UniPolModule::new();
        
        let mut q = UniPolModule::new();
        q.push(None);

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_zero_m() {
        let f = FlatMatrix::<UniPolRing<F2>>::zero(0, 1);
        let g = FlatMatrix::<UniPolRing<F2>>::zero(1, 1);

        let mut n = UniPolModule::new();
        n.push(None);
        
        let mut q = UniPolModule::new();
        q.push(None);

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0], None);
    }

    #[test]
    fn test_cohom_zero_q_1() {
        let f = FlatMatrix::<UniPolRing<F2>>::zero(1, 1);
        let g = FlatMatrix::<UniPolRing<F2>>::zero(1, 0);

        let mut n = UniPolModule::new();
        n.push(None);
        
        let q = UniPolModule::new();

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0], None);
    }

    #[test]
    fn test_cohom_zero_q_2() {
        let mut f = FlatMatrix::<UniPolRing<F2>>::zero(1, 1);
        let g = FlatMatrix::<UniPolRing<F2>>::zero(1, 0);

        f.set(0, 0, UniPolRing::<F2>::one());

        let mut n = UniPolModule::new();
        n.push(None);
        
        let q = UniPolModule::new();

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_zero_q_3() {
        let mut f = FlatMatrix::<UniPolRing<F2>>::zero(2, 1);
        let g = FlatMatrix::<UniPolRing<F2>>::zero(1, 0);

        f.set(1, 0, UniPolRing::<F2>::one());

        let mut n = UniPolModule::new();
        n.push(Some(1));
        
        let q = UniPolModule::new();

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_zero_q_m() {
        let f = FlatMatrix::<UniPolRing<F2>>::zero(0, 1);
        let g = FlatMatrix::<UniPolRing<F2>>::zero(1, 0);

        let mut n = UniPolModule::new();
        n.push(Some(1));
        
        let q = UniPolModule::new();

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0], Some(1));
    }
}