#[cfg(test)]
mod tests {
    use crate::{matrices::flat_matrix::FlatMatrix, rings::{finite_fields::F2, univariate_polynomial_ring::UniPolRing}, unipol::{UniPolModule, cohomology::internal_cohomology}, ring::CRing, matrix::Matrix};

    #[test]
    fn test_cohom_1() {
        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = UniPolModule::new();
        n.push(Some(9));
        
        let mut q = UniPolModule::new();
        q.push(Some(2));

        f.set(0, 0, UniPolRing(F2::one(), 6));
        g.set(0, 0, UniPolRing::<F2>::one());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0], Some(4));
    }

    #[test]
    fn test_cohom_2() {
        let mut f = FlatMatrix::zero(1, 1);
        let mut g = FlatMatrix::zero(1, 1);

        let mut n = UniPolModule::new();
        n.push(Some(9));
        
        let mut q = UniPolModule::new();
        q.push(Some(2));

        f.set(0, 0, UniPolRing(F2::one(), 6));
        g.set(0, 0, UniPolRing::<F2>::one());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0], Some(4));
    }

    #[test]
    fn test_cohom_3() {
        let mut f = FlatMatrix::zero(1, 2);
        let mut g = FlatMatrix::zero(2, 1);

        let mut n = UniPolModule::new();
        n.push(Some(1));
        n.push(Some(2));
        
        let mut q = UniPolModule::new();
        q.push(Some(1));

        f.set(0, 0, UniPolRing(F2::one(), 0));
        f.set(0, 1, UniPolRing(F2::one(), 0));
        g.set(0, 0, UniPolRing::<F2>::one());
        g.set(1, 0, UniPolRing::<F2>::one());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }
    
    #[test]
    fn test_cohom_4() {
        let mut f = FlatMatrix::zero(2, 3);
        let mut g = FlatMatrix::zero(3, 2);

        let mut n = UniPolModule::new();
        n.push(Some(1));
        n.push(Some(1));
        n.push(Some(2));
        
        let mut q = UniPolModule::new();
        q.push(Some(1));
        q.push(Some(1));

        f.set(0, 1, UniPolRing(F2::one(), 0));
        f.set(1, 0, UniPolRing(F2::one(), 0));
        f.set(1, 2, UniPolRing(F2::one(), 0));

        g.set(0, 0, UniPolRing::<F2>::one());
        g.set(2, 0, UniPolRing::<F2>::one());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_5() {
        let mut f = FlatMatrix::zero(2, 3);
        let mut g = FlatMatrix::zero(3, 2);

        let mut n = UniPolModule::new();
        n.push(Some(1));
        n.push(Some(1));
        n.push(Some(2));
        
        let mut q = UniPolModule::new();
        q.push(Some(1));
        q.push(Some(1));

        f.set(0, 1, UniPolRing(F2::one(), 0));
        f.set(1, 0, UniPolRing(F2::one(), 0));
        f.set(1, 2, UniPolRing(F2::one(), 0));

        g.set(0, 0, UniPolRing::<F2>::one());
        g.set(2, 0, UniPolRing::<F2>::one());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_6() {
        let mut f = FlatMatrix::zero(2, 2);
        let mut g = FlatMatrix::zero(2, 1);

        let mut n = UniPolModule::new();
        n.push(Some(1));
        n.push(Some(2));
        
        let mut q = UniPolModule::new();
        q.push(Some(1));

        f.set(0, 0, UniPolRing(F2::one(), 0));
        f.set(1, 0, UniPolRing(F2::one(), 0));
        f.set(0, 1, UniPolRing(F2::one(), 1));
        f.set(1, 1, UniPolRing(F2::one(), 1));

        g.set(1, 0, UniPolRing::<F2>::one());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0], Some(1));
    }

    #[test]
    fn test_cohom_7() {
        let mut f = FlatMatrix::zero(2, 3);
        let mut g = FlatMatrix::zero(3, 3);

        let mut n = UniPolModule::new();
        n.push(Some(1));
        n.push(Some(2));
        n.push(Some(3));
        
        let mut q = UniPolModule::new();
        q.push(Some(1));
        q.push(Some(2));
        q.push(Some(3));

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

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_8() {
        let mut f = FlatMatrix::zero(3, 2);
        let mut g = FlatMatrix::zero(2, 3);

        let mut n = UniPolModule::new();
        n.push(Some(2));
        n.push(Some(1));
        
        let mut q = UniPolModule::new();
        q.push(None);
        q.push(Some(2));
        q.push(Some(1));

        f.set(0, 0, UniPolRing(F2::zero(), 0));
        f.set(1, 0, UniPolRing(F2::one(), 1));
        f.set(2, 0, UniPolRing(F2::one(), 0));
        f.set(2, 1, UniPolRing(F2::one(), 0));
        
        g.set(0, 1, UniPolRing(F2::one(), 1));
        g.set(1, 1, UniPolRing(F2::one(), 1));
        g.set(0, 2, UniPolRing::<F2>::one());
        g.set(1, 2, UniPolRing::<F2>::one());

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 0);
    }

    #[test]
    fn test_cohom_9() {
        let mut f = FlatMatrix::zero(2, 4);
        let mut g = FlatMatrix::zero(4, 7);

        let mut n = UniPolModule::new();
        n.push(Some(1));
        n.push(Some(2));
        n.push(Some(3));
        n.push(None);
        
        let mut q = UniPolModule::new();
        q.push(Some(1));
        q.push(Some(1));
        q.push(Some(2));
        q.push(Some(2));
        q.push(Some(3));
        q.push(Some(3));
        q.push(Some(4));

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

        let (_, cohom) = internal_cohomology(&f, &g, &n, &q);

        assert_eq!(cohom.len(), 1);
        assert_eq!(cohom[0], None);
    }
}