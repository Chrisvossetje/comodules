#[cfg(test)]
mod tests {
    use crate::{
        matrices::flat_matrix::FlatMatrix,
        rings::{finite_fields::F2, univariate_polynomial_ring::UniPolRing},
        unipol::UniPolModule,
        matrix::Matrix,
        ring::CRing,
    };

    #[test]
    fn test_module_creation() {
        let mut module = UniPolModule::new();
        
        // Add some elements with quotient relations
        module.push(None);           // Free element
        module.push(Some(3));        // Element with quotient t^3
        module.push(Some(2));        // Element with quotient t^2
        
        assert_eq!(module.len(), 3);
        assert_eq!(module[0], None);
        assert_eq!(module[1], Some(3));
        assert_eq!(module[2], Some(2));
    }

    #[test]
    fn test_zero_matrix_operations() {
        let zero_matrix: FlatMatrix<UniPolRing<F2>> = FlatMatrix::zero(2, 3);
        
        // Verify dimensions
        assert_eq!(zero_matrix.domain, 2);
        assert_eq!(zero_matrix.codomain, 3);
        
        // Verify all entries are zero
        for i in 0..2 {
            for j in 0..3 {
                assert!(zero_matrix.get(i, j).is_zero());
            }
        }
    }

    #[test]
    fn test_basic_matrix_operations() {
        let mut matrix = FlatMatrix::zero(2, 2);
        
        // Set some entries
        matrix.set(0, 0, UniPolRing(F2::one(), 0));  // 1
        matrix.set(0, 1, UniPolRing(F2::one(), 1));  // t
        matrix.set(1, 0, UniPolRing(F2::one(), 2));  // t^2
        matrix.set(1, 1, UniPolRing::zero());        // 0
        
        // Verify entries
        assert_eq!(matrix.get(0, 0), UniPolRing(F2::one(), 0));
        assert_eq!(matrix.get(0, 1), UniPolRing(F2::one(), 1));
        assert_eq!(matrix.get(1, 0), UniPolRing(F2::one(), 2));
        assert!(matrix.get(1, 1).is_zero());
    }

    #[test]
    fn test_empty_module() {
        let module = UniPolModule::new();
        assert_eq!(module.len(), 0);
        assert!(module.is_empty());
    }

    #[test]
    fn test_module_with_various_quotients() {
        let mut module = UniPolModule::new();
        
        // Add elements with different quotient relations
        module.push(None);        // Free element k[t]
        module.push(Some(1));     // k[t]/(t)
        module.push(Some(4));     // k[t]/(t^4)
        module.push(None);        // Another free element
        module.push(Some(10));    // k[t]/(t^10)
        
        assert_eq!(module.len(), 5);
        
        // Verify quotient relations
        assert_eq!(module[0], None);
        assert_eq!(module[1], Some(1));
        assert_eq!(module[2], Some(4));
        assert_eq!(module[3], None);
        assert_eq!(module[4], Some(10));
    }

    #[test]
    fn test_polynomial_ring_operations() {
        let zero = UniPolRing::<F2>::zero();
        let one = UniPolRing::<F2>::one();
        let t = UniPolRing(F2::one(), 1);
        let t_squared = UniPolRing(F2::one(), 2);
        
        // Test basic properties
        assert!(zero.is_zero());
        assert!(!one.is_zero());
        assert!(!t.is_zero());
        assert!(!t_squared.is_zero());
        
        // Test degree
        assert_eq!(one.1, 0);  // degree of 1
        assert_eq!(t.1, 1);    // degree of t
        assert_eq!(t_squared.1, 2);  // degree of t^2
    }

    #[test]
    fn test_matrix_polynomial_entries() {
        let mut matrix = FlatMatrix::zero(1, 3);
        
        // Set polynomial entries
        matrix.set(0, 0, UniPolRing(F2::one(), 0));   // 1
        matrix.set(0, 1, UniPolRing(F2::one(), 1));   // t
        matrix.set(0, 2, UniPolRing(F2::one(), 3));   // t^3
        
        // Verify polynomial degrees
        assert_eq!(matrix.get(0, 0).1, 0);
        assert_eq!(matrix.get(0, 1).1, 1);
        assert_eq!(matrix.get(0, 2).1, 3);
        
        // Verify coefficients
        assert_eq!(matrix.get(0, 0).0, F2::one());
        assert_eq!(matrix.get(0, 1).0, F2::one());
        assert_eq!(matrix.get(0, 2).0, F2::one());
    }
}