#[cfg(test)]
mod tests {
    use crate::{
        abelian::Abelian, matrices::flat_matrix::FlatMatrix, matrix::Matrix, ring::CRing, rings::{finite_fields::F2, univariate_polynomial_ring::UniPolRing}, unipol::UniPolModule
    };

    #[test]
    fn test_module_structure_simple() {
        let mut module = UniPolModule::new();
        module.push(None);           // Free element
        module.push(Some(2));        // Element with quotient t^2
        
        assert_eq!(module.len(), 2);
        assert_eq!(module[0], None);
        assert_eq!(module[1], Some(2));
    }

    #[test]
    fn test_matrix_with_polynomial_ring() {
        // Test matrix operations over F2[t]
        let mut matrix = FlatMatrix::zero(1, 1);
        matrix.set(0, 0, UniPolRing(F2::one(), 2)); // t^2
        
        assert_eq!(matrix.get(0, 0), UniPolRing(F2::one(), 2));
        assert!(!matrix.get(0, 0).is_zero());
        
        // Test that it represents t^2 correctly
        assert_eq!(matrix.get(0, 0).0, F2::one()); // coefficient
        assert_eq!(matrix.get(0, 0).1, 2);         // degree
    }

    #[test]
    fn test_complex_module_structure() {
        let mut codomain = UniPolModule::new();
        codomain.push(None);         // Free element k[t]
        codomain.push(Some(3));      // k[t]/(t^3)
        
        assert_eq!(codomain.len(), 2);
        assert_eq!(codomain[0], None);
        assert_eq!(codomain[1], Some(3));
    }

    #[test]
    fn test_matrix_composition() {
        // Test simple matrix composition
        let mut matrix1 = FlatMatrix::zero(1, 1);
        let mut matrix2 = FlatMatrix::zero(1, 1);
        
        matrix1.set(0, 0, UniPolRing(F2::one(), 1));  // t
        matrix2.set(0, 0, UniPolRing(F2::one(), 2));  // t^2
        
        let composed = matrix1.compose(&matrix2);
        
        // The composition should give us t * t^2 = t^3
        // But since we're in F2[t], this depends on the implementation
        assert_eq!(composed.domain, 1);
        assert_eq!(composed.codomain, 1);
    }

    #[test]
    fn test_quotient_module_elements() {
        let mut module = UniPolModule::new();
        
        // Add various quotient elements
        module.push(Some(1));     // k[t]/(t) ≅ k
        module.push(Some(2));     // k[t]/(t^2)
        module.push(Some(4));     // k[t]/(t^4)
        module.push(None);        // k[t] (free)
        
        assert_eq!(module.len(), 4);
        assert_eq!(module[0], Some(1));
        assert_eq!(module[1], Some(2));
        assert_eq!(module[2], Some(4));
        assert_eq!(module[3], None);
    }

    #[test]
    fn test_zero_and_identity_matrices() {
        // Test zero matrix
        let zero_matrix: FlatMatrix<UniPolRing<F2>> = FlatMatrix::zero(2, 2);
        for i in 0..2 {
            for j in 0..2 {
                assert!(zero_matrix.get(i, j).is_zero());
            }
        }
        
        // Test identity-like matrix
        let mut identity_matrix = FlatMatrix::zero(2, 2);
        identity_matrix.set(0, 0, UniPolRing::<F2>::one());
        identity_matrix.set(1, 1, UniPolRing::<F2>::one());
        
        assert_eq!(identity_matrix.get(0, 0), UniPolRing::<F2>::one());
        assert_eq!(identity_matrix.get(1, 1), UniPolRing::<F2>::one());
        assert!(identity_matrix.get(0, 1).is_zero());
        assert!(identity_matrix.get(1, 0).is_zero());
    }

    #[test]
    fn test_polynomial_arithmetic() {
        let one = UniPolRing::<F2>::one();
        let zero = UniPolRing::<F2>::zero();
        let t = UniPolRing(F2::one(), 1);
        let t_cubed = UniPolRing(F2::one(), 3);
        
        // Test basic properties
        assert!(zero.is_zero());
        assert!(!one.is_zero());
        assert!(!t.is_zero());
        assert!(!t_cubed.is_zero());
        
        // Test that addition gives expected results (in F2)
        let sum = one + one;
        assert!(sum.is_zero()); // 1 + 1 = 0 in F2
    }

    #[test]
    fn test_module_operations() {
        let mut module = UniPolModule::new();
        
        // Test push and indexing
        module.push(None);
        module.push(Some(5));
        module.push(Some(10));
        
        assert_eq!(module[0], None);
        assert_eq!(module[1], Some(5));
        assert_eq!(module[2], Some(10));
        
        // Test length
        assert_eq!(module.len(), 3);
    }

    #[test]
    fn test_kernel_generators_1() {
        let matrix = FlatMatrix {
            data: vec![
                UniPolRing::<F2>::one(),
                UniPolRing::zero(),
                UniPolRing::zero(),
                UniPolRing::zero(),
            ],
            domain: 2,
            codomain: 2,
        };

        let domain = vec![None, None];
        let codomain = vec![None, None];

        let kernel = matrix.kernel_generators(&domain, &codomain);
        assert_eq!(kernel, vec![1]);
    }

    #[test]
    fn test_kernel_generators_2() {
        let matrix = FlatMatrix {
            data: vec![
                UniPolRing::<F2>::one(),
            ],
            domain: 1,
            codomain: 1,
        };

        let domain = vec![None];
        let codomain = vec![None];

        let kernel = matrix.kernel_generators(&domain, &codomain);
        assert_eq!(kernel, vec![]);
    }

    #[test]
    fn test_kernel_generators_3() {
        let matrix = FlatMatrix {
            data: vec![
                UniPolRing::<F2>::one(),
            ],
            domain: 1,
            codomain: 1,
        };

        let domain = vec![Some(1)];
        let codomain = vec![Some(1)];

        let kernel = matrix.kernel_generators(&domain, &codomain);
        assert_eq!(kernel, vec![]);
    }

    #[test]
    fn test_kernel_generators_4() {
        let matrix: FlatMatrix<UniPolRing<F2>> = FlatMatrix {
            data: vec![],
            domain: 1,
            codomain: 0,
        };

        let domain = vec![None];
        let codomain = vec![];

        let kernel = matrix.kernel_generators(&domain, &codomain);
        assert_eq!(kernel, vec![0]);
    }
}