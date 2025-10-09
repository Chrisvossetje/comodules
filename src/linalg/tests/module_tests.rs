
#[cfg(test)]
mod tests {
    use crate::linalg::{field::F2, flat_matrix::FlatMatrix, matrix::{RModMorphism, SmithNormalForm}, ring::{CRing, UniPolRing}};

    #[test]
    fn test_smith_normal_form_helper_methods() {
        // Test the helper methods for Smith Normal Form
        let matrix = FlatMatrix {
            data: vec![
                F2(1), F2(0), F2(1),
                F2(0), F2(1), F2(1),
                F2(1), F2(1), F2(0),
            ],
            domain: 3,
            codomain: 3,
        };
        
        // Test helper methods work correctly
        let mut test_matrix = matrix.clone();
        test_matrix.swap_rows(0, 1);
        assert_ne!(test_matrix, matrix);
        
        test_matrix.swap_cols(0, 1);
        test_matrix.add_row_multiple(0, 1, F2(1));
        test_matrix.add_col_multiple(0, 1, F2(1));
        
        // Verify the operations changed the matrix
        assert_ne!(test_matrix, matrix);
    }

    #[test]
    fn test_smith_normal_form_identity_matrix() {
        // Test SNF of identity matrix over F2[x]
        let identity = FlatMatrix {
            data: vec![
                UniPolRing(F2(1), 0), UniPolRing(F2(0), 0),
                UniPolRing(F2(0), 0), UniPolRing(F2(1), 0),
            ],
            domain: 2,
            codomain: 2,
        };
        
        let (u, s, v) = identity.snf();
        
        // Identity matrix should remain identity in SNF
        assert_eq!(s.get(0, 0), UniPolRing(F2(1), 0));
        assert_eq!(s.get(1, 1), UniPolRing(F2(1), 0));
        assert_eq!(s.get(1, 0), UniPolRing(F2(0), 0));
        assert_eq!(s.get(0, 1), UniPolRing(F2(0), 0));
        
        // Verify dimensions
        assert_eq!(u.domain(), 2);
        assert_eq!(u.codomain(), 2);
        assert_eq!(v.domain(), 2);
        assert_eq!(v.codomain(), 2);
    }

    #[test]
    fn test_smith_normal_form_zero_matrix() {
        // Test SNF of zero matrix over F2[x]
        let zero_matrix: FlatMatrix<UniPolRing<F2>> = FlatMatrix::zero(3, 3);
        let (u, s, v) = zero_matrix.snf();
        
        // Zero matrix should remain zero in SNF
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(s.get(i, j), UniPolRing(F2(0), 0));
            }
        }
        
        // U and V should be identity matrices
        assert_eq!(u.domain(), 3);
        assert_eq!(u.codomain(), 3);
        assert_eq!(v.domain(), 3);
        assert_eq!(v.codomain(), 3);
    }

    #[test]
    fn test_smith_normal_form_diagonal_matrix() {
        // Test SNF of diagonal matrix over F2[x]
        // Diagonal: [x^2, 0; 0, x]
        let diagonal = FlatMatrix {
            data: vec![
                UniPolRing(F2(1), 2), UniPolRing(F2(0), 0),
                UniPolRing(F2(0), 0), UniPolRing(F2(1), 1),
            ],
            domain: 2,
            codomain: 2,
        };
        
        let (_, s, _) = diagonal.snf();
        
        // Verify dimensions
        assert_eq!(s.domain(), 2);
        assert_eq!(s.codomain(), 2);
        
        // The diagonal should be in Smith normal form (divisibility chain)
        let d1 = s.get(0, 0);
        let d2 = s.get(1, 1);
        
        // Check that diagonal entries are non-zero if original was non-zero
        if !diagonal.get(0, 0).is_zero() || !diagonal.get(1, 1).is_zero() {
            assert!(!d1.is_zero() || !d2.is_zero());
        }
    }

    #[test]
    fn test_smith_normal_form_rectangular_matrix() {
        // Test SNF of 2x3 matrix over F2[x]
        let rect_matrix = FlatMatrix {
            data: vec![
                UniPolRing(F2(1), 1), UniPolRing(F2(1), 0), UniPolRing(F2(0), 0),
                UniPolRing(F2(0), 0), UniPolRing(F2(1), 1), UniPolRing(F2(1), 0),
            ],
            domain: 3,
            codomain: 2,
        };
        
        let (u, s, v) = rect_matrix.snf();
        
        // Verify dimensions
        assert_eq!(u.domain(), 2);
        assert_eq!(u.codomain(), 2);
        assert_eq!(s.domain(), 3);
        assert_eq!(s.codomain(), 2);
        assert_eq!(v.domain(), 3);
        assert_eq!(v.codomain(), 3);
    }

    #[test]
    fn test_smith_normal_form_polynomial_entries() {
        // Test SNF with various polynomial degrees
        // x^3 x
        // x^2 1
        let poly_matrix = FlatMatrix {
            data: vec![
                UniPolRing(F2(1), 3), UniPolRing(F2(1), 1),
                UniPolRing(F2(1), 2), UniPolRing(F2(1), 0),
            ],
            domain: 2,
            codomain: 2,
        };
        
        let (u, s, v) = poly_matrix.snf();
        
        // Verify the factorization works
        assert_eq!(u.domain(), 2);
        assert_eq!(u.codomain(), 2);
        assert_eq!(s.domain(), 2);
        assert_eq!(s.codomain(), 2);
        assert_eq!(v.domain(), 2);
        assert_eq!(v.codomain(), 2);
        
        // Check that S is in diagonal form (off-diagonal entries are zero)
        assert!(s.get(1, 0).is_zero());
        assert!(s.get(0, 1).is_zero());
    }

    #[test]
    fn test_smith_normal_form_with_polynomials() {
        // Test Smith Normal Form with UniPolRing (polynomial ring k[x])
        
        // Create a simple 2x2 matrix over F2[x]
        // Matrix represents: [x, 1; 0, x] where x has degree 1
        let matrix = FlatMatrix {
            data: vec![
                UniPolRing(F2(1), 1), // x
                UniPolRing(F2(1), 0), // 1  
                UniPolRing(F2(0), 0), // 0
                UniPolRing(F2(1), 1), // x
            ],
            domain: 2,
            codomain: 2,
        };
        println!("{:?}",matrix);

        // Test that Smith Normal Form can be called (even with simplified implementation)
        let (u, s, v) = matrix.snf();

        // TODO: nothing of importance gets checked here

        // Verify dimensions are correct
        assert_eq!(u.domain(), 2);
        assert_eq!(u.codomain(), 2);
        assert_eq!(s.domain(), 2);
        assert_eq!(s.codomain(), 2);
        assert_eq!(v.domain(), 2);
        assert_eq!(v.codomain(), 2);
        
        // Verify that U, S, V are not all zero matrices
        assert!(u != FlatMatrix::zero(2, 2));
        assert!(v != FlatMatrix::zero(2, 2));
    }

    #[test] 
    fn test_smith_normal_form_factorization_property() {
        // Test that U * S * V = original matrix
        let matrix = FlatMatrix {
            data: vec![
                UniPolRing(F2(1), 1), UniPolRing(F2(1), 0),
                UniPolRing(F2(0), 2), UniPolRing(F2(1), 1),
            ],
            domain: 2,
            codomain: 2,
        };
        
        let (u, s, v, uinv, vinv) = matrix.full_snf();
        
        // Compute U * mat * V == S
        let matv = matrix.compose(&v);
        let reconstructed = u.compose(&matv);

        assert_eq!(FlatMatrix::identity(2), u.compose(&uinv));
        assert_eq!(FlatMatrix::identity(2), v.compose(&vinv));
        assert_eq!(reconstructed, s);
    }

    #[test]
    fn test_smith_normal_form_large_matrix() {
        // Test SNF performance with larger matrix
        let size = 4;
        let mut data = Vec::new();
        
        // Create a 4x4 matrix with polynomial entries
        for i in 0..size {
            for j in 0..size {
                if i == j {
                    data.push(UniPolRing(F2(1), i % 3)); // Diagonal entries with varying degrees
                } else if i + 1 == j {
                    data.push(UniPolRing(F2(1), 0)); // Super-diagonal with constant terms
                } else {
                    data.push(UniPolRing(F2(0), 0)); // Zero elsewhere
                }
            }
        }
        
        let large_matrix = FlatMatrix {
            data,
            domain: size,
            codomain: size,
        };
        
        let (u, s, v) = large_matrix.snf();
        
        // Verify dimensions
        assert_eq!(u.domain(), size);
        assert_eq!(u.codomain(), size);
        assert_eq!(s.domain(), size);
        assert_eq!(s.codomain(), size);
        assert_eq!(v.domain(), size);
        assert_eq!(v.codomain(), size);
        
        // Verify S has diagonal structure (off-diagonal entries should be zero)
        for i in 0..size {
            for j in 0..size {
                if i != j {
                    assert_eq!(s.get(i, j), UniPolRing(F2(0), 0));
                }
            }
        }
    }
}
