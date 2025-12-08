#[cfg(test)]
mod tests {
    use crate::{field::Field, ring::CRing, rings::finite_fields::Fp};
    
    type Fp23 = Fp<23>;

    #[test]
    fn test_fp23_addition() {
        let a = Fp23 { 0: 10 };
        let b = Fp23 { 0: 15 };
        assert_eq!(a + b, Fp23 { 0: 2 }); // 10 + 15 = 25 mod 23 = 2
    }

    #[test]
    fn test_fp23_subtraction() {
        let a = Fp23 { 0: 10 };
        let b = Fp23 { 0: 15 };
        assert_eq!(a - b, Fp23 { 0: 18 }); // 10 - 15 mod 23 = -5 mod 23 = 18
    }

    #[test]
    fn test_fp23_multiplication() {
        let a = Fp23 { 0: 7 };
        let b = Fp23 { 0: 6 };
        assert_eq!(a * b, Fp23 { 0: 19 }); // 7 * 6 = 42 mod 23 = 19
    }

    #[test]
    fn test_fp23_negation() {
        let a = Fp23 { 0: 5 };
        assert_eq!(-a, Fp23 { 0: 18 }); // -5 mod 23 = 23 - 5 = 18
    }

    #[test]
    fn test_fp23_add_assign() {
        let mut a = Fp23 { 0: 20 };
        a += Fp23 { 0: 5 };
        assert_eq!(a, Fp23 { 0: 2 }); // 20 + 5 = 25 mod 23 = 2
    }

    #[test]
    fn test_fp23_mul_assign() {
        let mut a = Fp23 { 0: 4 };
        a *= Fp23 { 0: 6 };
        assert_eq!(a, Fp23 { 0: 1 }); // 4 * 6 = 24 mod 23 = 1
    }

    #[test]
    fn test_fp23_sub_assign() {
        let mut a = Fp23 { 0: 3 };
        a -= Fp23 { 0: 10 };
        assert_eq!(a, Fp23 { 0: 16 }); // 3 - 10 mod 23 = -7 mod 23 = 16
    }

    #[test]
    fn test_fp23_inversion() {
        let a = Fp23 { 0: 3 };
        assert_eq!(a.inv(), Some(Fp23 { 0: 8 })); // 3 * 8 mod 23 = 1
    }

    #[test]
    fn test_fp23_inversion_none_for_zero() {
        let a = Fp23 { 0: 0 };
        assert_eq!(a.inv(), None); // 0 has no multiplicative inverse
    }

    #[test]
    fn test_fp23_one_and_zero() {
        assert_eq!(Fp23::one(), Fp23 { 0: 1 });
        assert_eq!(Fp23::zero(), Fp23 { 0: 0 });
    }

    #[test]
    fn test_fp23_is_zero() {
        assert!(Fp23 { 0: 0 }.is_zero());
        assert!(!Fp23 { 0: 5 }.is_zero());
    }

    #[test]
    fn test_fp23_characteristic() {
        assert_eq!(Fp23::get_characteristic(), 23);
    }

    #[test]
    fn test_fp23_sum() {
        let values = vec![Fp23 { 0: 3 }, Fp23 { 0: 5 }, Fp23 { 0: 17 }];
        let sum: Fp23 = values.into_iter().sum();
        assert_eq!(sum, Fp23 { 0: 2 }); // (3 + 5 + 17) mod 23 = 25 mod 23 = 2
    }
}
