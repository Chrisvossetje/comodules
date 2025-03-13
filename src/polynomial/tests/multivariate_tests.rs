mod tests {
    use crate::{
        linalg::field::Fp,
        polynomial::multivariate::{Monomial, MultivariatePolynomial},
    };

    type TestField = Fp<23>;

    #[test]
    fn test_minimize() {
        let poly = MultivariatePolynomial::new(vec![
            (Monomial([1, 0]), TestField { 0: 2 }),
            (Monomial([1, 0]), TestField { 0: 3 }),
            (Monomial([0, 1]), TestField { 0: 0 }),
        ]);
        assert_eq!(poly.0, vec![(Monomial([1, 0]), TestField { 0: 5 })]);
    }

    #[test]
    fn test_is_minimal() {
        let minimal_poly = MultivariatePolynomial(vec![
            (Monomial([0, 1]), TestField { 0: 1 }),
            (Monomial([1, 0]), TestField { 0: 2 }),
        ]);
        assert!(minimal_poly.is_minimal());

        let non_minimal_poly = MultivariatePolynomial(vec![
            (Monomial([1, 0]), TestField { 0: 2 }),
            (Monomial([0, 1]), TestField { 0: 1 }),
        ]);
        assert!(!non_minimal_poly.is_minimal());

        let duplicate_monomial_poly = MultivariatePolynomial(vec![
            (Monomial([1, 0]), TestField { 0: 2 }),
            (Monomial([1, 0]), TestField { 0: 3 }),
        ]);
        assert!(!duplicate_monomial_poly.is_minimal());
    }

    #[test]
    fn test_addition_minimal() {
        let poly1 = MultivariatePolynomial::new(vec![
            (Monomial([1, 0]), TestField { 0: 2 }),
            (Monomial([0, 1]), TestField { 0: 3 }),
        ]);
        let poly2 = MultivariatePolynomial::new(vec![
            (Monomial([1, 0]), TestField { 0: 4 }),
            (Monomial([0, 1]), TestField { 0: 20 }),
        ]);
        let result = MultivariatePolynomial::new(vec![(Monomial([1, 0]), TestField { 0: 6 })]);
        let sum = poly1 + poly2;
        assert!((sum) == result);
        assert!((sum).is_minimal());
    }

    #[test]
    fn test_multiplication_minimal() {
        let poly1 = MultivariatePolynomial::new(vec![(Monomial([1, 0]), TestField { 0: 2 })]);
        let poly2 = MultivariatePolynomial::new(vec![(Monomial([0, 1]), TestField { 0: 3 })]);
        let result = MultivariatePolynomial::new(vec![(Monomial([1, 1]), TestField { 0: 6 })]);
        let mult = poly1 * poly2;
        assert!((mult) == result);
        assert!((mult).is_minimal());
    }
}
