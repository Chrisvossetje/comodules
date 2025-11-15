#[cfg(test)]
mod tests {
    use crate::grading::{BiGrading, Grading, OrderedGrading, Parse};

    #[test]
    fn test_add() {
        let a = BiGrading(1, 2);
        let b = BiGrading(3, 4);
        let result = a + b;
        assert_eq!(result, BiGrading(4, 6));
    }

    #[test]
    fn test_add_assign() {
        let mut a = BiGrading(1, 2);
        let b = BiGrading(3, 4);
        a += b;
        assert_eq!(a, BiGrading(4, 6));
    }

    #[test]
    fn test_sub() {
        let a = BiGrading(5, 7);
        let b = BiGrading(3, 4);
        let result = a - b;
        assert_eq!(result, BiGrading(2, 3));
    }

    #[test]
    fn test_sub_assign() {
        let mut a = BiGrading(5, 7);
        let b = BiGrading(3, 4);
        a -= b;
        assert_eq!(a, BiGrading(2, 3));
    }

    #[test]
    fn test_display() {
        let a = BiGrading(1, 2);
        assert_eq!(format!("{}", a), "(1, 2)");
    }

    #[test]
    fn test_from_str() {
        let a: BiGrading = BiGrading::parse("(1, 2)").unwrap();
        assert_eq!(a, BiGrading(1, 2));
    }

    #[test]
    fn test_sum() {
        let a = BiGrading(1, 2);
        let b = BiGrading(3, 4);
        let c = BiGrading(5, 6);
        let vec = vec![a, b, c];
        let result: BiGrading = vec.into_iter().sum();
        assert_eq!(result, BiGrading(9, 12));
    }

    #[test]
    fn test_ord() {
        let a = BiGrading(1, 2);
        let b = BiGrading(2, 2);
        let c = BiGrading(1, 3);
        assert!(a < b);
        assert!(a < c);
        assert!(c < b);
    }

    #[test]
    fn test_grading_trait() {
        let a = BiGrading(1, 2);
        assert_eq!(BiGrading::degree_names(), vec!['t', 's']);
        assert_eq!(
            BiGrading::default_formulas(),
            ("t-s".to_string(), "s".to_string())
        );
        assert_eq!(a.export_grade(), vec![1, 2]);
        assert_eq!(a.incr(), BiGrading(2, 3));
        assert_eq!(BiGrading::zero(), BiGrading(0, 0));
        assert_eq!(BiGrading::infty(), BiGrading(i32::MAX, i32::MAX));
    }
}
