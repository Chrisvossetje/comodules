use std::{
    cmp::Ordering,
    fmt::{Debug, Display, Formatter},
    hash::Hash,
    iter::Sum,
    ops::{Add, AddAssign, Sub, SubAssign},
    str::FromStr,
};

use deepsize::DeepSizeOf;
use serde::{Deserialize, Serialize};

use crate::grading::grading::{GradedIndexing, Grading, Parse};

#[derive(
    Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord, Default, Serialize, Deserialize, DeepSizeOf,
)]
pub struct BiGrading(pub i16, pub i16);

impl Add for BiGrading {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        BiGrading(self.0 + other.0, self.1 + other.1)
    }
}

impl AddAssign for BiGrading {
    fn add_assign(&mut self, other: Self) {
        self.0 += other.0;
        self.1 += other.1;
    }
}

impl Sub for BiGrading {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        BiGrading(self.0 - other.0, self.1 - other.1)
    }
}

impl SubAssign for BiGrading {
    fn sub_assign(&mut self, other: Self) {
        self.0 -= other.0;
        self.1 -= other.1;
    }
}

impl Display for BiGrading {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.0, self.1)
    }
}

impl Debug for BiGrading {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.0, self.1)
    }
}

impl Parse for BiGrading {
    fn parse(s: &str) -> Result<Self, String> {
        let s = s.trim();
        // Remove parentheses if present
        let s = if s.starts_with('(') && s.ends_with(')') {
            &s[1..s.len() - 1]
        } else {
            s
        };

        let parts: Vec<&str> = s.split(',').collect();
        if parts.len() != 2 {
            return Err(format!("Grade: {} could not be parsed", s));
        }
        let t = i16::from_str(parts[0].trim())
            .map_err(|_| format!("Grade: {} could not be parsed", s))?;
        let s_val = i16::from_str(parts[1].trim())
            .map_err(|_| format!("Grade: {} could not be parsed", s))?;
        Ok(BiGrading(t, s_val))
    }
}

impl Sum for BiGrading {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(BiGrading::zero(), |a, b| a + b)
    }
}

impl<A> GradedIndexing<A, (usize, usize)> for Vec<Vec<A>> {
    fn get(&self, index: (usize, usize)) -> &A {
        &self[index.0][index.1]
    }

    fn mut_get(&mut self, index: (usize, usize)) -> &mut A {
        &mut self[index.0][index.1]
    }
}

impl Grading for BiGrading {
    type ContiguousMemory<A: Send + Sync> = Vec<Vec<A>>;
    type ContiguousIndex = (usize, usize);

    fn to_index(&self) -> Self::ContiguousIndex {
        (self.0 as usize, self.1 as usize)
    }

    fn degree_names() -> Vec<char> {
        vec!['t', 's']
    }

    fn default_formulas() -> (String, String) {
        ("t-s".to_string(), "s".to_string())
    }

    fn export_grade(self) -> Vec<i64> {
        vec![self.0 as i64, self.1 as i64]
    }

    fn zero() -> Self {
        BiGrading(0, 0)
    }

    fn integer_multiplication(self, other: i16) -> Self {
        BiGrading(self.0 * other, self.1 * other)
    }

    fn infty() -> Self {
        BiGrading(i16::MAX, i16::MAX)
    }

    fn iterator_from_zero(&self, include_self: bool) -> Vec<Self> {
        if include_self {
            (0..=self.0)
                .flat_map(|x| (0..=self.1).map(move |y| BiGrading(x, y)))
                .collect()
        } else {
            (0..=self.0)
                .flat_map(|x| (0..=self.1).map(move |y| BiGrading(x, y)))
                .filter(|&BiGrading(x, y)| !(x == self.0 && y == self.1))
                .collect()
        }
    }

    fn init_memory<A: Send + Sync, T: Fn() -> A>(&self, map: T) -> Self::ContiguousMemory<A> {
        let mut v = vec![];
        for _ in 0..=self.0 {
            let mut w = vec![];
            for _ in 0..=self.1 {
                w.push(map());
            }
            v.push(w);
        }
        v
    }

    fn incr(self) -> Self {
        BiGrading(self.0 + 1, self.1 + 1)
    }

    fn compare(self, rhs: &Self) -> Ordering {
        self.cmp(rhs)
    }

    fn nexts(&self) -> Vec<Self> {
        let a = BiGrading(self.0 + 1, self.1);
        let b = BiGrading(self.0, self.1 + 1);
        vec![a, b]
    }

    fn directions() -> usize {
        2
    }

    fn incomings(&self) -> usize {
        let mut m = 0;
        if self.0 > 0 {
            m += 1;
        }
        if self.1 > 0 {
            m += 1;
        }
        m
    }
}
