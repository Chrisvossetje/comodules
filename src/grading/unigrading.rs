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
pub struct UniGrading(pub i16);

impl Add for UniGrading {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        UniGrading(self.0 + other.0)
    }
}

impl AddAssign for UniGrading {
    fn add_assign(&mut self, other: Self) {
        self.0 += other.0;
    }
}

impl Sub for UniGrading {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        UniGrading(self.0 - other.0)
    }
}

impl SubAssign for UniGrading {
    fn sub_assign(&mut self, other: Self) {
        self.0 -= other.0;
    }
}

impl Sum for UniGrading {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(UniGrading::zero(), |a, b| a + b)
    }
}

impl Display for UniGrading {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Debug for UniGrading {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Parse for UniGrading {
    fn parse(s: &str) -> Result<Self, String> {
        Ok(UniGrading(i16::from_str(s).map_err(|_| {
            format!("Grade: {} could not be parsed", s)
        })?))
    }
}

impl<A> GradedIndexing<A, usize> for Vec<A> {
    fn get(&self, index: usize) -> &A {
        &self[index]
    }

    fn mut_get(&mut self, index: usize) -> &mut A {
        &mut self[index]
    }
}

impl Grading for UniGrading {
    type ContiguousMemory<A: Send + Sync> = Vec<A>;
    type ContiguousIndex = usize;

    fn to_index(&self) -> Self::ContiguousIndex {
        self.0 as usize
    }

    fn degree_names() -> Vec<char> {
        vec!['t']
    }

    fn default_formulas() -> (String, String) {
        ("t-s".to_string(), "s".to_string())
    }

    fn export_grade(self) -> Vec<i64> {
        vec![self.0 as i64]
    }

    fn zero() -> Self {
        UniGrading(0)
    }

    fn integer_multiplication(self, other: i16) -> Self {
        UniGrading(self.0 * other)
    }

    fn infty() -> Self {
        UniGrading(i16::MAX)
    }

    fn iterator_from_zero(&self, include_self: bool) -> Vec<Self> {
        if include_self {
            (0..=self.0).map(|x| UniGrading(x)).collect()
        } else {
            (0..self.0).map(|x| UniGrading(x)).collect()
        }
    }

    fn init_memory<A: Send + Sync, T: Fn() -> A>(&self, map: T) -> Self::ContiguousMemory<A> {
        let mut v = vec![];
        for _ in 0..=self.0 {
            v.push(map());
        }
        v
    }

    fn incr(self) -> Self {
        UniGrading(self.0 + 1)
    }

    fn compare(self, rhs: &Self) -> Ordering {
        self.0.cmp(&rhs.0)
    }
    
    fn nexts(&self) -> Vec<Self> {
        vec![UniGrading(self.0 + 1)]
    }

    fn directions() -> usize {
        1    
    }

    fn incomings(&self) -> usize {
        if self.0 > 0 {
            1
        } else {
            0
        }
    }
}
