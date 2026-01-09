use std::{
    cmp::Ordering,
    fmt::{Debug, Display},
    hash::Hash,
    iter::Sum,
    ops::{Add, AddAssign, Sub, SubAssign},
};

use deepsize::DeepSizeOf;

pub trait Parse: Sized {
    fn parse(s: &str) -> Result<Self, String>;
}

pub trait GradedIndexing<A, I>: Sized {
    fn get(&self, index: I) -> &A;
    fn mut_get(&mut self, index: I) -> &mut A;
}

pub trait Grading:
    'static
    + Clone
    + Hash
    + Copy
    + Debug
    + Display
    + Sized
    + PartialEq
    + Eq
    + Default
    + Add<Output = Self>
    + Sub<Output = Self>
    + AddAssign
    + SubAssign
    + Sync
    + Send
    + Parse
    + Sum
    + PartialOrd
    + Ord
    + DeepSizeOf
{
    type ContiguousIndex: Send + Sync;
    type ContiguousMemory<A: Send + Sync>: GradedIndexing<A, Self::ContiguousIndex>
        + Send
        + Sync
        + Default;

    fn iterator_from_zero(&self, include_self: bool) -> Vec<Self>;
    fn init_memory<A: Send + Sync, T: Fn() -> A>(&self, map: T) -> Self::ContiguousMemory<A>;
    fn nexts(&self) -> Vec<Self>;
    fn directions() -> usize;
    fn incomings(&self) -> usize;

    fn degree_names() -> Vec<char>;
    fn default_formulas() -> (String, String);
    fn export_grade(self) -> Vec<i64>;

    fn zero() -> Self;

    fn integer_multiplication(self, other: i16) -> Self;

    fn infty() -> Self;

    fn to_index(&self) -> Self::ContiguousIndex;

    fn incr(self) -> Self;
    fn compare(self, rhs: &Self) -> Ordering;
}
