use std::{
    fmt::Debug,
    hash::Hash,
    ops::{Add, AddAssign, Sub, SubAssign},
    str::FromStr,
};

pub trait Grading:
    'static
    + Clone
    + Hash
    + Copy
    + Debug
    + Sized
    + Add<Output = Self>
    + Sub<Output = Self>
    + PartialEq
    + Eq
    + AddAssign
    + SubAssign
    + PartialOrd
    + Ord
    + Sync
    + Send
{
    fn degree_names() -> Vec<char>;
    fn default_formulas() -> (String, String);
    fn export_grade(self) -> Vec<i32>;

    fn incr(self) -> Self;
    fn zero() -> Self;

    fn parse(parse: &str) -> Result<Self, ()>;
}

impl Grading for i32 {
    fn degree_names() -> Vec<char> {
        vec!['t']
    }

    fn default_formulas() -> (String, String) {
        ("t-s".to_string(), "s".to_string())
    }

    fn export_grade(self) -> Vec<i32> {
        vec![self]
    }

    fn zero() -> Self {
        0
    }

    fn incr(self) -> Self {
        self + 1
    }

    fn parse(parse: &str) -> Result<Self, ()> {
        i32::from_str(parse).map_err(|_| ())
    }
}

pub type UniGrading = i32;

#[derive(Debug, Clone, Copy, Hash)]
pub struct BiGrading(i32, i32);
