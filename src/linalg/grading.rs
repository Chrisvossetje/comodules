use std::{
    fmt::{Debug, Display, Formatter},
    hash::Hash,
    iter::Sum,
    ops::{Add, AddAssign, Sub, SubAssign},
    str::FromStr,
};

pub trait Grading:
    'static
    + Clone
    + Hash
    + Copy
    + Debug
    + Display
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
    + FromStr
    + Sum
{
    fn degree_names() -> Vec<char>;
    fn default_formulas() -> (String, String);
    fn export_grade(self) -> Vec<i32>;

    fn incr(self) -> Self;
    fn zero() -> Self;
    fn infty() -> Self;

    fn integer_multiplication(self, other: i32) -> Self;

    fn parse(parse: &str) -> Result<Self, String>;
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

    fn infty() -> Self {
        i32::MAX
    }

    fn incr(self) -> Self {
        self + 1
    }

    fn parse(parse: &str) -> Result<Self, String> {
        i32::from_str(parse).map_err(|_| format!("Grade: {} could not be parsed", parse))
    }

    fn integer_multiplication(self, other: i32) -> Self {
        self * other
    }
}

pub type UniGrading = i32;

#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct BiGrading(pub i32, pub i32);

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

impl FromStr for BiGrading {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts: Vec<&str> = s
            .trim_matches(|p| p == '(' || p == ')')
            .split(',')
            .collect();
        if parts.len() != 2 {
            return Err(());
        }
        let x = parts[0].trim().parse::<i32>().map_err(|_| ())?;
        let y = parts[1].trim().parse::<i32>().map_err(|_| ())?;
        Ok(BiGrading(x, y))
    }
}

impl Sum for BiGrading {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(BiGrading::zero(), |a, b| a + b)
    }
}

impl Grading for BiGrading {
    fn degree_names() -> Vec<char> {
        vec!['t', 's']
    }

    fn default_formulas() -> (String, String) {
        ("t-s".to_string(), "s".to_string())
    }

    fn export_grade(self) -> Vec<i32> {
        vec![self.0, self.1]
    }

    fn zero() -> Self {
        BiGrading(0, 0)
    }

    fn infty() -> Self {
        BiGrading(i32::MAX, i32::MAX)
    }

    fn incr(self) -> Self {
        BiGrading(self.0 + 1, self.1 + 1)
    }

    fn parse(parse: &str) -> Result<Self, String> {
        let parts: Vec<&str> = parse.split(',').collect();
        if parts.len() != 2 {
            return Err(format!("Grade: {} could not be parsed", parse));
        }
        let t = i32::from_str(parts[0].trim())
            .map_err(|_| format!("Grade: {} could not be parsed", parse))?;
        let s = i32::from_str(parts[1].trim())
            .map_err(|_| format!("Grade: {} could not be parsed", parse))?;
        Ok(BiGrading(t, s))
    }

    fn integer_multiplication(self, other: i32) -> Self {
        BiGrading(self.0 * other, self.1 * other)
    }
}
