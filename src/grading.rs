use std::{
    cmp::Ordering, fmt::{Debug, Display, Formatter}, hash::Hash, iter::Sum, ops::{Add, AddAssign, Sub, SubAssign}, str::FromStr
};

use serde::{Deserialize, Serialize};

pub trait Parse: Sized {
    fn parse(s: &str) -> Result<Self, String>;
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
{
    fn degree_names() -> Vec<char>;
    fn default_formulas() -> (String, String);
    fn export_grade(self) -> Vec<i32>;

    fn zero() -> Self;

    fn integer_multiplication(self, other: i32) -> Self;

    fn infty() -> Self;
}

pub trait OrderedGrading:
    'static + Grading
{
    fn incr(self) -> Self;
    fn compare(self, rhs: &Self) -> Ordering;
}

pub trait PolyGrading: Grading {
    const TAU_GRADE: Self;
    fn gen_diff_power(source: Self, target: Self) -> Option<usize>;
}

impl PolyGrading for BiGrading {
    const TAU_GRADE: BiGrading = BiGrading(0,-1);


    fn gen_diff_power(source: Self, target: Self) -> Option<usize> {
        let diff = target-source;
        if diff.0 == 0 {
            if diff.1 > 0 {
                None
            } else {
                Some(diff.1.unsigned_abs() as usize) 
            }
        } else {
            None
        }
        
    }
}





#[derive(Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord, Default, Serialize, Deserialize)]
pub struct UniGrading(pub i32);

#[derive(Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord, Default, Serialize, Deserialize)]
pub struct BiGrading(pub i32, pub i32);




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


impl OrderedGrading for UniGrading {
    fn incr(self) -> Self {
        UniGrading(self.0 + 1)
    }
    
    fn compare(self, rhs: &Self) -> Ordering {
        self.0.cmp(&rhs.0)
    }
}

impl Parse for UniGrading {
    fn parse(s: &str) -> Result<Self, String> {
        Ok(UniGrading(i32::from_str(s).map_err(|_| format!("Grade: {} could not be parsed", s))?))
    }

}

impl Grading for UniGrading {
    fn degree_names() -> Vec<char> {
        vec!['t']
    }

    fn default_formulas() -> (String, String) {
        ("t-s".to_string(), "s".to_string())
    }

    fn export_grade(self) -> Vec<i32> {
        vec![self.0]
    }

    fn zero() -> Self {
        UniGrading(0)
    }

    fn integer_multiplication(self, other: i32) -> Self {
        UniGrading(self.0 * other)
    }
    
    fn infty() -> Self {
        UniGrading(i32::MAX)
    }
}













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
            &s[1..s.len()-1]
        } else {
            s
        };
        
        let parts: Vec<&str> = s.split(',').collect();
        if parts.len() != 2 {
            return Err(format!("Grade: {} could not be parsed", s));
        }
        let t = i32::from_str(parts[0].trim())
            .map_err(|_| format!("Grade: {} could not be parsed", s))?;
        let s_val = i32::from_str(parts[1].trim())
            .map_err(|_| format!("Grade: {} could not be parsed", s))?;
        Ok(BiGrading(t, s_val))
    }
}

impl Sum for BiGrading {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(BiGrading::zero(), |a, b| a + b)
    }
}

impl OrderedGrading for BiGrading {
    fn incr(self) -> Self {
        BiGrading(self.0 + 1, self.1 + 1)
    }

    fn compare(self, rhs: &Self) -> Ordering {
        self.cmp(rhs)
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

    fn integer_multiplication(self, other: i32) -> Self {
        BiGrading(self.0 * other, self.1 * other)
    }
    
    fn infty() -> Self {
        BiGrading(i32::MAX, i32::MAX)
    }
}
