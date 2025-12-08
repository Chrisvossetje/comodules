use std::ops::{Add, Sub, Neg, Mul, AddAssign, SubAssign, MulAssign};
use std::fmt::Debug;

pub trait CRing:
    Clone
    + Copy
    + Debug
    + Default
    + Sized
    + PartialEq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Neg<Output = Self>
    + Mul<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + std::iter::Sum
    + Sync
    + Send
{
    fn is_zero(&self) -> bool;
    fn one() -> Self;
    fn zero() -> Self {
        Self::default()
    }

    fn parse(input: &str) -> Result<Self, String>;

    fn dot_product(l: &Vec<Self>, r: &Vec<Self>) -> Self {
        l.iter().zip(r.iter()).map(|(x, y)| *x * *y).sum()
    }


    fn is_unit(&self) -> bool;

    // 
    fn try_inverse(&self) -> Option<Self>;

    // If you know an element is a unit
    fn unsafe_inverse(&self) -> Self;
}


pub(crate) trait ValuationRing: CRing {
    /// true iff self | other
    /// if false then other | self
    fn divides(self, other: &Self) -> bool;
    fn unsafe_divide(self, div: Self) -> Self;
}

// TODO : Do potentially more SNF stuff ? Prob remove
// pub(crate) trait PID: CRing {
//     fn gcd(&self, other: &Self) -> Self;
// }

