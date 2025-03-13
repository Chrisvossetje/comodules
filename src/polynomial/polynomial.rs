use crate::linalg::field::Field;
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

pub trait Polynomial<F: Field>:
    Clone
    + Debug
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
    fn eval(&self, scalar: F) -> F;
    fn scalar_mult(&mut self, scalar: F);
    fn is_zero(&self) -> bool;
    fn one() -> Self;
    fn zero() -> Self;

    fn parse(input: &str) -> Result<Self, String>;
}
