use std::ops::{Add, Sub, Neg, Mul, AddAssign, SubAssign, MulAssign};
use std::fmt::Debug;

use serde::{Deserialize, Serialize};

use crate::linalg::field::Field;



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
}

pub trait FieldAlgebra<F: Field>: CRing {
    fn scalar_mult(self, f: F) -> Self; 
}

pub trait ValuationRing: CRing {
    /// true iff self | other
    /// if false then other | self
    fn divides(self, other: &Self) -> bool;
    fn unsafe_divide(self, div: Self) -> Self;
}



#[derive(Debug, Clone, Copy, Deserialize, Serialize)]
pub struct UniPolRing<F: Field>(pub F, pub usize);


impl<F: Field> PartialEq for UniPolRing<F> {
    fn eq(&self, other: &Self) -> bool {
        (self.is_zero() && other.is_zero()) || self.0 == other.0 && self.1 == other.1
    }
}

impl<F: Field> CRing for UniPolRing<F> {
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    fn one() -> Self {
        Self(F::one(), 0)
    }

    fn zero() -> Self {
        Self(F::zero(), 0)
    }

    fn parse(input: &str) -> Result<Self, String> {
        if input.trim().is_empty() {
            return Ok(Self::zero());
        }
        
        // For now, parse as field element with degree 0
        let field_val = F::parse(input)?;
        Ok(Self(field_val, 0))
    }
    
    fn is_unit(&self) -> bool {
        self.0.is_unit() && (self.1 == 0)
    }
}

impl<F: Field> FieldAlgebra<F> for UniPolRing<F> {
    fn scalar_mult(self, f: F) -> Self {
        Self(self.0 * f, self.1)
    }
}

impl<F: Field> Default for UniPolRing<F> {
    fn default() -> Self {
        Self(F::zero(), 0)
    }
}

impl<F: Field> ValuationRing for UniPolRing<F> {
    /// For true, self | other
    /// For false, other | self
    fn divides(self, other: &Self) -> bool {
        if self.1 <= other.1 || other.0.is_zero() {
            true
        } else {
            false
        }
    }
    
    /// Self / div
    fn unsafe_divide(self, div: Self) -> Self {
        UniPolRing(self.0 * div.0.inv().unwrap(), self.1.checked_sub(div.1).unwrap())
    }
}  

impl<F: Field> Add for UniPolRing<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        // Add coefficients if same degree, otherwise need polynomial representation
        if self.1 == rhs.1 {
            Self(self.0 + rhs.0, self.1)
        } else {
            // For simplicity, this assumes we're only working with monomials
            // A full implementation would need polynomial representation
            if self.0.is_zero() {
                rhs
            } else if rhs.0.is_zero() {
                self
            } else {
                // Cannot add different degree terms without full polynomial
                panic!("Cannot add different degree monomials without polynomial representation")
            }
        }
    }
}

impl<F: Field> Sub for UniPolRing<F> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<F: Field> Neg for UniPolRing<F> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0, self.1)
    }
}

impl<F: Field> Mul for UniPolRing<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0, self.1 + rhs.1)
    }
}

impl<F: Field> AddAssign for UniPolRing<F> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<F: Field> SubAssign for UniPolRing<F> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<F: Field> MulAssign for UniPolRing<F> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<F: Field> std::iter::Sum for UniPolRing<F> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}