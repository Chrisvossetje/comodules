use std::fmt::{Debug, Formatter};
use std::ops::{Add, Sub, Neg, Mul, AddAssign, SubAssign, MulAssign};

use deepsize::DeepSizeOf;
use serde::{Deserialize, Serialize};

use crate::{field::Field, ring::{CRing, ValuationRing}};



#[derive(Clone, Copy, Deserialize, Serialize, DeepSizeOf)]
pub struct UniPolRing<F: Field>(pub F, pub u16);


impl<F: Field> PartialEq for UniPolRing<F> {
    fn eq(&self, other: &Self) -> bool {
        (self.is_zero() && other.is_zero()) || self.0 == other.0 && self.1 == other.1
    }
}

impl<F: Field> CRing for UniPolRing<F> {
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    fn try_inverse(&self) -> Option<Self> {
        if self.1 > 0 {
            None
        } else {
            match self.0.inv() {
                Some(inv) => Some(Self(inv, 0)),
                None => None,
            }
        }
    }

    fn one() -> Self {
        Self(F::one(), 0)
    }

    fn zero() -> Self {
        Self(F::zero(), 0)
    }

    fn parse(input: &str) -> Result<Self, String> {
        let input = input.trim();
        if input.is_empty() {
            return Ok(Self::zero());
        }
        let (field_val, power) =  match input.split_once('t') {
            Some((lhs,p)) => {

                let field_el = if lhs.is_empty() {
                    F::one()
                } else {
                    F::parse(lhs)?
                };
                
                match p.split_once('^') {
                    Some((_,a)) => (field_el, str::parse(a).map_err(|_| {
                        format!("{a}, was not able to parse as i32")
                    })?),
                    None => {
                        (field_el,1)
                    },
                }
            }
            ,
            None => {
                (F::parse(input)?, 0)
            },
        };
        Ok(Self(field_val, power))
    }
    
    fn is_unit(&self) -> bool {
        self.0.is_unit() && (self.1 == 0)
    }
    
    fn unsafe_inverse(&self) -> Self {
        Self(self.0.unsafe_inverse(), 0)
    }
}

impl<F: Field> Default for UniPolRing<F> {
    fn default() -> Self {
        Self(F::zero(), 0)
    }
}

impl<F: Field> Debug for UniPolRing<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        if self.is_zero() {
            write!(f, "0")
        } else {
            if self.1 == 0 {
                write!(f, "1")
            } else {
                write!(f, "t^{:?}", self.1)
            }
        }
    }
}


impl<F: Field> ValuationRing for UniPolRing<F> {
    /// For true, self | other
    /// For false, other | self
    fn divides(self, other: &Self) -> bool {
        if !self.is_zero() && (self.1 <= other.1 || other.is_zero()) {
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