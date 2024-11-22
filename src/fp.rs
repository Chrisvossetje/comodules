use std::{fmt::Debug, ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign}};


pub trait Field : Debug + Sized + PartialEq + Copy + Add<Output = Self> + Sub<Output = Self> + Neg<Output = Self> + Mul<Output = Self> + AddAssign + SubAssign + MulAssign {
    fn inv(self) -> Self;
    fn get_characteristic(&self) -> usize;
    fn is_zero(&self) -> bool;
    fn one() -> Self;
    fn zero() -> Self;
}


#[derive(Clone, Copy, PartialEq, Eq)]
pub struct F2 {
    value: u8
}

impl Add for F2 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        F2 { value: self.value ^ rhs.value }
    }
}

impl Sub for F2 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        F2 { value: self.value ^ rhs.value }
    }
    
}

impl Neg for F2 {
    type Output = Self;

    fn neg(self) -> Self {
        self
    }
    
}

impl Mul for F2 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        F2 { value: self.value & rhs.value }
    }
}

impl AddAssign for F2 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for F2 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl MulAssign for F2 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl From<u64> for F2 {
    fn from(value: u64) -> Self {
        F2 { value: (value % 2) as u8 }
    }
}

impl std::fmt::Debug for F2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value)
    }
}

impl Field for F2 {
    fn get_characteristic(&self) -> usize {
        2
    }

    
    fn inv(self) -> Self {
        self
    }
    
    fn is_zero(&self) -> bool {
        self.value == 0
    }

    fn one() -> F2 {
        F2 {
            value: 1
        }
    }

    fn zero() -> F2 {
        F2 {
            value: 0
        }
    }
}