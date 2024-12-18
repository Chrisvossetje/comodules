use core::panic;
use std::{fmt::Debug, ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign}};


pub trait Field : Clone + Copy + Debug + Sized + PartialEq + Add<Output = Self> + Sub<Output = Self> + Neg<Output = Self> + Mul<Output = Self> + AddAssign + SubAssign + MulAssign {
    fn inv(self) -> Option<Self>;
    fn get_characteristic(&self) -> usize;
    fn is_zero(&self) -> bool;
    fn one() -> Self;
    fn zero() -> Self;
}


impl Field for f64 {
    fn inv(self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(1.0f64 / self)
        }
    }

    fn get_characteristic(&self) -> usize {
        0
    }
    
    fn is_zero(&self) -> bool {
        return self.is_normal();
    }
    
    fn one() -> Self {
        1.0f64
    }
    
    fn zero() -> Self {
        0.0f64
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Fp<const P: u8>(pub u8);

impl<const P: u8> Add for Fp<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Fp(((self.0 as u64 + rhs.0 as u64) % P as u64) as u8)
    }
} 
impl<const P: u8> Mul for Fp<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Fp(((self.0 as u64 * rhs.0 as u64) % P as u64) as u8)
    }
} 

impl<const P: u8> Sub for Fp<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Fp(((P as u64 + self.0 as u64 - rhs.0 as u64) % P as u64) as u8)
    }
}

impl<const P: u8> Neg for Fp<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Fp(((P as u64 - self.0 as u64) % P as u64) as u8)
    }
} 

impl<const P: u8> AddAssign for Fp<P> {
    fn add_assign(&mut self, rhs: Self) {
        self.0 = ((self.0 as u64 + rhs.0 as u64) % P as u64) as u8;
    }
} 
impl<const P: u8> MulAssign for Fp<P> {
    fn mul_assign(&mut self, rhs: Self) {
        self.0 = ((self.0 as u64 * rhs.0 as u64) % P as u64) as u8;
    }
} 

impl<const P: u8> SubAssign for Fp<P> {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 = ((P as u64 + self.0 as u64 - rhs.0 as u64) % P as u64) as u8;
    }
} 

impl<const P: u8> std::fmt::Debug for Fp<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<const P: u8> Field for Fp<P> {
    fn inv(self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(Fp((self.0 as u64).pow(P as u32 - 2) as u8)) // a^(n-2) = a^(-1) mod P, Fermat's little theorem
        }
    }

    fn get_characteristic(&self) -> usize {
        P as usize
    }

    fn is_zero(&self) -> bool {
        self.0 == 0
    }

    fn one() -> Self {
        Fp(1)
    }
    
    fn zero() -> Self {
        Fp(0)
    }
}


#[derive(Clone, Copy, PartialEq, Eq)]
pub struct F2(pub u8);

impl Add for F2 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        F2(self.0 ^ rhs.0)
    }
}

impl Sub for F2 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self.add(rhs)
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
        F2(self.0 & rhs.0)
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
        F2((value & 1) as u8)
    }
}

impl std::fmt::Debug for F2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Field for F2 {
    fn get_characteristic(&self) -> usize {
        2
    }

    
    fn inv(self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(self)
        }
    }
    
    fn is_zero(&self) -> bool {
        self.0 == 0
    }

    fn one() -> F2 {
        F2(1)
    }

    fn zero() -> F2 {
        F2(0)
    }
}