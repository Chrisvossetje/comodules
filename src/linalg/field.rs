use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

pub trait CRing:
    Clone
    + Copy
    + Debug
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
    fn zero() -> Self;

    fn parse(input: &str) -> Result<Self, String>;

    fn dot_product(l: &Vec<Self>, r: &Vec<Self>) -> Self {
        l.iter().zip(r.iter()).map(|(x, y)| *x * *y).sum()
    }
}

pub trait Field: CRing {
    fn inv(self) -> Option<Self>;
    fn get_characteristic() -> usize;

    fn as_usize(self) -> usize;
}

impl CRing for f64 {
    fn is_zero(&self) -> bool {
        *self == 0.0
    }

    fn one() -> Self {
        1.0
    }

    fn zero() -> Self {
        0.0
    }

    fn parse(input: &str) -> Result<Self, String> {
        input
            .parse()
            .map_err(|_| format!("Field: {} could not be parsed", input))
    }
}

impl Field for f64 {
    fn inv(self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(1.0 / self)
        }
    }

    fn get_characteristic() -> usize {
        0
    }

    fn as_usize(self) -> usize {
        self as usize
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

impl<const P: u8> std::iter::Sum for Fp<P> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut sum: u64 = 0;
        for a in iter {
            sum += a.0 as u64;
        }
        Fp((sum % (P as u64)) as u8)
    }
}

impl<const P: u8> CRing for Fp<P> {
    fn is_zero(&self) -> bool {
        self.0 == 0
    }

    fn one() -> Self {
        Fp(1)
    }

    fn zero() -> Self {
        Fp(0)
    }

    fn parse(input: &str) -> Result<Self, String> {
        let chr: i32 = input
            .parse()
            .map_err(|_| format!("Field: {} could not be parsed", input))?;
        Ok(Self((chr % (P as i32)) as u8))
    }
}

impl<const P: u8> Field for Fp<P> {
    fn inv(self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            match P {
                2 | 3 => Some(self),
                _ => {
                    let mut result: u64 = 1;
                    let mut exp = P - 2;
                    let mut base = self.0 as u64;
                    while exp > 0 {
                        if exp % 2 == 1 {
                            result = (result * base) % P as u64;
                        }
                        base = (base * base) % P as u64;
                        exp /= 2;
                    }
                    Some(Fp(result as u8))
                }
            }
        }
    }

    fn get_characteristic() -> usize {
        P as usize
    }

    fn as_usize(self) -> usize {
        self.0 as usize
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

impl std::iter::Sum for F2 {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut sum: u32 = 0;
        for a in iter {
            sum += a.0 as u32;
        }
        F2((sum & 1) as u8)
    }
}

impl CRing for F2 {
    fn is_zero(&self) -> bool {
        self.0 == 0
    }

    fn one() -> F2 {
        F2(1)
    }

    fn zero() -> F2 {
        F2(0)
    }

    fn parse(input: &str) -> Result<Self, String> {
        let chr: i32 = input
            .parse()
            .map_err(|_| format!("Field: {} could not be parsed", input))?;
        Ok(Self((chr & 1) as u8))
    }
}

impl Field for F2 {
    fn inv(self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(self)
        }
    }

    fn get_characteristic() -> usize {
        2
    }

    fn as_usize(self) -> usize {
        self.0 as usize
    }
}
