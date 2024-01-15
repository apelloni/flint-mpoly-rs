use crate::flint_sys::flint::*;
use crate::flint_sys::fmpz::*;
use num::BigInt;
use regex::Regex;
use std::ffi::CString;
use std::fmt;
use std::mem::MaybeUninit;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::result::Result;
use std::str::FromStr;

/// Container for arbitrary precision rational numbers
#[derive(Debug)]
pub struct ZCoeff {
    pub raw: fmpz, // coefficient container
}

impl ZCoeff {
    /// Initialize coefficient to to zero
    /// ```
    /// use flint_mpoly::ZCoeff;
    /// let coeff = ZCoeff::zero();
    /// assert_eq!("0",coeff.to_str());
    /// ```
    pub fn zero() -> Self {
        unsafe {
            let mut c = MaybeUninit::uninit();
            fmpz_init(c.as_mut_ptr());
            fmpz_zero(c.as_mut_ptr());
            ZCoeff {
                raw: c.assume_init(),
            }
        }
    }

    /// Initialize coefficient to one
    /// ```
    /// use flint_mpoly::ZCoeff;
    /// let coeff = ZCoeff::one();
    /// assert_eq!("1",coeff.to_str());
    /// ```
    pub fn one() -> Self {
        unsafe {
            let mut c = MaybeUninit::uninit();
            fmpz_init(c.as_mut_ptr());
            fmpz_one(c.as_mut_ptr());
            ZCoeff {
                raw: c.assume_init(),
            }
        }
    }

    /// Initialize coefficient to the canonical form of the fraction p / q
    /// ```
    /// use flint_mpoly::ZCoeff;
    /// let coeff = ZCoeff::from_int(12);
    /// assert_eq!("12",coeff.to_str());
    /// ```
    pub fn from_int(n: i64) -> Self {
        unsafe {
            let mut c = MaybeUninit::uninit();
            fmpz_init(c.as_mut_ptr());
            fmpz_set_si(c.as_mut_ptr(), n);
            ZCoeff {
                raw: c.assume_init(),
            }
        }
    }

    /// Set a coefficient to the canonical form of the fraction p / q.
    /// ```
    /// use flint_mpoly::ZCoeff;
    /// let mut coeff = ZCoeff::default();
    /// coeff.set_from_int(-12);
    /// assert_eq!("-12",coeff.to_str());
    /// ```
    pub fn set_from_int(&mut self, n: i64) {
        unsafe {
            fmpz_set_si(&self.raw as *const _ as *mut _, n);
        }
    }

    /// Set a coefficient from string for arbitrary fractions
    /// ```
    /// use flint_mpoly::ZCoeff;
    /// let mut coeff = ZCoeff::default();
    /// coeff.set_from_str("-1283719293715117894283698").expect("bad string");
    /// assert_eq!("-1283719293715117894283698",coeff.to_str());
    /// ```
    pub fn set_from_str(&mut self, rat: &str) -> Result<(), String> {
        // Convert to C str
        let c_str = CString::new(rat).unwrap();
        unsafe {
            if fmpz_set_str(&mut self.raw as *mut _, c_str.as_ptr(), 10) == 0 {
                Ok(())
            } else {
                Err(format!("Failed to parse string {rat}"))
            }
        }
    }

    /// Check if the coefficient is zero
    /// ```
    /// use flint_mpoly::ZCoeff;
    /// let coeff = ZCoeff::from_int(0);
    /// assert!(coeff.is_zero());
    /// ```
    pub fn is_zero(&self) -> bool {
        unsafe { fmpz_is_zero(&self.raw as *const _ as *mut _) == 1 }
    }

    /// Check if the coefficient is one
    /// ```
    /// use flint_mpoly::ZCoeff;
    /// let coeff = ZCoeff::from_int(1);
    /// assert!(coeff.is_one());
    /// ```
    pub fn is_one(&self) -> bool {
        unsafe { fmpz_is_one(&self.raw as *const _ as *mut _) == 1 }
    }

    /// Rise to integer power
    /// ```
    /// use std::str::FromStr;
    /// use flint_mpoly::ZCoeff;
    /// let mut coeff = ZCoeff::from_str("-2").expect("bad string");
    /// coeff.pown(12);
    /// assert_eq!("4096",coeff.to_str());
    /// ```
    pub fn pown(&mut self, n: u64) {
        let mut res = ZCoeff::default();
        unsafe {
            fmpz_pow_ui(&mut res.raw as *mut _, &mut self.raw as *mut _, n);
            fmpz_swap(&mut res.raw as *mut _, &mut self.raw as *mut _);
        }
    }

    /// Format to human readable string
    pub fn to_str(&self) -> String {
        format!("{}", self).to_string()
    }

    /// Swap content with another ZCoeff
    pub fn swap(&mut self, other: &mut ZCoeff) {
        unsafe {
            fmpz_swap(&mut self.raw as *mut _, &mut other.raw as *mut _);
        }
    }
}

impl Default for ZCoeff {
    /// Default value for ZCoeff is zero
    fn default() -> Self {
        unsafe {
            let mut c = MaybeUninit::uninit();
            fmpz_init(c.as_mut_ptr());
            fmpz_zero(c.as_mut_ptr());
            ZCoeff {
                raw: c.assume_init(),
            }
        }
    }
}

impl FromStr for ZCoeff {
    type Err = String;
    /// Import the Q value from a str
    ///Set a coefficient to the canonical form a big fraction p / q.
    fn from_str(rat: &str) -> Result<Self, Self::Err> {
        // TODO better check for the string sanity
        let rg_clean = Regex::new(r"^[\ ]*\+").unwrap();
        // Convert to C str
        let c_str = CString::new(rg_clean.replace(rat, "").to_string()).unwrap();
        unsafe {
            let mut c = MaybeUninit::uninit();
            fmpz_init(c.as_mut_ptr());
            if fmpz_set_str(c.as_mut_ptr(), c_str.as_ptr(), 10) == 0 {
                Ok(ZCoeff {
                    raw: c.assume_init(),
                })
            } else {
                Err(format!("Failed to parse string {rat}"))
            }
        }
    }
}

/// Clone for ZCoeff
/// ```
/// use flint_mpoly::ZCoeff;
/// use std::str::FromStr;
/// let c1 = ZCoeff::from_str("123").unwrap();
/// let c2 = c1.clone();
/// assert_eq!(c1.to_str(),c2.to_str());
/// ```
impl Clone for ZCoeff {
    fn clone(&self) -> Self {
        let mut other = ZCoeff::default();
        unsafe {
            fmpz_set(&mut other.raw as *mut _, &self.raw as *const _ as *mut _);
        }
        other
    }
}

///Implement addition with operator `+`
impl<'a, 'b> Add<&'b ZCoeff> for &'a ZCoeff {
    type Output = ZCoeff;

    fn add(self, other: &'b ZCoeff) -> ZCoeff {
        let mut this = ZCoeff::default();
        unsafe {
            fmpz_add(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
            );
        }
        this
    }
}

///Implement addition with operator `-`
impl<'a, 'b> Sub<&'b ZCoeff> for &'a ZCoeff {
    type Output = ZCoeff;

    fn sub(self, other: &'b ZCoeff) -> ZCoeff {
        let mut this = ZCoeff::default();
        unsafe {
            fmpz_sub(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
            );
        }
        this
    }
}

///Implement addition with operator `*`
impl<'a, 'b> Mul<&'b ZCoeff> for &'a ZCoeff {
    type Output = ZCoeff;

    fn mul(self, other: &'b ZCoeff) -> ZCoeff {
        let mut this = ZCoeff::default();
        unsafe {
            fmpz_mul(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
            );
        }
        this
    }
}

///Implement negative sign
/// ```
/// use flint_mpoly::ZCoeff;
/// use std::str::FromStr;
/// // Define new polynomial
/// let c = ZCoeff::from_str("13").unwrap();
/// let c_neg = -&c;
/// assert_eq!("-13", c_neg.to_str());
/// ```
impl<'a> Neg for &'a ZCoeff {
    type Output = ZCoeff;
    fn neg(self) -> ZCoeff {
        let mut this = ZCoeff::default();
        unsafe {
            fmpz_neg(&mut this.raw as *mut _, &self.raw as *const _ as *mut _);
        }
        this
    }
}

///Implement comparison ==
/// ```
/// use flint_mpoly::ZCoeff;
/// use std::str::FromStr;
/// // Define new polynomial
/// let c1 = ZCoeff::from_str("7").unwrap();
/// let c2 = ZCoeff::from_str("7").unwrap();
/// assert_eq!(c1,c2);
/// assert!(c1==c2);
/// ```
impl PartialEq for ZCoeff {
    fn eq(&self, other: &Self) -> bool {
        unsafe {
            fmpz_equal(
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
            ) == 1
        }
    }
}

/// Clear the content of all raw pointers before dropping ZCoeff
impl Drop for ZCoeff {
    fn drop(&mut self) {
        unsafe {
            fmpz_clear(&mut self.raw as *mut _);
            //flint_cleanup();
        }
    }
}

// To use the `{}` for the structure ZCoeff
impl fmt::Display for ZCoeff {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            if f.sign_plus() {
                write!(f, "{:+}", 0)?;
            } else {
                write!(f, "{}", 0)?;
            }
            return Ok(());
        }
        let mut num = BigInt::default();
        // Get Sign
        let sgn = unsafe { (1 + fmpz_sgn(&self.raw as *const _ as *mut _)) / 2 };
        // Set numerator
        unsafe {
            for i in 0..fmpz_bits(&self.raw as *const _ as *mut _) {
                num.set_bit(i, fmpz_tstbit(&self.raw as *const _ as *mut _, i) == sgn);
            }
            // If negative add correction
            if sgn == 0 {
                num = -(num + 1u64);
            }
        }
        if f.sign_plus() {
            write!(f, "{:+}", num)?;
        } else {
            write!(f, "{}", num)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    #[should_panic]
    fn bad_string_1() {
        let coeff_str = "1-1";
        ZCoeff::from_str(coeff_str).expect("bad string");
    }
    #[test]
    #[should_panic]
    fn bad_string_2() {
        let coeff_str = "1+1";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        println!("{} -> {}", coeff_str, coeff);
    }
    #[test]
    fn sign() {
        let coeff_str = "+1";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, format!("{coeff:+}").as_str());
    }
    #[test]
    fn unsign() {
        let coeff_str = "1";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, format!("{coeff}").as_str());
    }
    #[test]
    fn zero() {
        let coeff_str = "0";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn positive_1() {
        let coeff_str = "1";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        println!("{:?}", coeff.raw);
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn positive_2() {
        let coeff_str = "2";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn positive_3() {
        let coeff_str = "+1";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        println!("{:?}", coeff.raw);
        assert_eq!(coeff_str, format!("{coeff:+}").as_str());
    }
    #[test]
    fn positive_4() {
        let coeff_str = "+2";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, format!("{coeff:+}").as_str());
    }
    #[test]
    fn negative_1() {
        let coeff_str = "-1";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn negative_2() {
        let coeff_str = "-2";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn positive_long_1() {
        let coeff_str = "12345678901234567890123456789";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn positive_long_2() {
        let coeff_str = "123456789012345678901234567890";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn negative_long_1() {
        let coeff_str = "-12345678901234567890123456789";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn negative_long_2() {
        let coeff_str = "-123456789012345678901234567890";
        let coeff = ZCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
}
