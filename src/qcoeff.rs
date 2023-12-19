use flint_sys::fmpq::*;
use flint_sys::fmpz::*;
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
pub struct QCoeff {
    pub raw: fmpq, // coefficient container
}

impl QCoeff {
    /// Initialize coefficient to to zero
    /// ```
    /// use flint_mpoly::QCoeff;
    /// let coeff = QCoeff::zero();
    /// assert_eq!("0",coeff.to_str());
    /// ```
    pub fn zero() -> Self {
        unsafe {
            let mut c = MaybeUninit::uninit();
            fmpq_init(c.as_mut_ptr());
            fmpq_zero(c.as_mut_ptr());
            QCoeff {
                raw: c.assume_init(),
            }
        }
    }

    /// Initialize coefficient to one
    /// ```
    /// use flint_mpoly::QCoeff;
    /// let coeff = QCoeff::one();
    /// assert_eq!("1",coeff.to_str());
    /// ```
    pub fn one() -> Self {
        unsafe {
            let mut c = MaybeUninit::uninit();
            fmpq_init(c.as_mut_ptr());
            fmpq_one(c.as_mut_ptr());
            QCoeff {
                raw: c.assume_init(),
            }
        }
    }

    /// Initialize coefficient to the canonical form of the fraction p / q
    /// ```
    /// use flint_mpoly::QCoeff;
    /// let coeff = QCoeff::from_int(12,8);
    /// assert_eq!("3/2",coeff.to_str());
    /// ```
    pub fn from_int(p: i64, q: u64) -> Self {
        unsafe {
            let mut c = MaybeUninit::uninit();
            fmpq_init(c.as_mut_ptr());
            fmpq_set_si(c.as_mut_ptr(), p, q);
            QCoeff {
                raw: c.assume_init(),
            }
        }
    }

    /// Set a coefficient to the canonical form of the fraction p / q.
    /// ```
    /// use flint_mpoly::QCoeff;
    /// let mut coeff = QCoeff::default();
    /// coeff.set_from_int(12,8);
    /// assert_eq!("3/2",coeff.to_str());
    /// ```
    pub fn set_from_int(&mut self, p: i64, q: u64) {
        unsafe {
            fmpq_set_si(&self.raw as *const _ as *mut _, p, q);
        }
    }

    /// Set a coefficient from string for arbitrary fractions
    /// ```
    /// use flint_mpoly::QCoeff;
    /// let mut coeff = QCoeff::default();
    /// coeff.set_from_str("-1283719293715117894283698/28166512").expect("bad string");
    /// assert_eq!("-641859646857558947141849/14083256",coeff.to_str());
    /// ```
    pub fn set_from_str(&mut self, rat: &str) -> Result<(), String> {
        // Convert to C str
        let c_str = CString::new(rat).unwrap();
        unsafe {
            if fmpq_set_str(&mut self.raw as *mut _, c_str.as_ptr(), 10) == 0 {
                fmpq_canonicalise(&mut self.raw as *mut _);
                Ok(())
            } else {
                Err(format!("Failed to parse string {rat}"))
            }
        }
    }

    /// Check if the coefficient is zero
    /// ```
    /// use flint_mpoly::QCoeff;
    /// let coeff = QCoeff::from_int(0,3);
    /// assert!(coeff.is_zero());
    /// ```
    pub fn is_zero(&self) -> bool {
        unsafe { fmpq_is_zero(&self.raw as *const _ as *mut _) == 1 }
    }

    /// Check if the coefficient is one
    /// ```
    /// use flint_mpoly::QCoeff;
    /// let coeff = QCoeff::from_int(3,3);
    /// assert!(coeff.is_one());
    /// ```
    pub fn is_one(&self) -> bool {
        unsafe { fmpq_is_one(&self.raw as *const _ as *mut _) == 1 }
    }

    /// Check if the coefficient is a whole number
    /// ```
    /// use flint_mpoly::QCoeff;
    /// let coeff = QCoeff::from_int(9,3);
    /// assert!(coeff.is_int());
    /// assert_eq!("3",coeff.to_str());
    /// ```
    pub fn is_int(&self) -> bool {
        unsafe { fmpz_is_one(&self.raw.den as *const _ as *mut _) == 1 }
    }

    /// Rise to integer power
    /// ```
    /// use std::str::FromStr;
    /// use flint_mpoly::QCoeff;
    /// let mut coeff = QCoeff::from_str("-1/2").expect("bad string");
    /// coeff.pown(-12);
    /// assert_eq!("4096",coeff.to_str());
    /// ```

    pub fn pown(&mut self, n: i64) {
        let mut res = QCoeff::default();
        unsafe {
            fmpq_pow_si(&mut res.raw as *mut _, &mut self.raw as *mut _, n);
            fmpq_swap(&mut res.raw as *mut _, &mut self.raw as *mut _);
        }
    }

    /// Format to human readable string
    pub fn to_str(&self) -> String {
        format!("{}", self).to_string()
    }
}

impl Default for QCoeff {
    /// Default value for QCoeff is zero
    fn default() -> Self {
        unsafe {
            let mut c = MaybeUninit::uninit();
            fmpq_init(c.as_mut_ptr());
            fmpq_zero(c.as_mut_ptr());
            QCoeff {
                raw: c.assume_init(),
            }
        }
    }
}

impl FromStr for QCoeff {
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
            fmpq_init(c.as_mut_ptr());
            if fmpq_set_str(c.as_mut_ptr(), c_str.as_ptr(), 10) == 0 {
                fmpq_canonicalise(c.as_mut_ptr());
                Ok(QCoeff {
                    raw: c.assume_init(),
                })
            } else {
                Err(format!("Failed to parse string {rat}"))
            }
        }
    }
}

/// Clone for QCoeff
/// ```
/// use flint_mpoly::QCoeff;
/// use std::str::FromStr;
/// let c1 = QCoeff::from_str("123/321").unwrap();
/// let c2 = c1.clone();
/// assert_eq!(c1.to_str(),c2.to_str());
/// ```
impl Clone for QCoeff {
    fn clone(&self) -> Self {
        let mut other = QCoeff::default();
        unsafe {
            fmpq_set(&mut other.raw as *mut _, &self.raw as *const _ as *mut _);
        }
        other
    }
}

///Implement addition with operator `+`
impl<'a, 'b> Add<&'b QCoeff> for &'a QCoeff {
    type Output = QCoeff;

    fn add(self, other: &'b QCoeff) -> QCoeff {
        let mut this = QCoeff::default();
        unsafe {
            fmpq_add(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
            );
        }
        this
    }
}

///Implement addition with operator `-`
impl<'a, 'b> Sub<&'b QCoeff> for &'a QCoeff {
    type Output = QCoeff;

    fn sub(self, other: &'b QCoeff) -> QCoeff {
        let mut this = QCoeff::default();
        unsafe {
            fmpq_sub(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
            );
        }
        this
    }
}

///Implement addition with operator `*`
impl<'a, 'b> Mul<&'b QCoeff> for &'a QCoeff {
    type Output = QCoeff;

    fn mul(self, other: &'b QCoeff) -> QCoeff {
        let mut this = QCoeff::default();
        unsafe {
            fmpq_mul(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
            );
        }
        this
    }
}

/// Implement addition with operator `/`
impl<'a, 'b> Div<&'b QCoeff> for &'a QCoeff {
    type Output = QCoeff;

    fn div(self, other: &'b QCoeff) -> QCoeff {
        let mut this = QCoeff::default();
        unsafe {
            fmpq_div(
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
/// use flint_mpoly::QCoeff;
/// use std::str::FromStr;
/// // Define new polynomial
/// let c = QCoeff::from_str("13/7").unwrap();
/// let c_neg = -&c;
/// assert_eq!("-13/7", c_neg.to_str());
/// ```
impl<'a> Neg for &'a QCoeff {
    type Output = QCoeff;
    fn neg(self) -> QCoeff {
        let mut this = QCoeff::default();
        unsafe {
            fmpq_neg(&mut this.raw as *mut _, &self.raw as *const _ as *mut _);
        }
        this
    }
}

/// Clear the content of all raw pointers before dropping QCoeff
impl Drop for QCoeff {
    fn drop(&mut self) {
        unsafe {
            fmpq_clear(&mut self.raw as *mut _);
        }
    }
}

// To use the `{}` for the structure QCoeff
impl fmt::Display for QCoeff {
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
        let mut den = BigInt::default();
        // Get Sign
        let sgn = unsafe { (1 + fmpz_sgn(&self.raw.num as *const _ as *mut _)) / 2 };
        // Set numerator
        unsafe {
            for i in 0..fmpz_bits(&self.raw.num as *const _ as *mut _) {
                num.set_bit(
                    i,
                    fmpz_tstbit(&self.raw.num as *const _ as *mut _, i) == sgn,
                );
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
        if !self.is_int() {
            // Set denominator
            unsafe {
                for i in 0..fmpz_bits(&self.raw.den as *const _ as *mut _) {
                    den.set_bit(i, fmpz_tstbit(&self.raw.den as *const _ as *mut _, i) == 1);
                }
            }
            write!(f, "/{}", den)?;
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
        QCoeff::from_str(coeff_str).expect("bad string");
    }
    #[test]
    #[should_panic]
    fn bad_string_2() {
        let coeff_str = "1+1";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        println!("{} -> {}", coeff_str, coeff);
    }
    #[test]
    fn sign() {
        let coeff_str = "+1";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, format!("{coeff:+}").as_str());
    }
    #[test]
    fn unsign() {
        let coeff_str = "1";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, format!("{coeff}").as_str());
    }
    #[test]
    fn zero() {
        let coeff_str = "0";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn positive_1() {
        let coeff_str = "1";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        println!("{:?}", coeff.raw);
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn positive_2() {
        let coeff_str = "2";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn positive_3() {
        let coeff_str = "+1";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        println!("{:?}", coeff.raw);
        assert_eq!(coeff_str, format!("{coeff:+}").as_str());
    }
    #[test]
    fn positive_4() {
        let coeff_str = "+2";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, format!("{coeff:+}").as_str());
    }
    #[test]
    fn negative_1() {
        let coeff_str = "-1";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn negative_2() {
        let coeff_str = "-2";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn positive_long_1() {
        let coeff_str = "12345678901234567890123456789";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn positive_long_2() {
        let coeff_str = "123456789012345678901234567890";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn negative_long_1() {
        let coeff_str = "-12345678901234567890123456789";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn negative_long_2() {
        let coeff_str = "-123456789012345678901234567890";
        let coeff = QCoeff::from_str(coeff_str).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
}
