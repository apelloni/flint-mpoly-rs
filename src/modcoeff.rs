use crate::{QCoeff, ZCoeff};
use flint_sys::fmpz::*;
use flint_sys::fmpz_mod::*;
use num::BigInt;
use std::fmt;
use std::mem::MaybeUninit;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::result::Result;
use std::str::FromStr;
use std::convert::TryInto;

/// Container for arbitrary precision rational numbers
#[derive(Debug)]
pub struct ModCoeff {
    pub raw: fmpz,                // coefficient container
    pub ctx: fmpz_mod_ctx_struct, // context container
    pub modulo: u64,
    _res: fmpz, // coefficient container
}

impl ModCoeff {
    /// Default value for ModCoeff is zero
    pub fn new(modulo: &u64) -> Self {
        // Set modulo
        let n = ZCoeff::from_int((*modulo).try_into().expect("modulo not fitting in i64"));
        unsafe {
            // Set raw and ctx
            let mut raw = MaybeUninit::uninit();
            let mut ctx = MaybeUninit::uninit();
            let mut res = MaybeUninit::uninit();
            fmpz_init(raw.as_mut_ptr());
            fmpz_init(res.as_mut_ptr());
            fmpz_mod_ctx_init(ctx.as_mut_ptr(), &n as *const _ as *mut _);
            ModCoeff {
                raw: raw.assume_init(),
                ctx: ctx.assume_init(),
                modulo: *modulo,
                _res: res.assume_init(),
            }
        }
    }
    /// Initialize coefficient to to zero
    /// ```
    /// use flint_mpoly::ModCoeff;
    /// let coeff = ModCoeff::zero(&7);
    /// assert_eq!("0",coeff.to_str());
    /// ```
    pub fn zero(modulo: &u64) -> Self {
        // Set modulo
        let n = ZCoeff::from_int((*modulo).try_into().expect("modulo not fitting in i64"));
        unsafe {
            // Set raw and ctx
            let mut raw = MaybeUninit::uninit();
            let mut ctx = MaybeUninit::uninit();
            let mut res = MaybeUninit::uninit();
            fmpz_init(raw.as_mut_ptr());
            fmpz_init(res.as_mut_ptr());
            fmpz_mod_ctx_init(ctx.as_mut_ptr(), &n as *const _ as *mut _);
            fmpz_zero(raw.as_mut_ptr());
            ModCoeff {
                raw: raw.assume_init(),
                ctx: ctx.assume_init(),
                modulo: *modulo,
                _res: res.assume_init(),
            }
        }
    }

    /// Initialize coefficient to one
    /// ```
    /// use flint_mpoly::ModCoeff;
    /// let coeff = ModCoeff::one(&7);
    /// assert_eq!("1",coeff.to_str());
    /// ```
    pub fn one(modulo: &u64) -> Self {
        // Set modulo
        let n = ZCoeff::from_int((*modulo).try_into().expect("modulo not fitting in i64"));
        unsafe {
            // Set raw and ctx
            let mut raw = MaybeUninit::uninit();
            let mut ctx = MaybeUninit::uninit();
            let mut res = MaybeUninit::uninit();
            fmpz_init(raw.as_mut_ptr());
            fmpz_init(res.as_mut_ptr());
            fmpz_mod_ctx_init(ctx.as_mut_ptr(), &n as *const _ as *mut _);
            fmpz_one(raw.as_mut_ptr());
            ModCoeff {
                raw: raw.assume_init(),
                ctx: ctx.assume_init(),
                modulo: *modulo,
                _res: res.assume_init(),
            }
        }
    }

    /// Import the Q value from a str
    /// ```
    /// use flint_mpoly::ModCoeff;
    /// let modulo = 17;
    /// let mut coeff = ModCoeff::from_str("-1283719293715117894283698/28166512",&modulo).expect("bad string");
    /// assert_eq!("4",coeff.to_str());
    /// ```
    pub fn from_str(rat: &str, modulo: &u64) -> Result<Self, String> {
        // Set modulo
        let n = ZCoeff::from_int((*modulo).try_into().expect("modulo not fitting in i64"));
        // Convert to QCoeff fraction
        let qres = QCoeff::from_str(rat)?;
        unsafe {
            let mut raw = MaybeUninit::uninit();
            let mut ctx = MaybeUninit::uninit();
            let mut res = MaybeUninit::uninit();
            fmpz_init(raw.as_mut_ptr());
            fmpz_init(res.as_mut_ptr());
            fmpz_mod_ctx_init(ctx.as_mut_ptr(), &n as *const _ as *mut _);
            fmpz_mod_divides(
                res.as_mut_ptr(),
                &qres.raw.num as *const _ as *mut _,
                &qres.raw.den as *const _ as *mut _,
                ctx.as_mut_ptr(),
            );
            // Make canonical
            fmpz_mod_set_fmpz(raw.as_mut_ptr(), res.as_mut_ptr(), ctx.as_mut_ptr());
            Ok(ModCoeff {
                raw: raw.assume_init(),
                ctx: ctx.assume_init(),
                modulo: *modulo,
                _res: res.assume_init(),
            })
        }
    }

    /// Initialize coefficient to the canonical form of the fraction p / q
    /// ```
    /// use flint_mpoly::ModCoeff;
    /// let modulo = 7;
    /// let mut coeff = ModCoeff::from_int(12, &modulo);
    /// assert_eq!("5",coeff.to_str());
    /// ```
    pub fn from_int(c: i64, modulo: &u64) -> Self {
        // Set modulo
        let n = ZCoeff::from_int((*modulo).try_into().expect("modulo not fitting in i64"));
        unsafe {
            // Set raw and ctx
            let mut raw = MaybeUninit::uninit();
            let mut ctx = MaybeUninit::uninit();
            let mut res = MaybeUninit::uninit();
            fmpz_init(raw.as_mut_ptr());
            fmpz_init(res.as_mut_ptr());
            fmpz_mod_ctx_init(ctx.as_mut_ptr(), &n as *const _ as *mut _);
            fmpz_set_si(res.as_mut_ptr(), c);
            // Make canonical
            fmpz_mod_set_fmpz(raw.as_mut_ptr(), res.as_mut_ptr(), ctx.as_mut_ptr());
            ModCoeff {
                raw: raw.assume_init(),
                ctx: ctx.assume_init(),
                modulo: *modulo,
                _res: res.assume_init(),
            }
        }
    }

    /// Set a coefficient to the canonical form of the fraction p / q.
    /// ```
    /// use flint_mpoly::ModCoeff;
    /// let mut coeff = ModCoeff::new(&7);
    /// coeff.set_from_int(12);
    /// assert_eq!("5",coeff.to_str());
    /// ```
    pub fn set_from_int(&mut self, c: i64) {
        unsafe {
            fmpz_set_si(&self.raw as *const _ as *mut _, c);
        }
        self.canonical();
    }

    /// Set a coefficient from string for arbitrary fractions
    /// ```
    /// use flint_mpoly::ModCoeff;
    /// let mut coeff = ModCoeff::new(&17);
    /// coeff.set_from_str("-1283719293715117894283698/28166512").expect("bad string");
    /// assert_eq!("4",coeff.to_str());
    /// ```
    pub fn set_from_str(&mut self, rat: &str) -> Result<(), String> {
        // Convert to QCoeff fraction
        let qres = QCoeff::from_str(rat)?;
        unsafe {
            fmpz_mod_divides(
                &mut self.raw as *mut _,
                &qres.raw.num as *const _ as *mut _,
                &qres.raw.den as *const _ as *mut _,
                &mut self.ctx as *mut _,
            );
            Ok(())
        }
    }

    /// Check if the coefficient is zero
    /// ```
    /// use flint_mpoly::ModCoeff;
    /// let modulo = 19;
    /// let coeff = ModCoeff::from_int(0,&modulo);
    /// assert!(coeff.is_zero());
    /// ```
    pub fn is_zero(&self) -> bool {
        unsafe { fmpz_is_zero(&self.raw as *const _ as *mut _) == 1 }
    }

    /// Check if the coefficient is one
    /// ```
    /// use flint_mpoly::ModCoeff;
    /// let modulo = 19;
    /// let coeff = ModCoeff::from_int(20,&modulo);
    /// assert!(coeff.is_one());
    /// ```
    pub fn is_one(&self) -> bool {
        unsafe {
            fmpz_mod_is_one(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }

    /// Map the value to fit in [0,modulo)
    fn canonical(&mut self) {
        unsafe {
            fmpz_mod_set_fmpz(
                &mut self._res as *mut _,
                &mut self.raw as *mut _,
                &mut self.ctx as *mut _,
            );
            fmpz_swap(&mut self.raw as *mut _, &mut self._res as *mut _);
        }
    }

    /// Rise to integer power
    /// ```
    /// use flint_mpoly::ModCoeff;
    /// let modulo = 7;
    /// let mut coeff = ModCoeff::from_str("-1/2",&modulo).expect("bad string");
    /// coeff.pown(12);
    /// assert_eq!("1",coeff.to_str());
    /// ```
    pub fn pown(&mut self, n: u64) {
        unsafe {
            fmpz_mod_pow_ui(
                &mut self._res as *mut _,
                &mut self.raw as *mut _,
                n,
                &mut self.ctx as *mut _,
            );
            fmpz_swap(&mut self.raw as *mut _, &mut self._res as *mut _);
        }
    }

    /// Map to inverse a^-1 = b mod n
    pub fn inv(&mut self) {
        unsafe {
            fmpz_mod_inv(
                &mut self._res as *mut _,
                &mut self.raw as *mut _,
                &mut self.ctx as *mut _,
            );
            fmpz_swap(&mut self._res as *mut _, &mut self.raw as *mut _);
        }
    }

    /// Format to human readable string
    pub fn to_str(&self) -> String {
        format!("{}", self).to_string()
    }

    /// Swap content with another ModCoeff
    pub fn swap(&mut self, other: &mut ModCoeff) {
        assert!(
            self.modulo == other.modulo,
            "Unmatch modulo while swapping expressions"
        );
        unsafe {
            fmpz_swap(&mut self.raw as *mut _, &mut other.raw as *mut _);
        }
    }
}

/// Clone for ModCoeff
/// ```
/// use flint_mpoly::ModCoeff;
/// use std::str::FromStr;
/// let modulo = 7;
/// let c1 = ModCoeff::from_str( "123/321",&modulo).unwrap();
/// let c2 = c1.clone();
/// assert_eq!(c1.to_str(),c2.to_str());
/// ```
impl Clone for ModCoeff {
    fn clone(&self) -> Self {
        let mut other = ModCoeff::new(&self.modulo);
        unsafe {
            fmpz_set(&mut other.raw as *mut _, &self.raw as *const _ as *mut _);
        }
        other
    }
}

///Implement addition with operator `+`
impl<'a, 'b> Add<&'b ModCoeff> for &'a ModCoeff {
    type Output = ModCoeff;

    fn add(self, other: &'b ModCoeff) -> ModCoeff {
        assert!(
            self.modulo == other.modulo,
            "Unmatch modulo while adding expressions"
        );
        let mut this = ModCoeff::new(&self.modulo);
        unsafe {
            fmpz_mod_add(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

///Implement addition with operator `-`
impl<'a, 'b> Sub<&'b ModCoeff> for &'a ModCoeff {
    type Output = ModCoeff;

    fn sub(self, other: &'b ModCoeff) -> ModCoeff {
        assert!(
            self.modulo == other.modulo,
            "Unmatch modulo while subtracting expressions"
        );
        let mut this = ModCoeff::new(&self.modulo);
        unsafe {
            fmpz_mod_sub(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

///Implement addition with operator `*`
impl<'a, 'b> Mul<&'b ModCoeff> for &'a ModCoeff {
    type Output = ModCoeff;

    fn mul(self, other: &'b ModCoeff) -> ModCoeff {
        assert!(
            self.modulo == other.modulo,
            "Unmatch modulo while multiplying expressions"
        );
        let mut this = ModCoeff::new(&self.modulo);
        unsafe {
            fmpz_mod_mul(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

/// Implement addition with operator `/`
impl<'a, 'b> Div<&'b ModCoeff> for &'a ModCoeff {
    type Output = ModCoeff;

    fn div(self, other: &'b ModCoeff) -> ModCoeff {
        assert!(
            self.modulo == other.modulo,
            "Unmatch modulo while dividing expressions"
        );
        let mut this = ModCoeff::new(&self.modulo);
        unsafe {
            assert!(
                fmpz_mod_divides(
                    &mut this.raw as *mut _,
                    &self.raw as *const _ as *mut _,
                    &other.raw as *const _ as *mut _,
                    &self.ctx as *const _ as *mut _,
                ) == 1,
                "Failed to divied mod {}",
                self.modulo
            );
        }
        this
    }
}

///Implement negative sign
/// ```
/// use flint_mpoly::ModCoeff;
/// use std::str::FromStr;
/// let modulo = 17;
/// let c = ModCoeff::from_str( "13/7",&modulo).unwrap();
/// let c_neg = -&c;
/// assert_eq!("3", c_neg.to_str());
/// ```
impl<'a> Neg for &'a ModCoeff {
    type Output = ModCoeff;
    fn neg(self) -> ModCoeff {
        let mut this = ModCoeff::new(&self.modulo);
        unsafe {
            fmpz_mod_neg(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

/// Clear the content of all raw pointers before dropping ModCoeff
impl Drop for ModCoeff {
    fn drop(&mut self) {
        unsafe {
            fmpz_clear(&mut self.raw as *mut _);
            fmpz_clear(&mut self._res as *mut _);
            fmpz_mod_ctx_clear(&mut self.ctx as *mut _);
        }
    }
}

// To use the `{}` for the structure ModCoeff
impl fmt::Display for ModCoeff {
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
        let modulo = 11311;
        let coeff_str = "1-1";
        ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
    }
    #[test]
    #[should_panic]
    fn bad_string_2() {
        let modulo = 11311;
        let coeff_str = "1+1";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        println!("{} -> {}", coeff_str, coeff);
    }
    #[test]
    fn sign() {
        let modulo = 11311;
        let coeff_str = "+1";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        assert_eq!(coeff_str, format!("{coeff:+}").as_str());
    }
    #[test]
    fn unsign() {
        let modulo = 11311;
        let coeff_str = "1";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        assert_eq!(coeff_str, format!("{coeff}").as_str());
    }
    #[test]
    fn zero() {
        let modulo = 11311;
        let coeff_str = "0";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn positive_1() {
        let modulo = 11311;
        let coeff_str = "1";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        println!("{:?}", coeff.raw);
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn positive_2() {
        let modulo = 11311;
        let coeff_str = "2";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        assert_eq!(coeff_str, coeff.to_str());
    }
    #[test]
    fn positive_3() {
        let modulo = 11311;
        let coeff_str = "+1";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        println!("{:?}", coeff.raw);
        assert_eq!(coeff_str, format!("{coeff:+}").as_str());
    }
    #[test]
    fn positive_4() {
        let modulo = 11311;
        let coeff_str = "+2";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        assert_eq!(coeff_str, format!("{coeff:+}").as_str());
    }
    #[test]
    fn negative_1() {
        let modulo = 11311;
        let coeff_str = "-1";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        assert_eq!(11310.to_string(), coeff.to_str());
    }
    #[test]
    fn negative_2() {
        let modulo = 11311;
        let coeff_str = "-2";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        assert_eq!(11309.to_string(), coeff.to_str());
    }
    #[test]
    fn positive_long_1() {
        let modulo = 11311;
        let coeff_str = "12345678901234567890123456789";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        assert_eq!(2352.to_string(), coeff.to_str());
    }
    #[test]
    fn positive_long_2() {
        let modulo = 11311;
        let coeff_str = "123456789012345678901234567890";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        assert_eq!(898.to_string(), coeff.to_str());
    }
    #[test]
    fn negative_long_1() {
        let modulo = 11311;
        let coeff_str = "-12345678901234567890123456789";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        assert_eq!(8959.to_string(), coeff.to_str());
    }
    #[test]
    fn negative_long_2() {
        let modulo = 11311;
        let coeff_str = "-123456789012345678901234567890";
        let coeff = ModCoeff::from_str(coeff_str, &modulo).expect("bad string");
        assert_eq!(10413.to_string(), coeff.to_str());
    }
}
