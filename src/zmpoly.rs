use crate::parser::parse_mpoly;
use crate::{QCoeff, QMPoly, ZCoeff};
use flint_sys::flint::*;
use flint_sys::fmpq_mpoly::*;
use flint_sys::fmpz::*;
use flint_sys::fmpz_mpoly::*;
use flint_sys::mpoly::ordering_t_ORD_DEGLEX;
use regex::Regex;
use std::fmt;
use std::mem::MaybeUninit;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Container for polynomials over rational numbers
#[derive(Debug)]
pub struct ZMPoly {
    pub raw: fmpz_mpoly_struct,     // Polynomial
    pub ctx: fmpz_mpoly_ctx_struct, // Polynomial context
    pub vars: Vec<String>,          // Name of variables
    _coeff1: ZCoeff,                // coefficient container
    _coeff2: ZCoeff,                // coefficient container
    _res: ZCoeff,                   // coefficient container
    _exp: [u64; 10],                // exponent container
}

/// Container for polynomial that uses Flint as backend
impl ZMPoly {
    /// New 0 polynomial with selected variables
    /// ```
    /// use flint_mpoly::ZMPoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly = ZMPoly::new(&vars);
    /// assert_eq!("0",mpoly.to_str());
    /// ```
    pub fn new(vars: &[String]) -> Self {
        unsafe {
            let mut raw = MaybeUninit::uninit();
            let mut ctx = MaybeUninit::uninit();
            let nvars = vars.len() as i64;
            fmpz_mpoly_ctx_init(ctx.as_mut_ptr(), nvars, ordering_t_ORD_DEGLEX);
            fmpz_mpoly_init(raw.as_mut_ptr(), ctx.as_mut_ptr());
            ZMPoly {
                //raw: raw,
                raw: raw.assume_init(),
                ctx: ctx.assume_init(),
                vars: vars.to_vec(),
                _coeff1: ZCoeff::default(),
                _coeff2: ZCoeff::default(),
                _res: ZCoeff::default(),
                _exp: [0; 10],
            }
        }
    }

    /// Initialize polynomial from string assuming the given variables
    /// ```
    /// use flint_mpoly::ZMPoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly= ZMPoly::from_str("(10*x1+10*x2)", &vars).unwrap();
    /// assert_eq!("+10*x1+10*x2",mpoly.to_str());
    /// ```
    pub fn from_str(expr_str: &str, vars: &[String]) -> Result<Self, &'static str> {
        parse_mpoly(expr_str, vars)
            .expect("parse polynomial")
            .to_zmpoly()
    }

    /// Clone will copy the pointers and we will have two objects pointing at the same memory
    /// address. With this function we can create a clone in a safe way
    /// ```
    /// use flint_mpoly::ZMPoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly_1= ZMPoly::from_str("(x1+x2)", &vars).unwrap();
    /// let mpoly_2= ZMPoly::clone_from(&mpoly_1);
    /// assert_eq!(mpoly_1.to_str(), mpoly_2.to_str());
    /// ```
    pub fn clone_from(other: &Self) -> Self {
        let mut this = ZMPoly::new(&other.vars);
        unsafe {
            fmpz_mpoly_set(
                &mut this.raw as *mut _,
                &other.raw as *const _ as *mut _,
                &other.ctx as *const _ as *mut _,
            );
        }
        this
    }

    /// Set polynomial to zero
    pub fn set_to_zero(&mut self) {
        unsafe {
            fmpz_mpoly_zero(&mut self.raw as *mut _, &mut self.ctx as *mut _);
        }
    }

    /// Set polynomial to one
    pub fn set_to_one(&mut self) {
        unsafe {
            fmpz_mpoly_one(&mut self.raw as *mut _, &mut self.ctx as *mut _);
        }
    }

    /// Set polynomial from string assuming the variables already stored in ZMPoly
    /// ```
    /// use flint_mpoly::ZMPoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mpoly= ZMPoly::new(&vars);
    /// mpoly.set_from_str("x1+x2").unwrap();
    /// assert_eq!("+x1+x2",mpoly.to_str());
    /// ```
    pub fn set_from_str(&mut self, expr_str: &str) -> Result<(), &'static str> {
        *self = ZMPoly::from_str(expr_str, &self.vars)?;
        Ok(())
    }

    /// Check if the polynomial is zero
    pub fn is_zero(&self) -> bool {
        unsafe {
            fmpz_mpoly_is_zero(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }

    /// Check if the polynomial is one
    pub fn is_one(&self) -> bool {
        unsafe {
            fmpz_mpoly_is_one(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }
    /// Check if the polynomial is const
    pub fn is_const(&self) -> bool {
        unsafe {
            fmpz_mpoly_is_fmpz(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }

    /// Add to the polynomial the coefficient with the specified powers
    /// ```
    /// use flint_mpoly::ZMPoly;
    /// use flint_mpoly::ZCoeff;
    /// // Define new polynomial
    /// let vars = [String::from("x1"), String::from("x2")];
    /// let mut mpoly= ZMPoly::new(&vars);
    /// // Set coefficient of monomial with power x1^2*x2^3 and coefficient 3
    /// let pows = [2,3];
    /// let coeff = ZCoeff::from_int(3);
    /// mpoly.add_coeff(&pows,&coeff);
    /// assert_eq!("+3*x1^2*x2^3", mpoly.to_str());
    /// ```
    pub fn add_coeff(&mut self, monimial_powers: &[u64], coeff: &ZCoeff) {
        unsafe {
            fmpz_mpoly_get_coeff_fmpz_ui(
                &mut self._coeff1.raw as *mut _,
                &mut self.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
            if self._coeff1.is_zero() {
                fmpz_set(
                    &mut self._res.raw as *mut _,
                    &coeff.raw as *const _ as *mut _,
                )
            } else {
                fmpz_add(
                    &mut self._res.raw as *mut _,
                    &mut self._coeff1.raw as *mut _,
                    &coeff.raw as *const _ as *mut _,
                );
            }
            fmpz_mpoly_set_coeff_fmpz_ui(
                &mut self.raw as *mut _,
                &mut self._res.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
        }
    }

    /// Add to the polynomial the coefficient p/q with the specified powers
    /// Add to the polynomial the coefficient with the specified powers
    /// ```
    /// use flint_mpoly::ZMPoly;
    /// // Define new polynomial
    /// let vars = [String::from("x1"), String::from("x2")];
    /// let mut mpoly= ZMPoly::new(&vars);
    /// // Set coefficient of monomial with power x1^2*x2^3 and coefficient 3
    /// let pows = [2,3];
    /// mpoly.add_coeff_int(&pows,3);
    /// assert_eq!("+3*x1^2*x2^3", mpoly.to_str());
    /// ```
    pub fn add_coeff_int(&mut self, monimial_powers: &[u64], n: i64) {
        unsafe {
            // Store current coefficient to _coeff1
            fmpz_mpoly_get_coeff_fmpz_ui(
                &mut self._coeff1.raw as *mut _,
                &mut self.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
            // Set coefficient to add to _coeff2
            self._coeff2.set_from_int(n);
            if self._coeff1.is_zero() {
                fmpz_swap(
                    &mut self._res.raw as *mut _,
                    &mut self._coeff2.raw as *mut _,
                )
            } else {
                fmpz_add(
                    &mut self._res.raw as *mut _,
                    &mut self._coeff1.raw as *mut _,
                    &mut self._coeff2.raw as *mut _,
                );
            }
            fmpz_mpoly_set_coeff_fmpz_ui(
                &mut self.raw as *mut _,
                &mut self._res.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
        }
    }

    /// Add to the polynomial the coefficient p/q with the specified powers
    /// ```
    /// use flint_mpoly::ZMPoly;
    /// // Define new polynomial
    /// let vars = [String::from("x1"), String::from("x2")];
    /// let mut mpoly= ZMPoly::new(&vars);
    /// // Set coefficient of monomial with power x1^2*x2^3 and coefficient 3
    /// let pows = [2,3];
    /// mpoly.add_coeff_str(&pows,"3");
    /// assert_eq!("+3*x1^2*x2^3", mpoly.to_str());
    /// ```
    pub fn add_coeff_str(&mut self, monimial_powers: &[u64], coeff: &str) {
        unsafe {
            // Store current coefficient to _coeff1
            fmpz_mpoly_get_coeff_fmpz_ui(
                &mut self._coeff1.raw as *mut _,
                &mut self.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
            // Set coefficient to add to _coeff2
            self._coeff2.set_from_str(coeff).unwrap();
            if self._coeff1.is_zero() {
                fmpz_swap(
                    &mut self._res.raw as *mut _,
                    &mut self._coeff2.raw as *mut _,
                );
            } else {
                fmpz_add(
                    &mut self._res.raw as *mut _,
                    &mut self._coeff1.raw as *mut _,
                    &mut self._coeff2.raw as *mut _,
                );
            }
            fmpz_mpoly_set_coeff_fmpz_ui(
                &mut self.raw as *mut _,
                &mut self._res.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
        }
    }
    /// Rise polynomial to integer power
    pub fn pown(&mut self, n: u64) {
        let mut pow_res = ZMPoly::new(&self.vars);
        unsafe {
            fmpz_mpoly_pow_ui(
                &mut pow_res.raw as *mut _,
                &mut self.raw as *mut _,
                n,
                &mut self.ctx as *mut _,
            );
            fmpz_mpoly_swap(
                &mut pow_res.raw as *mut _,
                &mut self.raw as *mut _,
                &mut self.ctx as *mut _,
            );
        }
    }
    /// Swap the structure content with another ZMPoly provided they have the same context
    pub fn swap(&mut self, other: &mut Self) {
        assert!(
            self.vars == other.vars,
            "Unmatch variables while swapping expressions"
        );
        unsafe {
            fmpz_mpoly_swap(
                &mut self.raw as *mut _,
                &mut other.raw as *mut _,
                &mut self.ctx as *mut _,
            );
        }
    }
    /// Set function to zero and return the previous value
    /// ```
    /// use flint_mpoly::ZMPoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mpoly_a = ZMPoly::from_str("+x1+x2",&vars).unwrap();
    /// assert_eq!("+x1+x2",mpoly_a.to_str());
    /// let mut mpoly_b = mpoly_a.zero_move();
    /// assert_eq!("0",mpoly_a.to_str());
    /// assert_eq!("+x1+x2",mpoly_b.to_str());
    /// ```
    pub fn zero_move(&mut self) -> Self {
        let mut out = Self::new(&self.vars);
        out.swap(self);
        out
    }
    /// Clear the function and set it to zero
    pub fn clear(&mut self) {
        self.set_to_zero();
    }
    /// Clear all memory allocated to ZMPoly
    pub fn free(&mut self) {
        self.clear();
        unsafe {
            fmpz_mpoly_realloc(&mut self.raw as *mut _, 0, &mut self.ctx as *mut _);
        }
    }

    /// Convert the polynomial to a human readable string
    pub fn to_str(&self) -> String {
        let rg = Regex::new(r"([+-])1\*").unwrap();
        rg.replace_all(format!("{}", self).as_str(), "$1")
            .to_string()
    }

    /// Convert the QMPoly
    /// ```
    /// use flint_mpoly::{ZMPoly,QMPoly};
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let zpoly = ZMPoly::from_str("+2*x1+4*x2",&vars).unwrap();
    /// let qpoly:QMPoly = zpoly.to_qmpoly();
    /// assert_eq!("+2*x1+4*x2",qpoly.to_str());
    /// ```
    pub fn to_qmpoly(&self) -> QMPoly {
        let mut res = QMPoly::new(&self.vars);
        let mut c = QCoeff::default();
        unsafe {
            let length = fmpz_mpoly_length(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );

            for ci in 0..length {
                // get coeff
                fmpz_mpoly_get_term_coeff_fmpz(
                    &self._coeff1.raw as *const _ as *mut _,
                    &self.raw as *const _ as *mut _,
                    ci,
                    &self.ctx as *const _ as *mut _,
                );
                // get exponent
                fmpz_mpoly_get_term_exp_ui(
                    &self._exp as *const _ as *mut _,
                    &self.raw as *const _ as *mut _,
                    ci,
                    &self.ctx as *const _ as *mut _,
                );

                if fmpz_is_zero(&self._coeff1.raw as *const _ as *mut _) == 0 {
                    // Cast to fmpq
                    fmpz_swap(
                        &mut c.raw.num as *mut _,
                        &self._coeff1.raw as *const _ as *mut _,
                    );
                    // add to polynomial
                    fmpq_mpoly_set_coeff_fmpq_ui(
                        &mut res.raw as *mut _,
                        &mut c.raw as *mut _,
                        self._exp.as_ptr(),
                        &mut res.ctx as *mut _,
                    );
                }
            }
        }
        res
    }
}

///Implement addition with operator `+`
/// ```
/// use flint_mpoly::ZMPoly;
/// // Define new polynomial
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly_a = ZMPoly::from_str("x1",&vars).unwrap();
/// let mpoly_b = ZMPoly::from_str("x2",&vars).unwrap();
/// let mpoly_ab = &mpoly_a + &mpoly_b;
/// assert_eq!("+x1+x2", mpoly_ab.to_str());
/// ```
impl<'a, 'b> Add<&'b ZMPoly> for &'a ZMPoly {
    type Output = ZMPoly;

    fn add(self, other: &'b ZMPoly) -> ZMPoly {
        assert!(
            self.vars == other.vars,
            "Unmatch variables in binary operation"
        );
        let mut this = ZMPoly::new(&self.vars);
        unsafe {
            fmpz_mpoly_add(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
                &other.ctx as *const _ as *mut _,
            );
        }
        this
    }
}
///Implement addition with operator `-`
/// ```
/// use flint_mpoly::ZMPoly;
/// use flint_mpoly::parse_mpoly;
/// // Define new polynomial
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly_a = ZMPoly::from_str("x1",&vars).unwrap();
/// let mpoly_b = ZMPoly::from_str("x2",&vars).unwrap();
/// let mpoly_ab = &mpoly_a - &mpoly_b;
/// assert_eq!("+x1-x2", mpoly_ab.to_str());
/// ```
impl<'a, 'b> Sub<&'b ZMPoly> for &'a ZMPoly {
    type Output = ZMPoly;

    fn sub(self, other: &'b ZMPoly) -> ZMPoly {
        assert!(
            self.vars == other.vars,
            "Unmatch variables in binary operation"
        );
        let mut this = ZMPoly::new(&self.vars);
        unsafe {
            fmpz_mpoly_sub(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
                &other.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

///Implement addition with operator `*`
/// ```
/// use flint_mpoly::ZMPoly;
/// // Define new polynomial
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly_a = ZMPoly::from_str("x1",&vars).unwrap();
/// let mpoly_b = ZMPoly::from_str("x2",&vars).unwrap();
/// let mpoly_ab = &mpoly_a * &mpoly_b;
/// assert_eq!("+x1*x2", mpoly_ab.to_str());
/// ```
impl<'a, 'b> Mul<&'b ZMPoly> for &'a ZMPoly {
    type Output = ZMPoly;

    fn mul(self, other: &'b ZMPoly) -> ZMPoly {
        assert!(
            self.vars == other.vars,
            "Unmatch variables in binary operation"
        );
        let mut this = ZMPoly::new(&self.vars);
        unsafe {
            fmpz_mpoly_mul(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
                &other.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

///Implement negative sign
/// ```
/// use flint_mpoly::ZMPoly;
/// // Define new polynomial
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly = ZMPoly::from_str("x1",&vars).unwrap();
/// let mpoly_neg = -&mpoly;
/// assert_eq!("-x1", mpoly_neg.to_str());
/// ```
impl<'a> Neg for &'a ZMPoly {
    type Output = ZMPoly;

    fn neg(self) -> ZMPoly {
        let mut this = ZMPoly::new(&self.vars);
        unsafe {
            fmpz_mpoly_neg(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

/// Clone for ZMPoly
/// ```
/// use flint_mpoly::ZMPoly;
/// use std::str::FromStr;
/// let vars = [String::from("x1"),String::from("x2")];
/// let mpoly_1 = ZMPoly::from_str("x1+x2",&vars).unwrap();
/// let mpoly_2 = mpoly_1.clone();
/// assert_eq!(mpoly_1.to_str(),mpoly_2.to_str());
/// ```
impl Clone for ZMPoly {
    fn clone(&self) -> Self {
        ZMPoly::clone_from(self)
    }
}

/// Clear the content of all raw pointers before droping ZMPoly
impl Drop for ZMPoly {
    fn drop(&mut self) {
        unsafe {
            fmpz_mpoly_clear(&mut self.raw as *mut _, &mut self.ctx as *mut _);
            fmpz_mpoly_ctx_clear(&mut self.ctx as *mut _);
            //flint_cleanup();
        }
    }
}

// To use the `{}` for the structure ZMPoly
impl fmt::Display for ZMPoly {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        unsafe {
            let length = fmpz_mpoly_length(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );

            let nvars = self.vars.len();
            for ci in 0..length {
                // get coeff
                fmpz_mpoly_get_term_coeff_fmpz(
                    &self._coeff1.raw as *const _ as *mut _,
                    &self.raw as *const _ as *mut _,
                    ci,
                    &self.ctx as *const _ as *mut _,
                );
                // get exponent
                fmpz_mpoly_get_term_exp_ui(
                    &self._exp as *const _ as *mut _,
                    &self.raw as *const _ as *mut _,
                    ci,
                    &self.ctx as *const _ as *mut _,
                );

                if fmpz_is_zero(&self._coeff1.raw as *const _ as *mut _) == 0 {
                    write!(f, "{:+}", self._coeff1)?;

                    for (n, p) in self._exp[..nvars].iter().enumerate() {
                        match *p {
                            0 => continue,
                            1 => {
                                write!(f, "*{}", self.vars[n])?;
                            }
                            _ => {
                                write!(f, "*{}^{}", self.vars[n], p)?;
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

// Send and Sync
unsafe impl Send for ZMPoly {}
unsafe impl Sync for ZMPoly {}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    #[should_panic]
    fn bad_string_1() {
        let vars = [String::from("x1"), String::from("x2")];
        ZMPoly::from_str("x1/x2", &vars).unwrap();
    }
    #[test]
    #[should_panic]
    fn bad_string_2() {
        let vars = [String::from("x1")];
        ZMPoly::from_str("+x1/2", &vars).unwrap();
    }
}
