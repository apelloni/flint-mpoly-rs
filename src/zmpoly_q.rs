use crate::parser::parse_mrat;
use crate::{ZCoeff, ZMPoly};
use flint_sys::flint::*;
use flint_sys::fmpq::*;
use flint_sys::fmpq_mpoly::*;
use flint_sys::fmpz::*;
use flint_sys::fmpz_mpoly::*;
use flint_sys::fmpz_mpoly_q::*;
use flint_sys::mpoly::ordering_t_ORD_DEGLEX;
use regex::Regex;
use std::fmt;
use std::mem::MaybeUninit;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Container for polynomials over rational numbers
#[derive(Debug)]
pub struct ZMPolyQ {
    pub raw: fmpz_mpoly_q_struct,   // Rational function
    pub ctx: fmpz_mpoly_ctx_struct, // Polynomial context
    pub vars: Vec<String>,          // Name of variables
    _coeff1: ZCoeff,                // coefficient container
    _coeff2: ZCoeff,                // coefficient container
    _res: ZCoeff,                   // coefficient container
    _exp: [u64; 10],                // exponent container
}

/// Container for polynomial that uses Flint as backend
impl ZMPolyQ {
    /// New 0 polynomial with selected variables
    /// ```
    /// use flint_mpoly::ZMPolyQ;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly_q  = ZMPolyQ::new(&vars);
    /// assert_eq!("0",mpoly_q.to_str());
    /// ```
    pub fn new(vars: &[String]) -> Self {
        unsafe {
            let mut raw = MaybeUninit::uninit();
            let mut ctx = MaybeUninit::uninit();
            let nvars = vars.len() as i64;
            fmpz_mpoly_ctx_init(ctx.as_mut_ptr(), nvars, ordering_t_ORD_DEGLEX);
            fmpz_mpoly_q_init(raw.as_mut_ptr(), ctx.as_mut_ptr());
            ZMPolyQ {
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
    /// use flint_mpoly::ZMPolyQ;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly_q = ZMPolyQ::from_str("(2-x1+3*x2)/(x2-1+x2^2)", &vars).unwrap();
    /// assert_eq!("(-x1+3*x2+2)/(+x2^2+x2-1)",mpoly_q.to_str());
    /// ```
    pub fn from_str(expr_str: &str, vars: &[String]) -> Result<Self, &'static str> {
        let tmp = parse_mrat(expr_str, vars)
            .expect("parse polynomial")
            .to_zmrat()?;
        let mut expr = ZMPolyQ::new(vars);
        unsafe {
            fmpz_mpoly_set(
                &mut expr.raw.num as *mut _,
                &tmp.num.raw as *const _ as *mut _,
                &expr.ctx as *const _ as *mut _,
            );
            fmpz_mpoly_set(
                &mut expr.raw.den as *mut _,
                &tmp.den.raw as *const _ as *mut _,
                &expr.ctx as *const _ as *mut _,
            );
        }
        Ok(expr)
    }

    /// Clone will copy the pointers and we will have two objects pointing at the same memory
    /// address. With this function we can create a clone in a safe way
    /// ```
    /// use flint_mpoly::ZMPolyQ;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly_q_1= ZMPolyQ::from_str("(x1+x2)", &vars).unwrap();
    /// let mpoly_q_2= ZMPolyQ::clone_from(&mpoly_q_1);
    /// assert_eq!(mpoly_q_1.to_str(), mpoly_q_2.to_str());
    /// ```
    pub fn clone_from(other: &Self) -> Self {
        let mut this = ZMPolyQ::new(&other.vars);
        unsafe {
            fmpz_mpoly_q_set(
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
            fmpz_mpoly_q_zero(&mut self.raw as *mut _, &mut self.ctx as *mut _);
        }
    }

    /// Set polynomial to one
    pub fn set_to_one(&mut self) {
        unsafe {
            fmpz_mpoly_q_one(&mut self.raw as *mut _, &mut self.ctx as *mut _);
        }
    }

    /// Set polynomial from string assuming the variables already stored in ZMPolyQ
    /// ```
    /// use flint_mpoly::ZMPolyQ;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mpoly_q = ZMPolyQ::new(&vars);
    /// mpoly_q.set_from_str("3*x1-5*x2").unwrap();
    /// assert_eq!("+3*x1-5*x2",mpoly_q.to_str());
    /// ```
    pub fn set_from_str(&mut self, expr_str: &str) -> Result<(), &'static str> {
        *self = ZMPolyQ::from_str(expr_str, &self.vars)?;
        Ok(())
    }

    /// Check if the polynomial is zero
    pub fn is_zero(&self) -> bool {
        unsafe {
            fmpz_mpoly_q_is_zero(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }

    /// Check if the polynomial is one
    pub fn is_one(&self) -> bool {
        unsafe {
            fmpz_mpoly_q_is_one(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }
    /// Check if the polynomial is const
    pub fn is_const(&self) -> bool {
        unsafe {
            fmpz_mpoly_q_is_fmpq(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }

    /// Rise polynomial to integer power
    pub fn pown(&mut self, n: u64) {
        let mut pow_res = ZMPoly::new(&self.vars);
        unsafe {
            fmpz_mpoly_pow_ui(
                &mut pow_res.raw as *mut _,
                &mut self.raw.num as *mut _,
                n,
                &mut self.ctx as *mut _,
            );
            fmpz_mpoly_swap(
                &mut pow_res.raw as *mut _,
                &mut self.raw.num as *mut _,
                &mut self.ctx as *mut _,
            );
            fmpz_mpoly_pow_ui(
                &mut pow_res.raw as *mut _,
                &mut self.raw.den as *mut _,
                n,
                &mut self.ctx as *mut _,
            );
            fmpz_mpoly_swap(
                &mut pow_res.raw as *mut _,
                &mut self.raw.den as *mut _,
                &mut self.ctx as *mut _,
            );
        }
    }
    /// Swap the structure content with another QMPoly provided they have the same context
    pub fn swap(&mut self, other: &mut Self) {
        assert!(
            self.vars == other.vars,
            "Unmatch variables while swapping expressions"
        );
        unsafe {
            fmpz_mpoly_swap(
                &mut self.raw.num as *mut _,
                &mut other.raw.num as *mut _,
                &mut self.ctx as *mut _,
            );
            fmpz_mpoly_swap(
                &mut self.raw.den as *mut _,
                &mut other.raw.den as *mut _,
                &mut self.ctx as *mut _,
            );
        }
    }
    /// Set function to zero and return the previous value
    /// ```
    /// use flint_mpoly::ZMPolyQ;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mpoly_a = ZMPolyQ::from_str("+x1+x2",&vars).unwrap();
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
    /// Clear all memory allocated to QMPoly
    pub fn free(&mut self) {
        self.clear();
        unsafe {
            fmpz_mpoly_realloc(&mut self.raw.num as *mut _, 0, &mut self.ctx as *mut _);
            fmpz_mpoly_realloc(&mut self.raw.den as *mut _, 0, &mut self.ctx as *mut _);
        }
    }

    /// Convert the polynomial to a human readable string
    pub fn to_str(&self) -> String {
        let rg = Regex::new(r"([+-])1\*").unwrap();
        rg.replace_all(format!("{}", self).as_str(), "$1")
            .to_string()
    }
}

///Implement addition with operator `+`
/// ```
/// use flint_mpoly::ZMPolyQ;
/// // Define new polynomial
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly_a = ZMPolyQ::from_str("x1",&vars).unwrap();
/// let mpoly_b = ZMPolyQ::from_str("x2",&vars).unwrap();
/// let mpoly_ab = &mpoly_a + &mpoly_b;
/// assert_eq!("+x1+x2", mpoly_ab.to_str());
/// ```
impl<'a, 'b> Add<&'b ZMPolyQ> for &'a ZMPolyQ {
    type Output = ZMPolyQ;

    fn add(self, other: &'b ZMPolyQ) -> ZMPolyQ {
        assert!(
            self.vars == other.vars,
            "Unmatch variables in binary operation"
        );
        let mut this = ZMPolyQ::new(&self.vars);
        unsafe {
            fmpz_mpoly_q_add(
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
/// use flint_mpoly::ZMPolyQ;
/// use flint_mpoly::parse_mpoly;
/// // Define new polynomial
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly_a = ZMPolyQ::from_str("x1",&vars).unwrap();
/// let mpoly_b = ZMPolyQ::from_str("x2",&vars).unwrap();
/// let mpoly_ab = &mpoly_a - &mpoly_b;
/// assert_eq!("+x1-x2", mpoly_ab.to_str());
/// ```
impl<'a, 'b> Sub<&'b ZMPolyQ> for &'a ZMPolyQ {
    type Output = ZMPolyQ;

    fn sub(self, other: &'b ZMPolyQ) -> ZMPolyQ {
        assert!(
            self.vars == other.vars,
            "Unmatch variables in binary operation"
        );
        let mut this = ZMPolyQ::new(&self.vars);
        unsafe {
            fmpz_mpoly_q_sub(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
                &other.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

/// Set function to zero and return the previous value
/// ```
/// use flint_mpoly::ZMPolyQ;
/// let vars = [String::from("x1"),String::from("x2")];
/// let mut mpoly_a = ZMPolyQ::from_str("(x1+x2)/(x2-1)",&vars).unwrap();
/// let mut mpoly_b = ZMPolyQ::from_str("(x1*x2-x1)/(x2+1)",&vars).unwrap();
/// let mpoly_res = &mpoly_a * &mpoly_b;
/// assert_eq!("(+x1^2+x1*x2)/(+x2+1)",mpoly_res.to_str());
/// ```

impl<'a, 'b> Mul<&'b ZMPolyQ> for &'a ZMPolyQ {
    type Output = ZMPolyQ;

    fn mul(self, other: &'b ZMPolyQ) -> ZMPolyQ {
        assert!(
            self.vars == other.vars,
            "Unmatch variables in binary operation"
        );
        let mut this = ZMPolyQ::new(&self.vars);
        unsafe {
            fmpz_mpoly_q_mul(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
                &other.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

/// Set function to zero and return the previous value
/// ```
/// use flint_mpoly::ZMPolyQ;
/// let vars = [String::from("x1"),String::from("x2")];
/// let mut mpoly_a = ZMPolyQ::from_str("(x1+x2)/(x2-1)",&vars).unwrap();
/// let mut mpoly_b = ZMPolyQ::from_str("(x1*x2-x1)/(x2+1)",&vars).unwrap();
/// let mpoly_res = &mpoly_a / &mpoly_b;
/// assert_eq!("(+x1*x2+x2^2+x1+x2)/(+x1*x2^2-2*x1*x2+x1)",mpoly_res.to_str());
/// ```
impl<'a, 'b> Div<&'b ZMPolyQ> for &'a ZMPolyQ {
    type Output = ZMPolyQ;

    fn div(self, other: &'b ZMPolyQ) -> ZMPolyQ {
        assert!(
            self.vars == other.vars,
            "Unmatch variables in binary operation"
        );
        let mut this = ZMPolyQ::new(&self.vars);
        unsafe {
            fmpz_mpoly_q_div(
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
/// use flint_mpoly::ZMPolyQ;
/// // Define new polynomial
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly = ZMPolyQ::from_str("x1",&vars).unwrap();
/// let mpoly_neg = -&mpoly;
/// assert_eq!("-x1", mpoly_neg.to_str());
/// ```
impl<'a> Neg for &'a ZMPolyQ {
    type Output = ZMPolyQ;

    fn neg(self) -> ZMPolyQ {
        let mut this = ZMPolyQ::new(&self.vars);
        unsafe {
            fmpz_mpoly_q_neg(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

/// Clone for ZMPolyQ
/// ```
/// use flint_mpoly::ZMPolyQ;
/// use std::str::FromStr;
/// let vars = [String::from("x1"),String::from("x2")];
/// let mpoly_1 = ZMPolyQ::from_str("x1+x2",&vars).unwrap();
/// let mpoly_2 = mpoly_1.clone();
/// assert_eq!(mpoly_1.to_str(),mpoly_2.to_str());
/// ```
impl Clone for ZMPolyQ {
    fn clone(&self) -> Self {
        ZMPolyQ::clone_from(self)
    }
}

/// Clear the content of all raw pointers before droping ZMPolyQ
impl Drop for ZMPolyQ {
    fn drop(&mut self) {
        unsafe {
            fmpz_mpoly_q_clear(&mut self.raw as *mut _, &mut self.ctx as *mut _);
            fmpz_mpoly_ctx_clear(&mut self.ctx as *mut _);
            //flint_cleanup();
        }
    }
}

// To use the `{}` for the structure ZMPolyQ
impl fmt::Display for ZMPolyQ {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }

        let den_is_one = unsafe {
            fmpz_mpoly_is_one(
                &self.raw.den as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        };

        if !den_is_one {
            write!(f, "(")?;
        }

        for (i, &num_den) in [self.raw.num, self.raw.den].iter().enumerate() {
            unsafe {
                let length = fmpz_mpoly_length(
                    &num_den as *const _ as *mut _,
                    &self.ctx as *const _ as *mut _,
                );
                let nvars = self.vars.len();
                for ci in 0..length {
                    // get coeff
                    fmpz_mpoly_get_term_coeff_fmpz(
                        &self._coeff1.raw as *const _ as *mut _,
                        &num_den as *const _ as *mut _,
                        ci,
                        &self.ctx as *const _ as *mut _,
                    );
                    // get exponent
                    fmpz_mpoly_get_term_exp_ui(
                        &self._exp as *const _ as *mut _,
                        &num_den as *const _ as *mut _,
                        ci,
                        &self.ctx as *const _ as *mut _,
                    );

                    if fmpq_is_zero(&self._coeff1.raw as *const _ as *mut _) == 0 {
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
            if i == 0 && !den_is_one {
                write!(f, ")/(")?;
            } else {
                break;
            }
        }
        if !den_is_one {
            write!(f, ")")?;
        }
        Ok(())
    }
}

// Send and Sync
unsafe impl Send for ZMPolyQ {}
unsafe impl Sync for ZMPolyQ {}
