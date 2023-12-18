use crate::parser::parse_mpoly;
use crate::qcoeff::QCoeff;
use flint_sys::fmpq::*;
use flint_sys::fmpq_mpoly::*;
use flint_sys::mpoly::ordering_t_ORD_DEGLEX;
use regex::Regex;
use std::fmt;
use std::mem::MaybeUninit;
use std::ops::{Add, Div, Mul, Neg, Sub};

pub struct QMPoly {
    pub raw: fmpq_mpoly_struct,     // Polynomial
    pub ctx: fmpq_mpoly_ctx_struct, // Polynomial context
    pub vars: Vec<String>,          // Name of variables
    _coeff1: QCoeff,                // coefficient container
    _coeff2: QCoeff,                // coefficient container
    _res: QCoeff,                   // coefficient container
    _exp: [u64; 10],                // exponent container
}

/// Container for polynomial that uses Flint as backend
impl QMPoly {
    /// New 0 polynomial with selected variables
    /// ```
    /// use flint_mpoly::qmpoly::QMPoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly = QMPoly::new(&vars);
    /// assert_eq!("0",mpoly.to_str());
    /// ```
    pub fn new(vars: &[String]) -> Self {
        unsafe {
            let mut raw = MaybeUninit::uninit();
            let mut ctx = MaybeUninit::uninit();
            let nvars = vars.len() as i64;
            fmpq_mpoly_ctx_init(ctx.as_mut_ptr(), nvars, ordering_t_ORD_DEGLEX);
            fmpq_mpoly_init(raw.as_mut_ptr(), ctx.as_mut_ptr());
            QMPoly {
                //raw: raw,
                raw: raw.assume_init(),
                ctx: ctx.assume_init(),
                vars: vars.to_vec(),
                _coeff1: QCoeff::default(),
                _coeff2: QCoeff::default(),
                _res: QCoeff::default(),
                _exp: [0; 10],
            }
        }
    }

    /// Initialize polynomial from string assuming the given variables
    /// ```
    /// use flint_mpoly::qmpoly::QMPoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly= QMPoly::from_str("(x1+x2)", &vars).unwrap();
    /// assert_eq!("+x1+x2",mpoly.to_str());
    /// ```
    pub fn from_str(expr_str: &str, vars: &[String]) -> Result<Self, &'static str> {
        Ok(parse_mpoly(expr_str, vars).expect("parse polynomial"))
    }

    /// Clone will copy the pointers and we will have two objects pointing at the same memory
    /// address. With this function we can create a clone in a safe way
    /// ```
    /// use flint_mpoly::qmpoly::QMPoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly_1= QMPoly::from_str("(x1+x2)", &vars).unwrap();
    /// let mpoly_2= QMPoly::clone_from(&mpoly_1);
    /// assert_eq!(mpoly_1.to_str(), mpoly_2.to_str());
    /// ```
    pub fn clone_from(other: &Self) -> Self {
        let mut this = QMPoly::new(&other.vars);
        unsafe {
            fmpq_mpoly_set(
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
            fmpq_mpoly_zero(&mut self.raw as *mut _, &mut self.ctx as *mut _);
        }
    }

    /// Set polynomial to one
    pub fn set_to_one(&mut self) {
        unsafe {
            fmpq_mpoly_one(&mut self.raw as *mut _, &mut self.ctx as *mut _);
        }
    }

    /// Set polynomial from string assuming the variables already stored in QMPoly
    /// ```
    /// use flint_mpoly::qmpoly::QMPoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mpoly= QMPoly::new(&vars);
    /// mpoly.set_from_str("x1+x2").unwrap();
    /// assert_eq!("+x1+x2",mpoly.to_str());
    /// ```
    pub fn set_from_str(&mut self, expr_str: &str) -> Result<(), &'static str> {
        *self = QMPoly::from_str(expr_str, &self.vars)?;
        Ok(())
    }

    /// Check if the polynomial is zero
    pub fn is_zero(&self) -> bool {
        unsafe {
            fmpq_mpoly_is_zero(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }

    /// Check if the polynomial is one
    pub fn is_one(&self) -> bool {
        unsafe {
            fmpq_mpoly_is_one(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }
    /// Check if the polynomial is const
    pub fn is_const(&self) -> bool {
        unsafe {
            fmpq_mpoly_is_fmpq(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }

    /// Add to the polynomial the coefficient with the specified powers
    /// ```
    /// use flint_mpoly::qmpoly::QMPoly;
    /// use flint_mpoly::qcoeff::QCoeff;
    /// // Define new polynomial
    /// let vars = [String::from("x1"), String::from("x2")];
    /// let mut mpoly= QMPoly::new(&vars);
    /// // Set coefficient of monomial with power x1^2*x2^3 and coefficient 3/2
    /// let pows = [2,3];
    /// let coeff = QCoeff::from_int(3, 2);
    /// mpoly.add_coeff(&pows,&coeff);
    /// assert_eq!("+3/2*x1^2*x2^3", mpoly.to_str());
    /// ```
    pub fn add_coeff(&mut self, monimial_powers: &[u64], coeff: &QCoeff) {
        unsafe {
            fmpq_mpoly_get_coeff_fmpq_ui(
                &mut self._coeff1.raw as *mut _,
                &mut self.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
            if self._coeff1.is_zero() {
                fmpq_set(
                    &mut self._res.raw as *mut _,
                    &coeff.raw as *const _ as *mut _,
                )
            } else {
                fmpq_add(
                    &mut self._res.raw as *mut _,
                    &mut self._coeff1.raw as *mut _,
                    &coeff.raw as *const _ as *mut _,
                );
            }
            fmpq_mpoly_set_coeff_fmpq_ui(
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
    /// use flint_mpoly::qmpoly::QMPoly;
    /// // Define new polynomial
    /// let vars = [String::from("x1"), String::from("x2")];
    /// let mut mpoly= QMPoly::new(&vars);
    /// // Set coefficient of monomial with power x1^2*x2^3 and coefficient 3/2
    /// let pows = [2,3];
    /// mpoly.add_coeff_pq(&pows,3,2);
    /// assert_eq!("+3/2*x1^2*x2^3", mpoly.to_str());
    /// ```
    pub fn add_coeff_pq(&mut self, monimial_powers: &[u64], p: i64, q: u64) {
        unsafe {
            // Store current coefficient to _coeff1
            fmpq_mpoly_get_coeff_fmpq_ui(
                &mut self._coeff1.raw as *mut _,
                &mut self.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
            // Set coefficient to add to _coeff2
            self._coeff2.set_from_int(p, q);
            if self._coeff1.is_zero() {
                fmpq_swap(
                    &mut self._res.raw as *mut _,
                    &mut self._coeff2.raw as *mut _,
                )
            } else {
                fmpq_add(
                    &mut self._res.raw as *mut _,
                    &mut self._coeff1.raw as *mut _,
                    &mut self._coeff2.raw as *mut _,
                );
            }
            fmpq_mpoly_set_coeff_fmpq_ui(
                &mut self.raw as *mut _,
                &mut self._res.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
        }
    }

    /// Add to the polynomial the coefficient p/q with the specified powers
    /// ```
    /// use flint_mpoly::qmpoly::QMPoly;
    /// // Define new polynomial
    /// let vars = [String::from("x1"), String::from("x2")];
    /// let mut mpoly= QMPoly::new(&vars);
    /// // Set coefficient of monomial with power x1^2*x2^3 and coefficient 3/2
    /// let pows = [2,3];
    /// mpoly.add_coeff_str(&pows,"3/2");
    /// assert_eq!("+3/2*x1^2*x2^3", mpoly.to_str());
    /// ```
    pub fn add_coeff_str(&mut self, monimial_powers: &[u64], coeff: &str) {
        unsafe {
            // Store current coefficient to _coeff1
            fmpq_mpoly_get_coeff_fmpq_ui(
                &mut self._coeff1.raw as *mut _,
                &mut self.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
            // Set coefficient to add to _coeff2
            self._coeff2.set_from_str(coeff).unwrap();
            if self._coeff1.is_zero() {
                fmpq_swap(
                    &mut self._res.raw as *mut _,
                    &mut self._coeff2.raw as *mut _,
                );
            } else {
                fmpq_add(
                    &mut self._res.raw as *mut _,
                    &mut self._coeff1.raw as *mut _,
                    &mut self._coeff2.raw as *mut _,
                );
            }
            fmpq_mpoly_set_coeff_fmpq_ui(
                &mut self.raw as *mut _,
                &mut self._res.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
        }
    }
    /// Rise polynomial to integer power
    pub fn pown(&mut self, n: u64) {
        let mut pow_res = QMPoly::new(&self.vars);
        unsafe {
            fmpq_mpoly_pow_ui(
                &mut pow_res.raw as *mut _,
                &mut self.raw as *mut _,
                n,
                &mut self.ctx as *mut _,
            );
            fmpq_mpoly_swap(
                &mut pow_res.raw as *mut _,
                &mut self.raw as *mut _,
                &mut self.ctx as *mut _,
            );
        }
    }
    /// Clear the function and set it to zero
    pub fn clear(&mut self) {
        self.set_to_zero();
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
/// use flint_mpoly::qmpoly::QMPoly;
/// // Define new polynomial
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly_a = QMPoly::from_str("x1",&vars).unwrap();
/// let mpoly_b = QMPoly::from_str("x2",&vars).unwrap();
/// let mpoly_ab = &mpoly_a + &mpoly_b;
/// assert_eq!("+x1+x2", mpoly_ab.to_str());
/// ```
impl<'a, 'b> Add<&'b QMPoly> for &'a QMPoly {
    type Output = QMPoly;

    fn add(self, other: &'b QMPoly) -> QMPoly {
        // TODO: Check if the two contexts are the same
        let mut this = QMPoly::new(&self.vars);
        unsafe {
            fmpq_mpoly_add(
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
/// use flint_mpoly::qmpoly::QMPoly;
/// use flint_mpoly::parser::parse_mpoly;
/// // Define new polynomial
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly_a = QMPoly::from_str("x1",&vars).unwrap();
/// let mpoly_b = QMPoly::from_str("x2",&vars).unwrap();
/// let mpoly_ab = &mpoly_a - &mpoly_b;
/// assert_eq!("+x1-x2", mpoly_ab.to_str());
/// ```
impl<'a, 'b> Sub<&'b QMPoly> for &'a QMPoly {
    type Output = QMPoly;

    fn sub(self, other: &'b QMPoly) -> QMPoly {
        // TODO: Check if the two contexts are the same
        let mut this = QMPoly::new(&self.vars);
        unsafe {
            fmpq_mpoly_sub(
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
/// use flint_mpoly::qmpoly::QMPoly;
/// // Define new polynomial
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly_a = QMPoly::from_str("x1",&vars).unwrap();
/// let mpoly_b = QMPoly::from_str("x2",&vars).unwrap();
/// let mpoly_ab = &mpoly_a * &mpoly_b;
/// assert_eq!("+x1*x2", mpoly_ab.to_str());
/// ```
impl<'a, 'b> Mul<&'b QMPoly> for &'a QMPoly {
    type Output = QMPoly;

    fn mul(self, other: &'b QMPoly) -> QMPoly {
        // TODO: Check if the two contexts are the same
        let mut this = QMPoly::new(&self.vars);
        unsafe {
            fmpq_mpoly_mul(
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
/// use flint_mpoly::qmpoly::QMPoly;
/// // Define new polynomial
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly = QMPoly::from_str("x1",&vars).unwrap();
/// let mpoly_neg = -&mpoly;
/// assert_eq!("-x1", mpoly_neg.to_str());
/// ```
impl<'a> Neg for &'a QMPoly {
    type Output = QMPoly;

    fn neg(self) -> QMPoly {
        let mut this = QMPoly::new(&self.vars);
        unsafe {
            fmpq_mpoly_neg(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

/// Clone for QMPoly
/// ```
/// use flint_mpoly::qmpoly::QMPoly;
/// use std::str::FromStr;
/// let vars = [String::from("x1"),String::from("x2")];
/// let mpoly_1 = QMPoly::from_str("x1+x2",&vars).unwrap();
/// let mpoly_2 = mpoly_1.clone();
/// assert_eq!(mpoly_1.to_str(),mpoly_2.to_str());
/// ```
impl Clone for QMPoly {
    fn clone(&self) -> Self {
        QMPoly::clone_from(self)
    }
}

/// Clear the content of all raw pointers before droping QMPoly
impl Drop for QMPoly {
    fn drop(&mut self) {
        unsafe {
            fmpq_mpoly_clear(&mut self.raw as *mut _, &mut self.ctx as *mut _);
            fmpq_mpoly_ctx_clear(&mut self.ctx as *mut _);
        }
    }
}

// To use the `{}` for the structure QMPoly
impl fmt::Display for QMPoly {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        unsafe {
            let length = fmpq_mpoly_length(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );

            let nvars = self.vars.len();
            for ci in 0..length {
                // get coeff
                fmpq_mpoly_get_term_coeff_fmpq(
                    &self._coeff1.raw as *const _ as *mut _,
                    &self.raw as *const _ as *mut _,
                    ci,
                    &self.ctx as *const _ as *mut _,
                );
                // get exponent
                fmpq_mpoly_get_term_exp_ui(
                    &self._exp as *const _ as *mut _,
                    &self.raw as *const _ as *mut _,
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
        Ok(())
    }
}
