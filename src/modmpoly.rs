use crate::parser::parse_mod_mpoly;
use crate::{ModCoeff, ZCoeff};
use flint_sys::fmpz::*;
use flint_sys::fmpz_mod_mpoly::*;
use flint_sys::mpoly::ordering_t_ORD_DEGLEX;
use regex::Regex;
use std::fmt;
use std::mem::MaybeUninit;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Container for polynomials over rational numbers
#[derive(Debug)]
pub struct ModMPoly {
    pub raw: fmpz_mod_mpoly_struct,     // Polynomial
    pub ctx: fmpz_mod_mpoly_ctx_struct, // Polynomial context
    pub vars: Vec<String>,              // Name of variables
    pub modulo: u64,                    // mod n size
    _coeff1: ModCoeff,                  // coefficient container
    _coeff2: ModCoeff,                  // coefficient container
    _res: ModCoeff,                     // coefficient container
    _exp: [u64; 10],                    // exponent container
}

/// Container for polynomial that uses Flint as backend
impl ModMPoly {
    /// New 0 polynomial with selected variables
    /// ```
    /// use flint_mpoly::ModMPoly;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly = ModMPoly::new(&vars,&modulo);
    /// assert_eq!("0",mpoly.to_str());
    /// ```
    pub fn new(vars: &[String], modulo: &u64) -> Self {
        // Set modulo
        let n = ZCoeff::from_int((*modulo).try_into().expect("modulo not fitting in i64"));
        unsafe {
            let mut raw = MaybeUninit::uninit();
            let mut ctx = MaybeUninit::uninit();
            let nvars = vars.len() as i64;
            fmpz_mod_mpoly_ctx_init(
                ctx.as_mut_ptr(),
                nvars,
                ordering_t_ORD_DEGLEX,
                &n as *const _ as *mut _,
            );
            fmpz_mod_mpoly_init(raw.as_mut_ptr(), ctx.as_mut_ptr());
            ModMPoly {
                //raw: raw,
                raw: raw.assume_init(),
                ctx: ctx.assume_init(),
                vars: vars.to_vec(),
                modulo: *modulo,
                _coeff1: ModCoeff::new(modulo),
                _coeff2: ModCoeff::new(modulo),
                _res: ModCoeff::new(modulo),
                _exp: [0; 10],
            }
        }
    }

    /// Initialize polynomial from string assuming the given variables
    /// ```
    /// use flint_mpoly::ModMPoly;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly= ModMPoly::from_str("(x1+x2)", &vars, &modulo).unwrap();
    /// assert_eq!("+x1+x2",mpoly.to_str());
    /// ```
    pub fn from_str(expr_str: &str, vars: &[String], modulo: &u64) -> Result<Self, &'static str> {
        Ok(parse_mod_mpoly(expr_str, vars, modulo).expect("parse polynomial"))
    }

    /// Clone will copy the pointers and we will have two objects pointing at the same memory
    /// address. With this function we can create a clone in a safe way
    /// ```
    /// use flint_mpoly::ModMPoly;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly_1= ModMPoly::from_str("(x1+x2)", &vars, &modulo).unwrap();
    /// let mpoly_2= ModMPoly::clone_from(&mpoly_1);
    /// assert_eq!(mpoly_1.to_str(), mpoly_2.to_str());
    /// ```
    pub fn clone_from(other: &Self) -> Self {
        let mut this = ModMPoly::new(&other.vars, &other.modulo);
        unsafe {
            fmpz_mod_mpoly_set(
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
            fmpz_mod_mpoly_zero(&mut self.raw as *mut _, &mut self.ctx as *mut _);
        }
    }

    /// Set polynomial to one
    pub fn set_to_one(&mut self) {
        unsafe {
            fmpz_mod_mpoly_one(&mut self.raw as *mut _, &mut self.ctx as *mut _);
        }
    }

    /// Set polynomial from string assuming the variables already stored in ModMPoly
    /// ```
    /// use flint_mpoly::ModMPoly;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mpoly= ModMPoly::new(&vars,&modulo);
    /// mpoly.set_from_str("x1+x2").unwrap();
    /// assert_eq!("+x1+x2",mpoly.to_str());
    /// ```
    pub fn set_from_str(&mut self, expr_str: &str) -> Result<(), &'static str> {
        *self = ModMPoly::from_str(expr_str, &self.vars, &self.modulo)?;
        Ok(())
    }

    /// Check if the polynomial is zero
    pub fn is_zero(&self) -> bool {
        unsafe {
            fmpz_mod_mpoly_is_zero(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }

    /// Check if the polynomial is one
    pub fn is_one(&self) -> bool {
        unsafe {
            fmpz_mod_mpoly_is_one(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }
    /// Check if the polynomial is const
    pub fn is_const(&self) -> bool {
        unsafe {
            fmpz_mod_mpoly_is_fmpz(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }

    /// Add to the polynomial the coefficient with the specified powers
    /// ```
    /// use flint_mpoly::ModMPoly;
    /// use flint_mpoly::ModCoeff;
    /// // Define new polynomial
    /// let modulo = 11311;
    /// let vars = [String::from("x1"), String::from("x2")];
    /// let mut mpoly= ModMPoly::new(&vars, &modulo);
    /// // Set coefficient of monomial with power x1^2*x2^3 and coefficient 3/2
    /// let pows = [2,3];
    /// let coeff = ModCoeff::from_str("3/2", &modulo).unwrap();
    /// mpoly.add_coeff(&pows,&coeff);
    /// assert_eq!("+5657*x1^2*x2^3", mpoly.to_str());
    /// ```
    pub fn add_coeff(&mut self, monimial_powers: &[u64], coeff: &ModCoeff) {
        unsafe {
            fmpz_mod_mpoly_get_coeff_fmpz_ui(
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
            fmpz_mod_mpoly_set_coeff_fmpz_ui(
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
    /// use flint_mpoly::ModMPoly;
    /// // Define new polynomial
    /// let modulo = 11311;
    /// let vars = [String::from("x1"), String::from("x2")];
    /// let mut mpoly= ModMPoly::new(&vars,&modulo);
    /// // Set coefficient of monomial with power x1^2*x2^3 and coefficient 3
    /// let pows = [2,3];
    /// mpoly.add_coeff_int(&pows,3);
    /// assert_eq!("+3*x1^2*x2^3", mpoly.to_str());
    /// ```
    pub fn add_coeff_int(&mut self, monimial_powers: &[u64], n: i64) {
        unsafe {
            // Store current coefficient to _coeff1
            fmpz_mod_mpoly_get_coeff_fmpz_ui(
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
            fmpz_mod_mpoly_set_coeff_fmpz_ui(
                &mut self.raw as *mut _,
                &mut self._res.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
        }
    }

    /// Add to the polynomial the coefficient p/q with the specified powers
    /// ```
    /// use flint_mpoly::ModMPoly;
    /// // Define new polynomial
    /// let modulo = 11311;
    /// let vars = [String::from("x1"), String::from("x2")];
    /// let mut mpoly= ModMPoly::new(&vars,&modulo);
    /// // Set coefficient of monomial with power x1^2*x2^3 and coefficient 3/2
    /// let pows = [2,3];
    /// mpoly.add_coeff_str(&pows,"3/2");
    /// assert_eq!("+5657*x1^2*x2^3", mpoly.to_str());
    /// ```
    pub fn add_coeff_str(&mut self, monimial_powers: &[u64], coeff: &str) {
        unsafe {
            // Store current coefficient to _coeff1
            fmpz_mod_mpoly_get_coeff_fmpz_ui(
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
            fmpz_mod_mpoly_set_coeff_fmpz_ui(
                &mut self.raw as *mut _,
                &mut self._res.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
        }
    }
    /// Rise polynomial to integer power
    pub fn pown(&mut self, n: u64) {
        let mut pow_res = ModMPoly::new(&self.vars, &self.modulo);
        unsafe {
            fmpz_mod_mpoly_pow_ui(
                &mut pow_res.raw as *mut _,
                &mut self.raw as *mut _,
                n,
                &mut self.ctx as *mut _,
            );
            fmpz_mod_mpoly_swap(
                &mut pow_res.raw as *mut _,
                &mut self.raw as *mut _,
                &mut self.ctx as *mut _,
            );
        }
    }
    /// Swap the structure content with another ModMPoly provided they have the same context
    pub fn swap(&mut self, other: &mut Self) {
        assert!(
            self.vars == other.vars,
            "Unmatch variables while swapping expressions"
        );
        unsafe {
            fmpz_mod_mpoly_swap(
                &mut self.raw as *mut _,
                &mut other.raw as *mut _,
                &mut self.ctx as *mut _,
            );
        }
    }
    /// Set function to zero and return the previous value
    /// ```
    /// use flint_mpoly::ModMPoly;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mpoly_a = ModMPoly::from_str("+x1+x2",&vars,&modulo).unwrap();
    /// assert_eq!("+x1+x2",mpoly_a.to_str());
    /// let mut mpoly_b = mpoly_a.zero_move();
    /// assert_eq!("0",mpoly_a.to_str());
    /// assert_eq!("+x1+x2",mpoly_b.to_str());
    /// ```
    pub fn zero_move(&mut self) -> Self {
        let mut out = Self::new(&self.vars, &self.modulo);
        out.swap(self);
        out
    }
    /// Clear the function and set it to zero
    pub fn clear(&mut self) {
        self.set_to_zero();
    }
    /// Clear all memory allocated to ModMPoly
    pub fn free(&mut self) {
        self.clear();
        unsafe {
            fmpz_mod_mpoly_init(&mut self.raw as *mut _, &mut self.ctx as *mut _);
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
/// use flint_mpoly::ModMPoly;
/// // Define new polynomial
/// let modulo = 11311;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly_a = ModMPoly::from_str("x1",&vars,&modulo).unwrap();
/// let mpoly_b = ModMPoly::from_str("x2",&vars,&modulo).unwrap();
/// let mpoly_ab = &mpoly_a + &mpoly_b;
/// assert_eq!("+x1+x2", mpoly_ab.to_str());
/// ```
impl<'a, 'b> Add<&'b ModMPoly> for &'a ModMPoly {
    type Output = ModMPoly;

    fn add(self, other: &'b ModMPoly) -> ModMPoly {
        assert!(
            self.vars == other.vars,
            "Unmatch variables in binary operation"
        );
        let mut this = ModMPoly::new(&self.vars, &self.modulo);
        unsafe {
            fmpz_mod_mpoly_add(
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
/// use flint_mpoly::ModMPoly;
/// use flint_mpoly::parse_mpoly;
/// // Define new polynomial
/// let modulo = 11311;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly_a = ModMPoly::from_str("x1",&vars,&modulo).unwrap();
/// let mpoly_b = ModMPoly::from_str("x2",&vars,&modulo).unwrap();
/// let mpoly_ab = &mpoly_a - &mpoly_b;
/// assert_eq!("+x1+11310*x2", mpoly_ab.to_str());
/// ```
impl<'a, 'b> Sub<&'b ModMPoly> for &'a ModMPoly {
    type Output = ModMPoly;

    fn sub(self, other: &'b ModMPoly) -> ModMPoly {
        assert!(
            self.vars == other.vars,
            "Unmatch variables in binary operation"
        );
        let mut this = ModMPoly::new(&self.vars, &self.modulo);
        unsafe {
            fmpz_mod_mpoly_sub(
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
/// use flint_mpoly::ModMPoly;
/// // Define new polynomial
/// let modulo = 11311;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly_a = ModMPoly::from_str("x1",&vars,&modulo).unwrap();
/// let mpoly_b = ModMPoly::from_str("x2",&vars,&modulo).unwrap();
/// let mpoly_ab = &mpoly_a * &mpoly_b;
/// assert_eq!("+x1*x2", mpoly_ab.to_str());
/// ```
impl<'a, 'b> Mul<&'b ModMPoly> for &'a ModMPoly {
    type Output = ModMPoly;

    fn mul(self, other: &'b ModMPoly) -> ModMPoly {
        assert!(
            self.vars == other.vars,
            "Unmatch variables in binary operation"
        );
        let mut this = ModMPoly::new(&self.vars, &self.modulo);
        unsafe {
            fmpz_mod_mpoly_mul(
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
/// use flint_mpoly::ModMPoly;
/// // Define new polynomial
/// let modulo = 11311;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mpoly = ModMPoly::from_str("x1",&vars,&modulo).unwrap();
/// let mpoly_neg = -&mpoly;
/// assert_eq!("+11310*x1", mpoly_neg.to_str());
/// ```
impl<'a> Neg for &'a ModMPoly {
    type Output = ModMPoly;

    fn neg(self) -> ModMPoly {
        let mut this = ModMPoly::new(&self.vars, &self.modulo);
        unsafe {
            fmpz_mod_mpoly_neg(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

/// Clone for ModMPoly
/// ```
/// use flint_mpoly::ModMPoly;
/// use std::str::FromStr;
/// let modulo = 11311;
/// let vars = [String::from("x1"),String::from("x2")];
/// let mpoly_1 = ModMPoly::from_str("x1+x2",&vars,&modulo).unwrap();
/// let mpoly_2 = mpoly_1.clone();
/// assert_eq!(mpoly_1.to_str(),mpoly_2.to_str());
/// ```
impl Clone for ModMPoly {
    fn clone(&self) -> Self {
        ModMPoly::clone_from(self)
    }
}

/// Clear the content of all raw pointers before droping ModMPoly
impl Drop for ModMPoly {
    fn drop(&mut self) {
        unsafe {
            fmpz_mod_mpoly_clear(&mut self.raw as *mut _, &mut self.ctx as *mut _);
            fmpz_mod_mpoly_ctx_clear(&mut self.ctx as *mut _);
        }
    }
}

// To use the `{}` for the structure ModMPoly
impl fmt::Display for ModMPoly {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        unsafe {
            let length = fmpz_mod_mpoly_length(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );

            let nvars = self.vars.len();
            for ci in 0..length {
                // get coeff
                fmpz_mod_mpoly_get_term_coeff_fmpz(
                    &self._coeff1.raw as *const _ as *mut _,
                    &self.raw as *const _ as *mut _,
                    ci,
                    &self.ctx as *const _ as *mut _,
                );
                // get exponent
                fmpz_mod_mpoly_get_term_exp_ui(
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
unsafe impl Send for ModMPoly {}
unsafe impl Sync for ModMPoly {}
