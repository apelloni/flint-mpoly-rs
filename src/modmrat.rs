use crate::parse_mod_mrat;
use crate::ModCoeff;
use crate::ModMPoly;
use crate::flint_sys::fmpz_mod_mpoly::*;
use regex::Regex;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Container for rational polynomials over rational numbers
#[derive(Debug)]
pub struct ModMRat {
    pub num: ModMPoly,     // numerator
    pub den: ModMPoly,     // Polynomial context
    pub vars: Vec<String>, // Name of variables
    pub modulo: u64,       // mod n size
    _reduced: bool,
    _coeff: ModCoeff, // coefficient container
    _gcd: ModMPoly,   // gcd container
    _res: ModMPoly,   // partial operations container
}

impl ModMRat {
    /// Create a rational function equal to 0
    /// ```
    /// use flint_mpoly::ModMRat;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mrat = ModMRat::new(&vars,&modulo);
    /// assert_eq!("0",mrat.to_str());
    /// ```
    pub fn new(vars: &[String], modulo: &u64) -> Self {
        let num = ModMPoly::new(vars, modulo);
        let mut den = ModMPoly::new(vars, modulo);
        den.set_to_one();
        ModMRat {
            num,
            den,
            vars: vars.to_vec(),
            modulo: *modulo,
            _reduced: false,
            _coeff: ModCoeff::new(modulo),
            _gcd: ModMPoly::new(vars, modulo),
            _res: ModMPoly::new(vars, modulo),
        }
    }

    /// Create a rational function from a polynomial
    /// ```
    /// use flint_mpoly::ModMPoly;
    /// use flint_mpoly::ModMRat;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly= ModMPoly::from_str("x1+3/4 x2",&vars,&modulo).unwrap();
    /// let mrat = ModMRat::from_mpoly(&mpoly,&modulo);
    /// assert_eq!("+x1+8484*x2",mrat.to_str());
    /// ```
    pub fn from_mpoly(mpoly: &ModMPoly, modulo: &u64) -> Self {
        let mut den = ModMPoly::new(&mpoly.vars, modulo);
        den.set_to_one();
        ModMRat {
            num: ModMPoly::clone_from(mpoly),
            den,
            vars: mpoly.vars.clone(),
            modulo: *modulo,
            _reduced: false,
            _coeff: ModCoeff::new(modulo),
            _gcd: ModMPoly::new(&mpoly.vars, modulo),
            _res: ModMPoly::new(&mpoly.vars, modulo),
        }
    }

    /// Create a rational function from two polynomials defining the numerator and denominator,
    /// respectively
    /// ```
    /// use flint_mpoly::ModMPoly;
    /// use flint_mpoly::ModMRat;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let num= ModMPoly::from_str("x1+3/4 x2",&vars,&modulo).unwrap();
    /// let den= ModMPoly::from_str("x1-x2",&vars,&modulo).unwrap();
    /// let mrat = ModMRat::from_mpolys(&num,&den);
    /// assert_eq!("(+x1+8484*x2)/(+x1+11310*x2)",mrat.to_str());
    /// ```
    pub fn from_mpolys(num: &ModMPoly, den: &ModMPoly) -> Self {
        assert_eq!(
            num.modulo, den.modulo,
            "Different modulo when calling from_mpolys"
        );
        assert_eq!(
            num.vars, den.vars,
            "Different variables when calling from_mpolys"
        );
        ModMRat {
            num: ModMPoly::clone_from(num),
            den: ModMPoly::clone_from(den),
            vars: num.vars.clone(),
            modulo: num.modulo,
            _reduced: false,
            _coeff: ModCoeff::new(&num.modulo),
            _gcd: ModMPoly::new(&num.vars, &num.modulo),
            _res: ModMPoly::new(&num.vars, &num.modulo),
        }
    }

    /// Initialize rational from string assuming the given variables
    /// ```
    /// use flint_mpoly::ModMRat;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mrat = ModMRat::from_str("(-1x1+x2)/x1", &vars,&modulo).unwrap();
    /// assert_eq!("(+11310*x1+x2)/(+x1)",mrat.to_str());
    /// ```
    pub fn from_str(expr_str: &str, vars: &[String], modulo: &u64) -> Result<Self, &'static str> {
        Ok(parse_mod_mrat(expr_str, vars, modulo).expect("parse rational"))
    }

    /// Set rational from string assuming the variables already stored in ModMRat
    /// ```
    /// use flint_mpoly::ModMRat;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat = ModMRat::new(&vars,&modulo);
    /// mrat.set_from_str("(x1+x2)/x1").unwrap();
    /// assert_eq!("(+x1+x2)/(+x1)",mrat.to_str());
    /// ```
    pub fn set_from_str(&mut self, expr_str: &str) -> Result<(), &'static str> {
        *self = ModMRat::from_str(expr_str, &self.vars, &self.modulo)?;
        Ok(())
    }

    /// Clone will copy the pointers and we will have two objects pointing at the same memory
    /// address. With this function we can create a clone in a safe way
    /// ```
    /// use flint_mpoly::ModMRat;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mrat_1 = ModMRat::from_str("x1/(x2+1)",&vars,&modulo).unwrap();
    /// let mrat_2 = ModMRat::clone_from(&mrat_1);
    /// assert_eq!(mrat_1.to_str(),mrat_2.to_str());
    /// ```
    pub fn clone_from(other: &Self) -> Self {
        let mut this = ModMRat::new(&other.vars, &other.modulo);
        unsafe {
            // copy numerator
            fmpz_mod_mpoly_set(
                &mut this.num.raw as *mut _,
                &other.num.raw as *const _ as *mut _,
                &other.num.ctx as *const _ as *mut _,
            );
            // copy denominator
            fmpz_mod_mpoly_set(
                &mut this.den.raw as *mut _,
                &other.den.raw as *const _ as *mut _,
                &other.den.ctx as *const _ as *mut _,
            );
        }
        this
    }
    /// Swap the structure content with another ModMRat provided they have the same context
    /// ```
    /// use flint_mpoly::ModMRat;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat_a = ModMRat::from_str("(x1+x2)/x1",&vars,&modulo).unwrap();
    /// let mut mrat_b = ModMRat::from_str("x2/(x1-x2)",&vars,&modulo).unwrap();
    /// assert_eq!("(+x1+x2)/(+x1)",mrat_a.to_str());
    /// assert_eq!("(+x2)/(+x1+11310*x2)",mrat_b.to_str());
    /// mrat_a.swap(&mut mrat_b);
    /// assert_eq!("(+x1+x2)/(+x1)",mrat_b.to_str());
    /// assert_eq!("(+x2)/(+x1+11310*x2)",mrat_a.to_str());
    /// ```
    pub fn swap(&mut self, other: &mut Self) {
        assert!(
            self.vars == other.vars,
            "Unmatch variables while swapping expressions"
        );
        unsafe {
            fmpz_mod_mpoly_swap(
                &mut self.num.raw as *mut _,
                &mut other.num.raw as *mut _,
                &mut self.num.ctx as *mut _,
            );
            fmpz_mod_mpoly_swap(
                &mut self.den.raw as *mut _,
                &mut other.den.raw as *mut _,
                &mut self.den.ctx as *mut _,
            );
        }
    }
    /// Set function to zero and return the previous value
    /// ```
    /// use flint_mpoly::ModMRat;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat_a = ModMRat::from_str("(x1+x2)/x1",&vars,&modulo).unwrap();
    /// assert_eq!("(+x1+x2)/(+x1)",mrat_a.to_str());
    /// let mut mrat_b = mrat_a.zero_move();
    /// assert_eq!("0",mrat_a.to_str());
    /// assert_eq!("(+x1+x2)/(+x1)",mrat_b.to_str());
    /// ```
    pub fn zero_move(&mut self) -> Self {
        let mut out = Self::new(&self.vars, &self.modulo);
        out.swap(self);
        out
    }
    /// Chec if the function is zero
    pub fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
    /// Check if the function is const
    pub fn is_const(&self) -> bool {
        self.num.is_const() && self.den.is_const()
    }
    /// Compute gcd of num end den
    fn gcd(&self) {
        unsafe {
            fmpz_mod_mpoly_gcd(
                &self._gcd.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
        }
    }
    /// Get GCD bewtween numerator and denominator
    /// ```
    /// use flint_mpoly::ModMRat;
    /// use flint_mpoly::parse_mpoly;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat = ModMRat::new(&vars,&modulo);
    /// mrat.num.set_from_str("(1+x1)*(x2-x1)").unwrap();
    /// mrat.den.set_from_str("(x1+x2)*(x2-x1)").unwrap();
    /// let gcd = mrat.get_gcd();
    /// assert_eq!("+x1+11310*x2", gcd.to_str());
    /// ```
    pub fn get_gcd(&self) -> ModMPoly {
        self.gcd();
        ModMPoly::clone_from(&self._gcd)
    }
    /// Reduce by dividing both numerator and denominators by the gcd.
    /// This expression is not guaranteed to return a denominaor or numerator in monic form.
    /// ```
    /// use flint_mpoly::ModMRat;
    /// use flint_mpoly::parse_mpoly;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat = ModMRat::new(&vars,&modulo);
    /// mrat.num.set_from_str("(1+x1)*(x2-x1)").unwrap();
    /// mrat.den.set_from_str("(x1+x2)*(x2-x1)").unwrap();
    /// mrat.reduce();
    /// assert_eq!("(+11310*x1+11310)/(+11310*x1+11310*x2)", mrat.to_str());
    /// ```
    pub fn reduce(&mut self) {
        self.gcd();
        if self._gcd.is_one() {
            return;
        }
        let mut exact = true;
        unsafe {
            // Simplify the numerator
            exact &= fmpz_mod_mpoly_divides(
                &self._res.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self._gcd.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            ) == 1;
            fmpz_mod_mpoly_swap(
                &self._res.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
            // Simplify the denominator
            exact &= fmpz_mod_mpoly_divides(
                &self._res.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self._gcd.raw as *const _ as *mut _,
                &self.den.ctx as *const _ as *mut _,
            ) == 1;
            fmpz_mod_mpoly_swap(
                &self._res.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self.den.ctx as *const _ as *mut _,
            );
        }
        assert!(exact, "Failed to reduce!");
        self._reduced = true;
    }

    /// Normalize to have a monic polynomial at the denominator by dividing both
    /// numerator and denominator by the denominator leading coefficient
    /// ```
    /// use flint_mpoly::ModMRat;
    /// use flint_mpoly::parse_mpoly;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat = ModMRat::new(&vars,&modulo);
    /// mrat.num.set_from_str("(-3-6*x1)").unwrap();
    /// mrat.den.set_from_str("(-2*x2-3*x1)").unwrap();
    /// mrat.normalize();
    /// assert_eq!("(+2*x1+1)/(+x1+3771*x2)", mrat.to_str());
    /// ```
    pub fn normalize(&mut self) {
        unsafe {
            let length =
                fmpz_mod_mpoly_length(&mut self.den.raw as *mut _, &mut self.den.ctx as *mut _);
            // Get leading coefficient
            if length > 0 {
                fmpz_mod_mpoly_get_term_coeff_fmpz(
                    &mut self._coeff.raw as *mut _,
                    &mut self.den.raw as *mut _,
                    0,
                    &mut self.den.ctx as *const _ as *mut _,
                );
                if !self._coeff.is_zero() && !self._coeff.is_one() {
                    // Get inverse of the coefficient
                    self._coeff.inv();
                    // scale numerator
                    fmpz_mod_mpoly_scalar_mul_fmpz(
                        &mut self._res.raw as *mut _,
                        &mut self.num.raw as *mut _,
                        &mut self._coeff.raw as *mut _,
                        &mut self.num.ctx as *mut _,
                    );
                    fmpz_mod_mpoly_swap(
                        &mut self._res.raw as *mut _,
                        &mut self.num.raw as *mut _,
                        &mut self.num.ctx as *mut _,
                    );
                    // scale denominator
                    fmpz_mod_mpoly_scalar_mul_fmpz(
                        &mut self._res.raw as *mut _,
                        &mut self.den.raw as *mut _,
                        &mut self._coeff.raw as *mut _,
                        &mut self.den.ctx as *mut _,
                    );
                    fmpz_mod_mpoly_swap(
                        &mut self._res.raw as *mut _,
                        &mut self.den.raw as *mut _,
                        &mut self.den.ctx as *mut _,
                    );
                }
            }
        }
    }

    /// Reduce the expression and normalize to have a monic polynomial at the denominator
    /// ```
    /// use flint_mpoly::ModMRat;
    /// use flint_mpoly::parse_mpoly;
    /// let modulo = 11311;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat = ModMRat::new(&vars,&modulo);
    /// mrat.num.set_from_str("(1+x1)*(x2-x1)").unwrap();
    /// mrat.den.set_from_str("(x1+x2)*(x2-x1)").unwrap();
    /// mrat.full_reduce();
    /// assert_eq!("(+x1+1)/(+x1+x2)", mrat.to_str());
    /// ```
    pub fn full_reduce(&mut self) {
        self.reduce();
        self.normalize();
    }
    /// Rise function to integer power
    pub fn pown(&mut self, n: u64) {
        // Bring franction to the simplest form before applying the power
        if !self._reduced {
            self.full_reduce();
        }
        unsafe {
            // Rise numerator
            fmpz_mod_mpoly_pow_ui(
                &mut self._res.raw as *mut _,
                &mut self.num.raw as *mut _,
                n,
                &mut self.num.ctx as *mut _,
            );
            fmpz_mod_mpoly_swap(
                &mut self._res.raw as *mut _,
                &mut self.num.raw as *mut _,
                &mut self.num.ctx as *mut _,
            );
            // Rise denominator
            fmpz_mod_mpoly_pow_ui(
                &mut self._res.raw as *mut _,
                &mut self.den.raw as *mut _,
                n,
                &mut self.den.ctx as *mut _,
            );
            fmpz_mod_mpoly_swap(
                &mut self._res.raw as *mut _,
                &mut self.den.raw as *mut _,
                &mut self.den.ctx as *mut _,
            );
        }
    }
    /// Clear the function and set it to zero
    pub fn clear(&mut self) {
        self.num.set_to_zero();
        self.den.set_to_one();
        // Free the coefficient of the helper functions
        self._gcd.set_to_zero();
        self._res.set_to_zero();
    }
    /// Clear all memory allocated to ModMRat
    pub fn free(&mut self) {
        self.num.free();
        self.den.free();
        // Free the coefficient of the helper functions
        self._gcd.free();
        self._res.free();
    }
    /// Convert the function to a human readable string
    pub fn to_str(&self) -> String {
        let rg = Regex::new(r"([+-])1\*").unwrap();
        rg.replace_all(format!("{}", self).as_str(), "$1")
            .to_string()
    }
}

///Implement addition with operator `+`
/// ```
/// use flint_mpoly::ModMRat;
/// let modulo = 11311;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat_a = ModMRat::from_str("1/x1",&vars,&modulo).unwrap();
/// let mrat_b = ModMRat::from_str("1/x2",&vars,&modulo).unwrap();
/// let mrat_ab = &mrat_a + &mrat_b;
/// assert_eq!("(+x1+x2)/(+x1*x2)", mrat_ab.to_str());
/// ```
impl<'a, 'b> Add<&'b ModMRat> for &'a ModMRat {
    type Output = ModMRat;

    fn add(self, other: &'b ModMRat) -> ModMRat {
        let mut this = ModMRat::new(&self.vars, &self.modulo);
        this.num = &(&self.num * &other.den) + &(&other.num * &self.den);
        this.den = &self.den * &other.den;
        this.full_reduce();
        this
    }
}

///Implement addition with operator `-`
/// ```
/// use flint_mpoly::ModMRat;
/// let modulo = 11311;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat_a = ModMRat::from_str("1/x1",&vars,&modulo).unwrap();
/// let mrat_b = ModMRat::from_str("1/x2",&vars,&modulo).unwrap();
/// let mrat_ab = &mrat_a - &mrat_b;
/// assert_eq!("(+11310*x1+x2)/(+x1*x2)", mrat_ab.to_str());
/// ```
impl<'a, 'b> Sub<&'b ModMRat> for &'a ModMRat {
    type Output = ModMRat;

    fn sub(self, other: &'b ModMRat) -> ModMRat {
        let mut this = ModMRat::new(&self.vars, &self.modulo);
        this.num = &(&self.num * &other.den) - &(&other.num * &self.den);
        this.den = &self.den * &other.den;
        this.full_reduce();
        this
    }
}
///Implement addition with operator `*`
/// ```
/// use flint_mpoly::ModMRat;
/// let modulo = 11311;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat_a = ModMRat::from_str("1/x1",&vars,&modulo).unwrap();
/// let mrat_b = ModMRat::from_str("1/x2",&vars,&modulo).unwrap();
/// let mrat_ab = &mrat_a * &mrat_b;
/// assert_eq!("(+1)/(+x1*x2)", mrat_ab.to_str());
/// ```
impl<'a, 'b> Mul<&'b ModMRat> for &'a ModMRat {
    type Output = ModMRat;

    fn mul(self, other: &'b ModMRat) -> ModMRat {
        let mut this = ModMRat::new(&self.vars, &self.modulo);
        this.num = &self.num * &other.num;
        this.den = &self.den * &other.den;
        this.full_reduce();
        this
    }
}
///Implement addition with operator `/`
/// ```
/// use flint_mpoly::ModMRat;
/// let modulo = 11311;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat_a = ModMRat::from_str("1/x1",&vars,&modulo).unwrap();
/// let mrat_b = ModMRat::from_str("1/x2",&vars,&modulo).unwrap();
/// let mrat_ab = &mrat_a / &mrat_b;
/// assert_eq!("(+x2)/(+x1)", mrat_ab.to_str());
/// ```
impl<'a, 'b> Div<&'b ModMRat> for &'a ModMRat {
    type Output = ModMRat;

    fn div(self, other: &'b ModMRat) -> ModMRat {
        let mut this = ModMRat::new(&self.vars, &self.modulo);
        this.num = &self.num * &other.den;
        this.den = &self.den * &other.num;
        this.full_reduce();
        this
    }
}

///Implement negative sign
/// ```
/// use flint_mpoly::ModMRat;
/// // Define new polynomial
/// let modulo = 11311;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat = ModMRat::from_str("x1/x2",&vars,&modulo).unwrap();
/// let mrat_neg = -&mrat;
/// assert_eq!("(+11310*x1)/(+x2)", mrat_neg.to_str());
/// ```
impl<'a> Neg for &'a ModMRat {
    type Output = ModMRat;
    fn neg(self) -> ModMRat {
        // TODO: Check if the two contexts are the same
        let mut this = ModMRat::new(&self.vars, &self.modulo);
        this.den = self.den.clone();
        unsafe {
            fmpz_mod_mpoly_neg(
                &mut this.num.raw as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

/// Clone for ModMRat
/// ```
/// use flint_mpoly::ModMRat;
/// use std::str::FromStr;
/// let modulo = 11311;
/// let vars = [String::from("x1"),String::from("x2"),String::from("x3")];
/// let mrat_1 = ModMRat::from_str("x1+x3/x2",&vars,&modulo).unwrap();
/// let mrat_2 = mrat_1.clone();
/// assert_eq!(mrat_1.to_str(),mrat_2.to_str());
/// ```
impl Clone for ModMRat {
    fn clone(&self) -> Self {
        ModMRat::clone_from(self)
    }
}

// To use the `{}` for the structure ModMRat
impl fmt::Display for ModMRat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            write!(f, "0")
        } else if self.den.is_one() {
            write!(f, "{}", self.num)
        } else {
            write!(f, "({})/({})", self.num, self.den)
        }
    }
}

// Send and Sync
unsafe impl Send for ModMRat {}
unsafe impl Sync for ModMRat {}
