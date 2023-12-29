use crate::parse_mrat;
use crate::{QCoeff, QMPoly, ZMPoly, ZMRat};
use flint_sys::fmpq_mpoly::*;
use regex::Regex;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Container for rational polynomials over rational numbers
#[derive(Debug)]
pub struct QMRat {
    pub num: QMPoly,       // numerator
    pub den: QMPoly,       // Polynomial context
    pub vars: Vec<String>, // Name of variables
    _reduced: bool,
    _coeff: QCoeff, // coefficient container
    _gcd: QMPoly,   // gcd container
    _res: QMPoly,   // partial operations container
}

impl QMRat {
    /// Create a rational function equal to 0
    /// ```
    /// use flint_mpoly::QMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mrat = QMRat::new(&vars);
    /// assert_eq!("0",mrat.to_str());
    /// ```
    pub fn new(vars: &[String]) -> Self {
        let num = QMPoly::new(vars);
        let mut den = QMPoly::new(vars);
        den.set_to_one();
        QMRat {
            num,
            den,
            vars: vars.to_vec(),
            _reduced: false,
            _coeff: QCoeff::default(),
            _gcd: QMPoly::new(vars),
            _res: QMPoly::new(vars),
        }
    }

    /// Create a rational function from a polynomial
    /// ```
    /// use flint_mpoly::QMPoly;
    /// use flint_mpoly::QMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly= QMPoly::from_str("x1+3/4 x2",&vars).unwrap();
    /// let mrat = QMRat::from_mpoly(&mpoly);
    /// assert_eq!("+x1+3/4*x2",mrat.to_str());
    /// ```
    pub fn from_mpoly(mpoly: &QMPoly) -> Self {
        let mut den = QMPoly::new(&mpoly.vars);
        den.set_to_one();
        QMRat {
            num: QMPoly::clone_from(mpoly),
            den,
            vars: mpoly.vars.clone(),
            _reduced: false,
            _coeff: QCoeff::default(),
            _gcd: QMPoly::new(&mpoly.vars),
            _res: QMPoly::new(&mpoly.vars),
        }
    }

    /// Create a rational function from two polynomials defining the numerator and denominator,
    /// respectively
    /// ```
    /// use flint_mpoly::QMPoly;
    /// use flint_mpoly::QMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let num= QMPoly::from_str("x1+3/4 x2",&vars).unwrap();
    /// let den= QMPoly::from_str("x1-x2",&vars).unwrap();
    /// let mrat = QMRat::from_mpolys(&num,&den);
    /// assert_eq!("(+x1+3/4*x2)/(+x1-x2)",mrat.to_str());
    /// ```
    pub fn from_mpolys(num: &QMPoly, den: &QMPoly) -> Self {
        QMRat {
            num: QMPoly::clone_from(num),
            den: QMPoly::clone_from(den),
            vars: num.vars.clone(),
            _reduced: false,
            _coeff: QCoeff::default(),
            _gcd: QMPoly::new(&num.vars),
            _res: QMPoly::new(&num.vars),
        }
    }

    /// Initialize rational from string assuming the given variables
    /// ```
    /// use flint_mpoly::QMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mrat = QMRat::from_str("(x1+x2)/x1", &vars).unwrap();
    /// assert_eq!("(+x1+x2)/(+x1)",mrat.to_str());
    /// ```
    pub fn from_str(expr_str: &str, vars: &[String]) -> Result<Self, &'static str> {
        Ok(parse_mrat(expr_str, vars).expect("parse rational"))
    }

    /// Set rational from string assuming the variables already stored in QMRat
    /// ```
    /// use flint_mpoly::QMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat = QMRat::new(&vars);
    /// mrat.set_from_str("(x1+x2)/x1").unwrap();
    /// assert_eq!("(+x1+x2)/(+x1)",mrat.to_str());
    /// ```
    pub fn set_from_str(&mut self, expr_str: &str) -> Result<(), &'static str> {
        *self = QMRat::from_str(expr_str, &self.vars)?;
        Ok(())
    }

    /// Clone will copy the pointers and we will have two objects pointing at the same memory
    /// address. With this function we can create a clone in a safe way
    /// ```
    /// use flint_mpoly::QMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mrat_1 = QMRat::from_str("x1/(x2+1)",&vars).unwrap();
    /// let mrat_2 = QMRat::clone_from(&mrat_1);
    /// assert_eq!(mrat_1.to_str(),mrat_2.to_str());
    /// ```
    pub fn clone_from(other: &Self) -> Self {
        let mut this = QMRat::new(&other.vars);
        unsafe {
            // copy numerator
            fmpq_mpoly_set(
                &mut this.num.raw as *mut _,
                &other.num.raw as *const _ as *mut _,
                &other.num.ctx as *const _ as *mut _,
            );
            // copy denominator
            fmpq_mpoly_set(
                &mut this.den.raw as *mut _,
                &other.den.raw as *const _ as *mut _,
                &other.den.ctx as *const _ as *mut _,
            );
        }
        this
    }
    /// Swap the structure content with another QMRat provided they have the same context
    /// ```
    /// use flint_mpoly::QMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat_a = QMRat::from_str("(x1+x2)/x1",&vars).unwrap();
    /// let mut mrat_b = QMRat::from_str("x2/(x1-x2)",&vars).unwrap();
    /// assert_eq!("(+x1+x2)/(+x1)",mrat_a.to_str());
    /// assert_eq!("(+x2)/(+x1-x2)",mrat_b.to_str());
    /// mrat_a.swap(&mut mrat_b);
    /// assert_eq!("(+x1+x2)/(+x1)",mrat_b.to_str());
    /// assert_eq!("(+x2)/(+x1-x2)",mrat_a.to_str());
    /// ```
    pub fn swap(&mut self, other: &mut Self) {
        assert!(
            self.vars == other.vars,
            "Unmatch variables while swapping expressions"
        );
        unsafe {
            fmpq_mpoly_swap(
                &mut self.num.raw as *mut _,
                &mut other.num.raw as *mut _,
                &mut self.num.ctx as *mut _,
            );
            fmpq_mpoly_swap(
                &mut self.den.raw as *mut _,
                &mut other.den.raw as *mut _,
                &mut self.den.ctx as *mut _,
            );
        }
    }
    /// Set function to zero and return the previous value
    /// ```
    /// use flint_mpoly::QMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat_a = QMRat::from_str("(x1+x2)/x1",&vars).unwrap();
    /// assert_eq!("(+x1+x2)/(+x1)",mrat_a.to_str());
    /// let mut mrat_b = mrat_a.zero_move();
    /// assert_eq!("0",mrat_a.to_str());
    /// assert_eq!("(+x1+x2)/(+x1)",mrat_b.to_str());
    /// ```
    pub fn zero_move(&mut self) -> Self {
        let mut out = Self::new(&self.vars);
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
            fmpq_mpoly_gcd(
                &self._gcd.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
        }
    }
    /// Get GCD bewtween numerator and denominator
    /// ```
    /// use flint_mpoly::QMRat;
    /// use flint_mpoly::parse_mpoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat = QMRat::new(&vars);
    /// mrat.num.set_from_str("(1+x1)*(x2-x1)").unwrap();
    /// mrat.den.set_from_str("(x1+x2)*(x2-x1)").unwrap();
    /// let gcd = mrat.get_gcd();
    /// assert_eq!("+x1-x2", gcd.to_str());
    /// ```
    pub fn get_gcd(&self) -> QMPoly {
        self.gcd();
        QMPoly::clone_from(&self._gcd)
    }
    /// Reduce by dividing both numerator and denominators by the gcd.
    /// This expression is not guaranteed to return a denominaor or numerator in monic form.
    /// ```
    /// use flint_mpoly::QMRat;
    /// use flint_mpoly::parse_mpoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat = QMRat::new(&vars);
    /// mrat.num.set_from_str("(1+x1)*(x2-x1)").unwrap();
    /// mrat.den.set_from_str("(x1+x2)*(x2-x1)").unwrap();
    /// mrat.reduce();
    /// assert_eq!("(-x1-1)/(-x1-x2)", mrat.to_str());
    /// ```
    pub fn reduce(&mut self) {
        self.gcd();
        if self._gcd.is_one() {
            return;
        }
        let mut exact = true;
        unsafe {
            // Simplify the numerator
            exact &= fmpq_mpoly_divides(
                &self._res.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self._gcd.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            ) == 1;
            fmpq_mpoly_swap(
                &self._res.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
            // Simplify the denominator
            exact &= fmpq_mpoly_divides(
                &self._res.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self._gcd.raw as *const _ as *mut _,
                &self.den.ctx as *const _ as *mut _,
            ) == 1;
            fmpq_mpoly_swap(
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
    /// use flint_mpoly::QMRat;
    /// use flint_mpoly::parse_mpoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat = QMRat::new(&vars);
    /// mrat.num.set_from_str("(-3-6*x1)").unwrap();
    /// mrat.den.set_from_str("(-2*x2-3*x1)").unwrap();
    /// mrat.normalize();
    /// assert_eq!("(+2*x1+1)/(+x1+2/3*x2)", mrat.to_str());
    /// ```
    pub fn normalize(&mut self) {
        unsafe {
            let length =
                fmpq_mpoly_length(&mut self.den.raw as *mut _, &mut self.den.ctx as *mut _);
            // Get leading coefficient
            if length > 0 {
                fmpq_mpoly_get_term_coeff_fmpq(
                    &mut self._coeff.raw as *mut _,
                    &mut self.den.raw as *mut _,
                    0,
                    &mut self.den.ctx as *const _ as *mut _,
                );
                if !self._coeff.is_zero() && !self._coeff.is_one() {
                    // scale numerator
                    fmpq_mpoly_scalar_div_fmpq(
                        &mut self._res.raw as *mut _,
                        &mut self.num.raw as *mut _,
                        &mut self._coeff.raw as *mut _,
                        &mut self.num.ctx as *mut _,
                    );
                    fmpq_mpoly_swap(
                        &mut self._res.raw as *mut _,
                        &mut self.num.raw as *mut _,
                        &mut self.num.ctx as *mut _,
                    );
                    // scale denominator
                    fmpq_mpoly_scalar_div_fmpq(
                        &mut self._res.raw as *mut _,
                        &mut self.den.raw as *mut _,
                        &mut self._coeff.raw as *mut _,
                        &mut self.den.ctx as *mut _,
                    );
                    fmpq_mpoly_swap(
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
    /// use flint_mpoly::QMRat;
    /// use flint_mpoly::parse_mpoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat = QMRat::new(&vars);
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
            fmpq_mpoly_pow_ui(
                &mut self._res.raw as *mut _,
                &mut self.num.raw as *mut _,
                n,
                &mut self.num.ctx as *mut _,
            );
            fmpq_mpoly_swap(
                &mut self._res.raw as *mut _,
                &mut self.num.raw as *mut _,
                &mut self.num.ctx as *mut _,
            );
            // Rise denominator
            fmpq_mpoly_pow_ui(
                &mut self._res.raw as *mut _,
                &mut self.den.raw as *mut _,
                n,
                &mut self.den.ctx as *mut _,
            );
            fmpq_mpoly_swap(
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
    /// Clear all memory allocated to QMRat
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
    /// Convert to ZMRat
    pub fn to_zmrat(&mut self) -> Result<ZMRat, &'static str> {
        unsafe {
            // Get numerical denominator of the QMPoly at the numerator
            fmpq_mpoly_get_denominator(
                &self._coeff.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
            // Rescale numerator
            fmpq_mpoly_scalar_mul_fmpq(
                &self._res.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self._coeff.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
            fmpq_mpoly_swap(
                &self._res.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
            // Rescale denominator
            fmpq_mpoly_scalar_mul_fmpq(
                &self._res.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self._coeff.raw as *const _ as *mut _,
                &self.den.ctx as *const _ as *mut _,
            );
            fmpq_mpoly_swap(
                &self._res.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self.den.ctx as *const _ as *mut _,
            );

            // Get numerical denominator of the QMPoly at the denominator
            fmpq_mpoly_get_denominator(
                &self._coeff.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self.den.ctx as *const _ as *mut _,
            );
            // Rescale numerator
            fmpq_mpoly_scalar_mul_fmpq(
                &self._res.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self._coeff.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
            fmpq_mpoly_swap(
                &self._res.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
            // Rescale denominator
            fmpq_mpoly_scalar_mul_fmpq(
                &self._res.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self._coeff.raw as *const _ as *mut _,
                &self.den.ctx as *const _ as *mut _,
            );
            fmpq_mpoly_swap(
                &self._res.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self.den.ctx as *const _ as *mut _,
            );
        }

        // Fill ZMRat
        let mut res = ZMRat::new(&self.vars);
        res.num = self.num.to_zmpoly()?;
        res.den = self.den.to_zmpoly()?;
        Ok(res)
    }
}

///Implement addition with operator `+`
/// ```
/// use flint_mpoly::QMRat;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat_a = QMRat::from_str("1/x1",&vars).unwrap();
/// let mrat_b = QMRat::from_str("1/x2",&vars).unwrap();
/// let mrat_ab = &mrat_a + &mrat_b;
/// assert_eq!("(+x1+x2)/(+x1*x2)", mrat_ab.to_str());
/// ```
impl<'a, 'b> Add<&'b QMRat> for &'a QMRat {
    type Output = QMRat;

    fn add(self, other: &'b QMRat) -> QMRat {
        let mut this = QMRat::new(&self.vars);
        this.num = &(&self.num * &other.den) + &(&other.num * &self.den);
        this.den = &self.den * &other.den;
        this.full_reduce();
        this
    }
}

///Implement addition with operator `-`
/// ```
/// use flint_mpoly::QMRat;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat_a = QMRat::from_str("1/x1",&vars).unwrap();
/// let mrat_b = QMRat::from_str("1/x2",&vars).unwrap();
/// let mrat_ab = &mrat_a - &mrat_b;
/// assert_eq!("(-x1+x2)/(+x1*x2)", mrat_ab.to_str());
/// ```
impl<'a, 'b> Sub<&'b QMRat> for &'a QMRat {
    type Output = QMRat;

    fn sub(self, other: &'b QMRat) -> QMRat {
        let mut this = QMRat::new(&self.vars);
        this.num = &(&self.num * &other.den) - &(&other.num * &self.den);
        this.den = &self.den * &other.den;
        this.full_reduce();
        this
    }
}
///Implement addition with operator `*`
/// ```
/// use flint_mpoly::QMRat;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat_a = QMRat::from_str("1/x1",&vars).unwrap();
/// let mrat_b = QMRat::from_str("1/x2",&vars).unwrap();
/// let mrat_ab = &mrat_a * &mrat_b;
/// assert_eq!("(+1)/(+x1*x2)", mrat_ab.to_str());
/// ```
impl<'a, 'b> Mul<&'b QMRat> for &'a QMRat {
    type Output = QMRat;

    fn mul(self, other: &'b QMRat) -> QMRat {
        let mut this = QMRat::new(&self.vars);
        this.num = &self.num * &other.num;
        this.den = &self.den * &other.den;
        this.full_reduce();
        this
    }
}
///Implement addition with operator `/`
/// ```
/// use flint_mpoly::QMRat;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat_a = QMRat::from_str("1/x1",&vars).unwrap();
/// let mrat_b = QMRat::from_str("1/x2",&vars).unwrap();
/// let mrat_ab = &mrat_a / &mrat_b;
/// assert_eq!("(+x2)/(+x1)", mrat_ab.to_str());
/// ```
impl<'a, 'b> Div<&'b QMRat> for &'a QMRat {
    type Output = QMRat;

    fn div(self, other: &'b QMRat) -> QMRat {
        let mut this = QMRat::new(&self.vars);
        this.num = &self.num * &other.den;
        this.den = &self.den * &other.num;
        this.full_reduce();
        this
    }
}

///Implement negative sign
/// ```
/// use flint_mpoly::QMRat;
/// // Define new polynomial
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat = QMRat::from_str("x1/x2",&vars).unwrap();
/// let mrat_neg = -&mrat;
/// assert_eq!("(-x1)/(+x2)", mrat_neg.to_str());
/// ```
impl<'a> Neg for &'a QMRat {
    type Output = QMRat;
    fn neg(self) -> QMRat {
        // TODO: Check if the two contexts are the same
        let mut this = QMRat::new(&self.vars);
        this.den = self.den.clone();
        unsafe {
            fmpq_mpoly_neg(
                &mut this.num.raw as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

/// Clone for QMRat
/// ```
/// use flint_mpoly::QMRat;
/// use std::str::FromStr;
/// let vars = [String::from("x1"),String::from("x2"),String::from("x3")];
/// let mrat_1 = QMRat::from_str("x1+x3/x2",&vars).unwrap();
/// let mrat_2 = mrat_1.clone();
/// assert_eq!(mrat_1.to_str(),mrat_2.to_str());
/// ```
impl Clone for QMRat {
    fn clone(&self) -> Self {
        QMRat::clone_from(self)
    }
}

// To use the `{}` for the structure QMRat
impl fmt::Display for QMRat {
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
unsafe impl Send for QMRat {}
unsafe impl Sync for QMRat {}
