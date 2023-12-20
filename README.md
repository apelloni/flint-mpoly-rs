# flint-mpoly-rs
multivariate polynomial for rust with FLINT backend


## Coefficients
 - [x] `QCoeff`: Coefficients that take values in $\mathbb{Q}$ (from: `fmpq.h`)
 - [ ] `ZCoeff`: Coefficients that take values in $\mathbb{Z}$ (from: `fmpz.h`)
 - [ ] `ModCoeff`: Coefficients that take values in $\mathbb{Z}/p\mathbb{Z}$ (from: `fmpz_mod.h`)

## Polynomials
#### Multivariate
 - [x] `QMPoly`: Multivariate Polynomials in $\mathbb{Q}[x_1,x_2,..]$ (from: `fmpq_mpoly.h`)
 - [ ] `ZMPoly`: Multivariate Polynomials in $\mathbb{Z}[x_1,x_2,..]$ (from: `fmpz_mpoly.h`)
 - [ ] `ModMPoly`: Multivariate Polynomials in $\mathbb{Z}/p\mathbb{Z}[x_1,x_2,..]$ (from: `fmpz_mod_mpoly.h`)
#### Univariate
 - [ ] `QPoly`: Univariate Polynomials in $\mathbb{Q}[x]$ (from: `fmpq_poly.h`)
 - [ ] `ZPoly`: Univariate Polynomials in $\mathbb{Z}[x]$ (from: `fmpz_poly.h`)
 - [ ] `ModPoly`: Univariate Polynomials in $\mathbb{Z}/p\mathbb{Z}[x]$ (from: `fmpz_mod_poly.h`)

## Rational Polynomials
#### Multivariate
 - [x] `QMRat`: Multivariate Rational in $\mathbb{Q}(x_1,x_2,..)$
 - [ ] `ZMRat`: Multivariate Rational in $\mathbb{Z}(x_1,x_2,..)$
 - [ ] `ModMRat`: Multivariate Rational in $\mathbb{Z}/p\mathbb{Z}(x_1,x_2,..)$
#### Univariate
 - [ ] `QRat`: Univariate Rational in $\mathbb{Q}(x)$
 - [ ] `ZRat`: Univariate Rational in $\mathbb{Z}(x)$
 - [ ] `ModRat`: Univariate Rational in $\mathbb{Z}/p\mathbb{Z}(x)$

 ---
 # Examples

### Parse an expression
```rust
 use flint_mpoly::QMRat;

 // Define variables
 let vars : Vec<String> = ["x1","x2"].iter().map(|x| x.to_string()).collect();

 // Parse expression from string
 let f = QMRat::from_str("(1+x1+x1^2)/(1-x1^3)*x2",&vars).unwrap();

 // The parsed result is already reduce to the canonical form
 assert_eq!("(-x2)/(+x1-1)",f.to_str());
```

Alternatively one can also add each coefficient individually, for example when
importing it from other structures.
This is defined for polynomials where one can simply add coefficients.

```rust
 use flint_mpoly::{QMPoly, QCoeff};
 use std::str::FromStr;

 // Define variables
 let vars : Vec<String> = ["x1","x2"].iter().map(|x| x.to_string()).collect();

 // Create QMPoly object
 let mut f = QMPoly::new(&vars);  // currently 0

 // Feed coefficients into MPoly usind different methods
 f.add_coeff_pq(&[2,0],1,2);                            // add 1/2*x1^2
 f.add_coeff_str(&[1,1],"5");                           // add 5*x1*x2
 f.add_coeff(&[1,1],&QCoeff::from_str("3/2").unwrap()); // add 3/2*x1*x2

 // The parsed result is already reduce to the canonical form
 assert_eq!("+1/2*x1^2+13/2*x1*x2",f.to_str());
```
