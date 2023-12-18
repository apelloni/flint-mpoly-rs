# flint-mpoly-rs
multivariate polynomial for rust with FLINT backend


## Coefficients
 - [x] `QCoeff`: Coefficients that take values in $\mathbb{Q}$ (from: `fmpq.h`)
 - [ ] `QCoeff`: Coefficients that take values in $\mathbb{Z}$ (from: `fmpz.h`)
 - [ ] `QCoeff`: Coefficients that take values in $\mathbb{Z}/p\mathbb{Z}$ (from: `fmpz_mod.h`)

## Polynomials
 - [x] `QMPoly`: Multivariate Polynomials in $\mathbb{Q}[x_1,x_2,..]$ (from: `fmpq_mpoly.h`)
 - [ ] `QMPoly`: Multivariate Polynomials in $\mathbb{Z}[x_1,x_2,..]$ (from: `fmpz_mpoly.h`)
 - [ ] `QMPoly`: Multivariate Polynomials in $\mathbb{Z}/p\mathbb{Z}[x_1,x_2,..]$ (from: `fmpz_mod_mpoly.h`)
 - [ ] `QMPoly`: Univariate Polynomials in $\mathbb{Q}[x]$ (from: `fmpq_poly.h`)
 - [ ] `QMPoly`: Univariate Polynomials in $\mathbb{Z}[x]$ (from: `fmpz_poly.h`)
 - [ ] `QMPoly`: Univariate Polynomials in $\mathbb{Z}/p\mathbb{Z}[x]$ (from: `fmpz_mod_poly.h`)

## Rational Polynomials
 - [x] `QMRat`: Multivariate Rational in $\mathbb{Q}(x_1,x_2,..)$
 - [ ] `QMRat`: Multivariate Rational in $\mathbb{Z}(x_1,x_2,..)$
 - [ ] `QMRat`: Multivariate Rational in $\mathbb{Z}/p\mathbb{Z}(x_1,x_2,..)$
 - [ ] `QMRat`: Univariate Rational in $\mathbb{Q}(x)$
 - [ ] `QMRat`: Univariate Rational in $\mathbb{Z}(x)$
 - [ ] `QMRat`: Univariate Rational in $\mathbb{Z}/p\mathbb{Z}(x)$
