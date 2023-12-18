use crate::qcoeff::QCoeff;
use crate::qmpoly::QMPoly;
use crate::qmrat::QMRat;
//use num::{BigInt, BigRational};
use pest::iterators::Pairs;
use pest::Parser;
use std::str::FromStr;

#[derive(Parser)]
#[grammar = "grammar.pest"]
pub struct INIParser;

pub fn parse_mrat(expr: &str, var_names: &[String]) -> Result<QMRat, &'static str> {
    let parsed = INIParser::parse(Rule::calculation, expr)
        .expect("Unsuccessful Parse") // unwrap the parse result
        .next()
        .unwrap(); // get and unwrap the `file` rule; never fails
    let parsed_mrat = parsed_eval_rat(parsed.into_inner(), var_names).unwrap();
    Ok(parsed_mrat)
}

pub fn parse_mpoly(expr: &str, var_names: &[String]) -> Result<QMPoly, &'static str> {
    let parsed = INIParser::parse(Rule::calculation, expr)
        .expect("Unsuccessful Parse") // unwrap the parse result
        .next()
        .unwrap(); // get and unwrap the `file` rule; never fails
    let parsed_mrat = parsed_eval_rat(parsed.into_inner(), var_names).unwrap();
    if parsed_mrat.den.is_one() {
        Ok(parsed_mrat.num)
    } else {
        Err("Expected polynomial found rational function")
    }
}

pub fn parsed_eval_rat(parser_rules: Pairs<Rule>, vars: &[String]) -> Result<QMRat, &'static str> {
    let mut block = QMRat::new(vars);
    let mut out = QMRat::new(vars);
    let mut res = QMRat::new(vars);
    let mut pows = vec![0; vars.len()];
    let mut exponent: i32;
    let mut last_rule = Rule::add;
    let mut start_of_input = true;
    let mut found_var;
    let mut dpower = false;
    for parsed_rule in parser_rules.into_iter() {
        //println!("{:?} \"{}\"", parsed_rule.as_rule(), parsed_rule.as_str());
        match parsed_rule.as_rule() {
            Rule::num | Rule::var | Rule::expr => {
                // Process adjacent expressions as if there is an
                // implicit multiplicaiton in between of them
                if last_rule == Rule::num || last_rule == Rule::var || last_rule == Rule::expr {
                    last_rule = Rule::multiply;
                }

                match (last_rule, parsed_rule.as_rule()) {
                    // Process the information in the expression in terms of polynomial
                    (Rule::power, Rule::num | Rule::expr) => {
                        // Process the power of an expression considering
                        // that one power its already stored in the block
                        // so we need to generate the remaining (power-1)
                        // contirbutions
                        exponent = parsed_rule
                            .clone()
                            .into_inner()
                            .as_str()
                            .parse::<i32>()
                            .expect("Exponent Parsing Error");
                        //assert!(exponent > 0, "Exponents must be positive!");
                        if dpower {
                            if exponent <= 0 {
                                res.pown((-exponent + 1) as u64);
                                last_rule = Rule::multiply;
                            }
                            if exponent > 0 {
                                res.pown((exponent - 1) as u64);
                                last_rule = Rule::divide;
                            }
                        } else {
                            if exponent <= 0 {
                                res.pown((-exponent + 1) as u64);
                                last_rule = Rule::divide;
                            }
                            if exponent > 0 {
                                res.pown((exponent - 1) as u64);
                                last_rule = Rule::multiply;
                            }
                        }
                    }
                    (_, Rule::num) => {
                        let coeff = QCoeff::from_str(parsed_rule.clone().into_inner().as_str())
                            .expect("bad string");
                        for p in pows.iter_mut() {
                            *p = 0;
                        }
                        res.clear();
                        res.num.add_coeff(&pows, &coeff);
                        dpower = last_rule == Rule::divide;
                    }
                    (_, Rule::var) => {
                        found_var = false;
                        let coeff = QCoeff::one();
                        for (vn, p) in pows.iter_mut().enumerate() {
                            if parsed_rule.clone().into_inner().as_str() == vars[vn] {
                                *p = 1;
                                found_var = true;
                            } else {
                                *p = 0;
                            }
                        }
                        if !found_var {
                            panic!(
                                "Unknown variable {:?} found. [ vars={:?} ]",
                                parsed_rule.into_inner().as_str(),
                                vars
                            );
                        }
                        res.clear();
                        res.num.add_coeff(&pows, &coeff);
                        dpower = last_rule == Rule::divide;
                    }
                    (_, Rule::expr) => {
                        res = parsed_eval_rat(parsed_rule.into_inner(), vars)?;
                        dpower = last_rule == Rule::divide;
                    }
                    _ => panic!("Impossible!"),
                };
                // append operation to the block based on `last_rule`
                match last_rule {
                    Rule::add => {
                        block = QMRat::clone_from(&res);
                    }
                    Rule::subtract => {
                        block.clear();
                        block = &block - &res;
                    }
                    Rule::multiply => {
                        block = &block * &res;
                    }
                    Rule::divide => {
                        block = &block / &res;
                    }
                    _ => {
                        panic!("Unknown operation sequence! .. {:?}  expr ..", last_rule)
                    }
                };
                last_rule = Rule::expr;
            }
            // When we have addition or subtraction we evaluate the block and add it to
            // the output polynomial
            Rule::add | Rule::subtract => {
                match last_rule {
                    Rule::num | Rule::var | Rule::expr => {
                        if start_of_input {
                            start_of_input = false;
                            out = QMRat::clone_from(&block);
                        } else {
                            out = &out + &block;
                        }
                        block.clear();
                        last_rule = parsed_rule.as_rule();
                    }
                    Rule::add => {
                        last_rule = if parsed_rule.as_rule() == Rule::add {
                            Rule::add
                        } else {
                            Rule::subtract
                        };
                    }
                    Rule::subtract => {
                        last_rule = if parsed_rule.as_rule() == Rule::subtract {
                            Rule::add
                        } else {
                            Rule::subtract
                        };
                    }
                    _ => {
                        panic!(
                            "Unknown operation sequence! .. {:?}  {:?} ..",
                            last_rule,
                            parsed_rule.as_rule()
                        )
                    }
                };
            }
            // The information about the binary operations *|/|^ that are performed within
            // the block are used only a posteriori as information in `last_rule`
            // The double multiplication is treated as a power operator
            Rule::multiply | Rule::divide | Rule::power => {
                match last_rule {
                    Rule::num | Rule::var | Rule::expr => {
                        last_rule = parsed_rule.as_rule();
                    }
                    Rule::multiply => last_rule = Rule::power,
                    _ => {
                        panic!(
                            "Unknown operation sequence! .. {:?}  {:?} ..",
                            last_rule,
                            parsed_rule.as_rule()
                        )
                    }
                };
            }
            x => {
                println!("{:?}", x);
                println!(" {:?}", parsed_rule.as_span());
                last_rule = x;
            }
        };
        //println!(
        //    " -> {} [{}]",
        //    out.to_str(var_names),
        //    block.to_str(var_names)
        //);
    }
    if !block.is_zero() {
        out = &out + &block;
    }
    //println!("{}", out);
    //println!("{}", out.to_str();
    Ok(out)
}

#[cfg(test)]
mod test {
    use super::*;
    //
    #[test]
    fn color_eyre() {
        color_eyre::install().unwrap();
    }
    #[test]
    #[should_panic]
    fn parsing_bad_string_1() {
        let var_names = vec![String::from("x")];
        let expr_str = "+10//1";
        parse_mpoly(expr_str, &var_names).expect("parsing polynomial");
    }
    #[test]
    #[should_panic]
    fn parsing_bad_string_2() {
        let var_names = vec![String::from("x")];
        let expr_str = "+10/1*y";
        parse_mpoly(expr_str, &var_names).expect("parsing polynomial");
    }
    #[test]
    #[should_panic]
    fn parsing_bad_string_3() {
        let var_names = vec![String::from("x")];
        let expr_str = "+10/1/x";
        parse_mpoly(expr_str, &var_names).expect("parsing polynomial");
    }
    #[test]
    fn parsing_test_polynomial_1() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "10+a^3*(7+a-a^2)*(w+1)^0";
        println!("input: {}", expr_str);
        let expr_parsed = parse_mpoly(expr_str, &var_names).expect("parsing polynomial");
        println!("output: {}", expr_parsed);
        assert_eq!("-a^5+a^4+7*a^3+10", expr_parsed.to_str());
    }
    #[test]
    fn parsing_test_polynomial_2() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "1/2*(1-w/4+w*(1/2-3/4))*4+w^2*(w^23*a^2*w^5)";
        let expr_parsed = parse_mpoly(expr_str, &var_names).expect("parsing polynomial");
        assert_eq!("+w^30*a^2-w+2", expr_parsed.to_str());
    }
    #[test]
    fn parsing_test_polynomial_3() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "1/(1-45^3-4*w*w+(2*w)^2)";
        let expr_parsed = parse_mpoly(expr_str, &var_names).expect("parsing polynomial");
        assert_eq!("-1/91124", expr_parsed.to_str());
    }
    #[test]
    fn parsing_test_polynomial_4() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "0+(w)*(-1)+(-w)*(-1)+(-w)+(w)+(-2)*(-2)+(-w)+(2 + w)";
        let expr_parsed = parse_mpoly(expr_str, &var_names).expect("parsing polynomial");
        assert_eq!("+6", expr_parsed.to_str());
    }
    #[test]
    fn parsing_test_polynomial_5() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "w^0";
        let expr_parsed = parse_mpoly(expr_str, &var_names).expect("parsing polynomial");
        assert_eq!("+1", expr_parsed.to_str());
    }

    #[test]
    fn parsing_test_mrat_1() {
        let var_names = vec![String::from("w"), String::from("a")];

        let expr_str = "(1-w)^2/(1-w)";
        let expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        assert_eq!("(-w+1)/(+1)", expr_parsed.to_str());
    }
    #[test]
    fn parsing_test_mrat_2() {
        let var_names = vec![String::from("w"), String::from("a")];
        //let expr_str = "w^0 ";
        let expr_str = "10+a^3*(7+a-a^2)*(w+1)^0";
        let expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        assert_eq!("(-a^5+a^4+7*a^3+10)/(+1)", expr_parsed.to_str());
    }
    #[test]
    fn parsing_test_mrat_3() {
        let var_names = vec![String::from("w"), String::from("a")];
        // TEST
        let expr_str = "1/2*(1-w/4+w*(1/2-3/4))*4+w^2*(w^23*a^2*w^5)";
        println!("input: {}", expr_str);
        let expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        assert_eq!("(+w^30*a^2-w+2)/(+1)", expr_parsed.to_str());
    }
    #[test]
    fn parsing_test_mrat_4() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "1/(1-45^3-4*w*w+(2*w)^2)";
        println!("input: {}", expr_str);
        let expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        println!("output: {}", expr_parsed.to_str());
        assert_eq!("(-1/91124)/(+1)", expr_parsed.to_str());
    }
    #[test]
    fn parsing_test_mrat_5() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "0+(w)*(-1)+(-w)*(-1)+(-w)+(w)+(-2)*(-2)+(-w)+(2 + w)";
        println!("input: {}", expr_str);
        let expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        println!("output: {}", expr_parsed.to_str());
        assert_eq!("(+6)/(+1)", expr_parsed.to_str());
    }
    #[test]
    fn parsing_test_mrat_6() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "w^0";
        println!("input: {}", expr_str);
        let mut expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        expr_parsed.reduce();
        println!("output: {}", expr_parsed.to_str());
    }

    #[test]
    fn minus_sign() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "-(1-w)";
        let expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        assert_eq!("+w-1", expr_parsed.num.to_str());
        assert_eq!("+1", expr_parsed.den.to_str());
    }

    #[test]
    fn whitespaces() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "( -     2  *           w  ^  3  + 4 )";
        let expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        println!("{}", expr_parsed.num);
        assert_eq!("-2*w^3+4", expr_parsed.num.to_str());
        assert_eq!("+1", expr_parsed.den.to_str());
    }

    #[test]
    fn implicit_multiplication_1() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "10a^3*2(7+a-a^2)*(w+1)^0";
        let expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        assert_eq!("-20*a^5+20*a^4+140*a^3", expr_parsed.num.to_str());
        assert_eq!("+1", expr_parsed.den.to_str());
    }
    #[test]
    fn implicit_multiplication_2() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "1/2(1-w/4+w(1/2-3/4))4+w^2*(w^23a^2w^5)";
        let expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        assert_eq!("+w^30*a^2-w+2", expr_parsed.num.to_str());
        assert_eq!("+1", expr_parsed.den.to_str());
    }

    #[test]
    fn power_as_double_star_1() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "10+a**3*(7+a-a**2)*(w+1)**0";
        let expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        assert_eq!("(-a^5+a^4+7*a^3+10)/(+1)", expr_parsed.to_str());
    }
    #[test]
    fn power_as_double_star_2() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "1/2*(1-w/4+w*(1/2-3/4))*4+w**2*(w**23*a**2*w**5)";
        let expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        assert_eq!("+w^30*a^2-w+2", expr_parsed.num.to_str());
        assert_eq!("+1", expr_parsed.den.to_str());
    }
    #[test]
    fn power_as_double_star_3() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "1/(1-45**3-4*w*w+(2*w)**2)";
        let expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        assert_eq!("-1/91124", expr_parsed.num.to_str());
        assert_eq!("+1", expr_parsed.den.to_str());
    }
    #[test]
    fn power_as_double_star_4() {
        let var_names = vec![String::from("w"), String::from("a")];
        let expr_str = "w**0";
        let expr_parsed = parse_mrat(expr_str, &var_names).expect("parsing rational");
        assert_eq!("+1", expr_parsed.num.to_str());
        assert_eq!("+1", expr_parsed.den.to_str());
    }

    #[test]
    fn parse_power_polynomial_num() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "x1/2^3";
        let res_str = "1/8*x1";
        let expr = parse_mrat(expr_str, &var_names).expect("parsing rational");
        let res = parse_mrat(res_str, &var_names).expect("parsing rational");
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn parse_dpower_polynomial_var() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "x1^4/x1^3";
        let res_str = "x1";
        let expr = parse_mrat(expr_str, &var_names).expect("parsing rational");
        let res = parse_mrat(res_str, &var_names).expect("parsing rational");
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn parse_dpower_polynomial_expr() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "x1/(2+1)^3";
        let res_str = "1/27*x1";
        let expr = parse_mrat(expr_str, &var_names).expect("parsing rational");
        let res = parse_mrat(res_str, &var_names).expect("parsing rational");
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn parse_dpower_polynomial_sum() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "x1/(2+1)^3+x1*9^2/3/3+x2^4";
        let res_str = "1/27*x1+9*x1+x2^4";
        let expr = parse_mrat(expr_str, &var_names).expect("parsing rational");
        let res = parse_mrat(res_str, &var_names).expect("parsing rational");
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn parse_dpower_expr_num() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "(1-x1)/(1+3)^3";
        let res_str = "1/4/4/4*(1-x1)";
        let expr = parse_mrat(expr_str, &var_names).expect("parsing rational");
        let res = parse_mrat(res_str, &var_names).expect("parsing rational");
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn parse_dpower_expr_var() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "(1-x1)/(x2)^2";
        let res_str = "(1-x1)/x2/x2";
        let expr = parse_mrat(expr_str, &var_names).expect("parsing rational");
        let res = parse_mrat(res_str, &var_names).expect("parsing rational");
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn parse_dpower_expr_expr() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        let expr_str = "(1-x1)^2/(x2+1)^2*(x2-1)^2";
        let res_str = "(1-x1)*(1-x1)/(x2+1)/(x2+1)*(x2-1)*(x2-1)";
        let expr = parse_mrat(expr_str, &var_names).expect("parsing rational");
        let res = parse_mrat(res_str, &var_names).expect("parsing rational");
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn geometric_division() {
        let var_names = vec![String::from("x1")];
        // Geometric Series
        let mpoly_string = "(1-x1^11)^2/(1-x1)^2";
        let res_str = "(1+x1+x1^2+x1^3+x1^4+x1^5+x1^6+x1^7+x1^8+x1^9+x1^10)^2";
        let expr = parse_mrat(mpoly_string, &var_names).expect("parsing rational");
        let res = parse_mrat(res_str, &var_names).expect("parsing rational");
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn univariate_rational() {
        let var_names = ["x1".to_string()];

        let mut expr = parse_mrat(
            "27/(29 *(37/(28 * x1)+(28 * x1)/2479))+37/(28 * x1)+(28 * x1)/2479",
            &var_names,
        )
        .expect("parsing rational");
        assert_eq!(
            "(+28/2479*x1^4+69079/812*x1^2+3393751/21952)/(+x1^3+91723/784*x1)",
            expr.to_str()
        );
    }

    #[test]
    fn multivariate_unreduceble() {
        let var_names = vec![String::from("x1"), String::from("x2")];
        // Multivariate
        let expr_str = "(+x1+x2)/(+x1-x2)";
        let mut expr = parse_mrat(expr_str, &var_names).expect("parsing rational");
        assert_eq!(expr_str, expr.to_str());
    }

    #[test]
    fn multivariate_big() {
        let var_names = [
            "x1".to_string(),
            "x2".to_string(),
            "x3".to_string(),
            "x4".to_string(),
            "a".to_string(),
        ];

        let mut input_string = String::new();
        input_string += "(((-2*x1+x3-x4)*(1-2*x1+x3-x4))/(1+(x1+a*x2-x3)/(-x1+x2/";
        input_string += "a+x3))+(x1+x2/a+x4)/(1+(x1+a*x2-x3)/(1-2*x1+x3-x4)))/((-";
        input_string += "2*x1+x3-x4)^2+(x1+a*x2-x3)*(-2*x1+x3-x4)^2+(-2*x1+x3-x4)";
        input_string += "^3+(x1+a*x2-x3)*(x1+x2/a+x4)+2*(-2*x1+x3-x4)*(x1+x2/a+x4";
        input_string += ")+(x1+a*x2-x3)*(-2*x1+x3-x4)*(x1+x2/a+x4)+(-2*x1+x3-x4)^";
        input_string += "2*(x1+x2/a+x4)+(x1+x2/a+x4)^2)";
        let expr = parse_mrat(&input_string, &var_names).expect("parsing rational");
        //expr.reduce();
        assert_eq!(
            "(-2*x1*a+x3*a-x4*a+a)/(+x2^2*a^3-x1*x2*a^2-x2*x4*a^2+x2^2*a+x2*a^2-x1*x2-x2*x4+x2)",
            expr.to_str()
        );
    }

    #[test]
    fn multivariate_rational_1() {
        let var_names = [
            "x1".to_string(),
            "x2".to_string(),
            "x3".to_string(),
            "x4".to_string(),
            "a".to_string(),
        ];
        let input_string = "11905/(1904*a)+(27368*a)/9849+x1+x2+x3^(-1)";
        let mut res_string = String::new();
        res_string += "(+1904*x3*a+11905*x3^2+7444096/1407*x3^2*a^2";
        res_string += "+1904*x2*x3^2*a+1904*x1*x3^2*a)/(+1904*x3^2*a)";
        let expr = parse_mrat(input_string, &var_names).expect("parsing rational");
        let res = parse_mrat(&res_string, &var_names).expect("parsing rational");
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn multivariate_rational_2() {
        let var_names = [
            "x1".to_string(),
            "x2".to_string(),
            "x3".to_string(),
            "x4".to_string(),
            "a".to_string(),
        ];

        let input_string = "+1/(x1^(-1)+x2^(-1)+x3^(-1)+x4^(-1))";
        let res_string = "(x1*x2*x3*x4)/(x1*x2*x3+x1*x2*x4+x1*x3*x4+x2*x3*x4)";
        let expr = parse_mrat(input_string, &var_names).expect("parsing rational");
        let res = parse_mrat(res_string, &var_names).expect("parsing rational");
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn multivariate_rational_3() {
        let var_names = [
            "x1".to_string(),
            "x2".to_string(),
            "x3".to_string(),
            "x4".to_string(),
            "a".to_string(),
        ];

        let input_string = "+((77*a)/16+x1^(-1)+x2+x3+x4^(-1))^(-1)";
        let mut res_string = String::new();
        res_string += "(+x1^5*x4^6+x1^6*x4^5+77/16*x1^6*x4^6*a+x1^6*x3*x4^6";
        res_string += "+x1^6*x2*x4^6)/(+x1^4*x4^6+2*x1^5*x4^5+77/8*x1^5*x4^6*";
        res_string += "a+2*x1^5*x3*x4^6+2*x1^5*x2*x4^6+x1^6*x4^4+77/8*x1^";
        res_string += "6*x4^5*a+5929/256*x1^6*x4^6*a^2+2*x1^6*x3*x4^5+77/8*";
        res_string += "x1^6*x3*x4^6*a+x1^6*x3^2*x4^6+2*x1^6*x2*x4^5+77/8*";
        res_string += "x1^6*x2*x4^6*a+2*x1^6*x2*x3*x4^6+x1^6*x2^2*x4^6)";
        let expr = parse_mrat(input_string, &var_names).expect("parsing rational");
        let res = parse_mrat(&res_string, &var_names).expect("parsing rational");
        println!("{}", expr.to_str());
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn multivariate_rational_4() {
        let var_names = [
            "x1".to_string(),
            "x2".to_string(),
            "x3".to_string(),
            "x4".to_string(),
            "a".to_string(),
        ];

        let input_string = "+1/(x1^(-1)+x2^(-1))";
        let res_string = "(x1*x2)/(x1+x2)";
        let expr = parse_mrat(input_string, &var_names).expect("parsing rational");
        let res = parse_mrat(res_string, &var_names).expect("parsing rational");
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn multivariate_rational_5() {
        let var_names = [
            "x1".to_string(),
            "x2".to_string(),
            "x3".to_string(),
            "x4".to_string(),
            "a".to_string(),
        ];

        let input_string = "+1/(x1^(-1)+x2^(-1)+x3^(-1))";
        let res_string = "(x1*x2*x3)/(x2*x3+x1*x3+x1*x2)";
        let expr = parse_mrat(input_string, &var_names).expect("parsing rational");
        let res = parse_mrat(res_string, &var_names).expect("parsing rational");
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn multivariate_rational_6() {
        let var_names = [
            "x1".to_string(),
            "x2".to_string(),
            "x3".to_string(),
            "x4".to_string(),
            "a".to_string(),
        ];

        let input_string = "+1/(x1^(-1))";
        let res_string = "x1";
        let expr = parse_mrat(input_string, &var_names).expect("parsing rational");
        let res = parse_mrat(res_string, &var_names).expect("parsing rational");
        assert_eq!(res.to_str(), expr.to_str());
    }

    #[test]
    fn parse_pown() {
        let var_names = vec![String::from("x1"), String::from("x2")];

        let expr_str = "(x1 + x2 )";
        let mut expr_1 = parse_mrat(expr_str, &var_names).expect("parsing rational");
        expr_1.pown(200);
        let expr_str_pow = "(x1 + x2 )^200";
        let expr_2 = parse_mrat(expr_str_pow, &var_names).expect("parsing rational");
        assert_eq!(expr_1.to_str(), expr_2.to_str());
    }
}
