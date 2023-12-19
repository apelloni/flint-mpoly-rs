// Import structures
use flint_mpoly::{QCoeff, QMPoly, QMRat};

// Allow import from string
use std::str::FromStr;

// Use parallel
use rayon::prelude::*;

// These tests will not even compile if one of the type doesn't support Send or/and Sync
mod send_sync {
    use super::*;
    use flint_sys::fmpq::fmpq;
    use flint_sys::fmpq_mpoly::{fmpq_mpoly_ctx_struct, fmpq_mpoly_struct};
    fn send_check<T: Send>() {}
    fn sync_check<T: Sync>() {}
    #[test]
    fn send_qcoeff_raw() {
        send_check::<fmpq>();
        sync_check::<fmpq>();
    }
    #[test]
    fn send_qcoeff() {
        send_check::<QCoeff>();
        sync_check::<QCoeff>();
    }
    //#[test]
    //fn send_qmpoly_raw() {
    //    send_check::<fmpq_mpoly_struct>();
    //    sync_check::<fmpq_mpoly_struct>();
    //}
    #[test]
    fn send_qmpoly_ctx() {
        send_check::<fmpq_mpoly_ctx_struct>();
        sync_check::<fmpq_mpoly_ctx_struct>();
    }
    #[test]
    fn send_qmpoly() {
        send_check::<QMPoly>();
        sync_check::<QMPoly>();
    }
    #[test]
    fn send_qmrat() {
        send_check::<QMRat>();
        sync_check::<QMRat>();
    }
}

mod parallel {
    use super::*;
    #[test]
    fn set_qcoeff() {
        let mut v_from = vec![
            QCoeff::from_str("1/2").unwrap(),
            QCoeff::from_str("2/3").unwrap(),
        ];
        let mut v_to = vec![QCoeff::default(); 2];
        v_to.par_iter_mut()
            .enumerate()
            .for_each(|(n, c)| *c = v_from[n].clone());
        for (from, to) in v_from.iter().zip(v_to.iter()) {
            assert_eq!(from.to_str(), to.to_str());
        }
    }
    #[test]
    fn set_qmpoly() {
        let vars: Vec<String> = ["x1", "x2", "x3"].iter().map(|x| x.to_string()).collect();
        let v_from = vec![
            QMPoly::from_str("x1/2", &vars).unwrap(),
            QMPoly::from_str("x2/3", &vars).unwrap(),
        ];
        let mut v_to = vec![QMPoly::new(&vars); 2];
        v_to.par_iter_mut()
            .enumerate()
            .for_each(|(n, c)| *c = v_from[n].clone());
        for (from, to) in v_from.iter().zip(v_to.iter()) {
            assert_eq!(from.to_str(), to.to_str());
        }
    }
    #[test]
    fn set_qmrat() {
        let vars: Vec<String> = ["x1", "x2", "x3"].iter().map(|x| x.to_string()).collect();
        let v_from = vec![
            QMRat::from_str("x1/x2", &vars).unwrap(),
            QMRat::from_str("x2/x3", &vars).unwrap(),
        ];
        let mut v_to = vec![QMRat::new(&vars); 2];
        v_to.par_iter_mut()
            .enumerate()
            .for_each(|(n, c)| *c = v_from[n].clone());
        for (from, to) in v_from.iter().zip(v_to.iter()) {
            assert_eq!(from.to_str(), to.to_str());
        }
    }
}
