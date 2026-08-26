//! Khinchin's constant (OEIS A002210) via the accelerated zeta series.
//!
//! Same mathematics and structure as ../../khinchin_fast.c, including
//! its two-region split:
//!
//! - Bernoulli region (n < n_direct): zeta(2n) is reconstructed from
//!   exact Bernoulli numbers streamed by FLINT's reverse Bernoulli
//!   iterator (`bernoulli_rev`, reached through the C shim in shim.c),
//!   one iterator per block, exactly like the C program.
//! - Direct region (n >= n_direct): the accelerated term is evaluated
//!   as the literal tail sum over k in [N, K], K ~ 2^(wp/2n) <= ~8N —
//!   no Bernoulli numbers at all.
//!
//! Arithmetic is rug (GMP/MPFR) at flat working precision; blocks fan
//! out over rayon, heaviest first.  Guard bits absorb the rounding
//! drift; unlike the C program, nothing here is interval-certified.
//!
//! Run:  cargo run --release -- DIGITS [OUTPUT_FILE]     (default 100)
//!
//! With OUTPUT_FILE the decimal value plus a newline is written there
//! (the same contract as ../../khinchin_fast.c); otherwise stdout.

use std::os::raw::{c_ulong, c_void};

use gmp_mpfr_sys::gmp::mpz_t;
use rayon::prelude::*;
use rug::float::Constant;
use rug::ops::Pow;
use rug::{Float, Integer};

extern "C" {
    fn k_rev_alloc(n: c_ulong) -> *mut c_void;
    fn k_rev_next(num: *mut mpz_t, den: *mut mpz_t, p: *mut c_void);
    fn k_rev_free(p: *mut c_void);
}

/// h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
fn alternating_harmonic(n: u64, wp: u32) -> Float {
    let mut h = Float::new(wp);
    for j in 1..=(2 * n - 1) {
        let t = Float::with_val(wp, 1u32) / j;
        if j % 2 == 1 {
            h += t;
        } else {
            h -= t;
        }
    }
    h
}

/// Accelerated terms for n in [first, last] via exact Bernoulli numbers
/// from FLINT's reverse iterator, descending like scp_bernoulli_range.
fn bernoulli_block(first: u64, last: u64, big_n: u64, wp: u32) -> Float {
    let iter = unsafe { k_rev_alloc(2 * last as c_ulong) };
    let mut num = Integer::new();
    let mut den = Integer::new();

    // factor = (2 pi)^(2n) / (2 * (2n)!) at n = last, updated downward.
    let pi = Float::with_val(wp, Constant::Pi);
    let tau2 = pi.square() * 4u32; // (2 pi)^2
    let inv_tau2 = Float::with_val(wp, 1u32) / &tau2;
    let mut factor = tau2.pow(last as u32)
        / Float::with_val(wp, Integer::from(Integer::factorial(2 * last as u32)))
        / 2u32;

    let mut powers: Vec<Float> = (2..big_n)
        .map(|k| Float::with_val(wp, k * k).pow(last as u32).recip())
        .collect();
    let mut h = alternating_harmonic(last, wp);
    let mut s = Float::new(wp);

    for n in (first..=last).rev() {
        unsafe { k_rev_next(num.as_raw_mut(), den.as_raw_mut(), iter) };
        let mut zeta = Float::with_val(wp, &num) / Float::with_val(wp, &den);
        zeta *= &factor;
        let mut tail = zeta.abs() - 1u32;
        for p in &powers {
            tail -= p;
        }
        s += tail * &h / n;
        if n > first {
            factor *= 2 * n * (2 * n - 1);
            factor *= &inv_tau2;
            for (i, k) in (2..big_n).enumerate() {
                powers[i] *= k * k;
            }
            h += Float::with_val(wp, 1u32) / ((2 * n - 2) * (2 * n - 1));
        }
    }
    unsafe { k_rev_free(iter) };
    s
}

/// Accelerated terms for n in [first, last] as literal tail sums over
/// k in [N, K], descending like scp_direct_range.  No Bernoulli numbers.
fn direct_block(first: u64, last: u64, big_n: u64, wp: u32) -> Float {
    let mut big_k = (2f64.powf((wp as f64 + 32.0) / (2.0 * first as f64)) as u64 + 1)
        .min(16 * big_n + 64);
    big_k = big_k.max(big_n);
    let mut powers: Vec<Float> = (big_n..=big_k)
        .map(|k| Float::with_val(wp, k * k).pow(last as u32).recip())
        .collect();
    let mut h = alternating_harmonic(last, wp);
    let mut s = Float::new(wp);

    for n in (first..=last).rev() {
        let mut tail = Float::new(wp);
        for p in &powers {
            tail += p;
        }
        s += tail * &h / n;
        if n > first {
            for (i, k) in (big_n..=big_k).enumerate() {
                powers[i] *= k * k;
            }
            h += Float::with_val(wp, 1u32) / ((2 * n - 2) * (2 * n - 1));
        }
    }
    s
}

/// K0 with enough guard precision for `digits` correct decimals.
fn khinchin(digits: usize) -> Float {
    let prec = ((digits as f64 + 2.0) * std::f64::consts::LOG2_10).ceil() as u32 + 32;
    let wp = prec + (prec as f64).sqrt().ceil() as u32 + 64;
    let big_n = std::cmp::max(3, (wp as f64).powf(0.35) as u64);
    let log2_n = (big_n as f64).log2();
    let m = (wp as f64 * std::f64::consts::LN_2 / (2.0 * (big_n as f64).ln())).ceil() as u64 + 1;
    let n_direct = ((wp as f64 / (2.0 * (log2_n + 3.0))).ceil() as u64).min(m + 1);

    // Finite logarithmic correction.
    let mut s = Float::new(wp);
    for k in 2..big_n {
        let a = (Float::with_val(wp, k - 1) / k).ln();
        let b = (Float::with_val(wp, k + 1) / k).ln();
        s -= a * b;
    }

    // Equal-count blocks per region, heaviest first: the Bernoulli
    // region's per-n cost grows with n, so its high-n blocks lead.
    let threads = rayon::current_num_threads() as u64;
    let mut blocks: Vec<(bool, u64, u64)> = Vec::new();
    let bern_terms = n_direct - 1;
    let bern_blocks = (4 * threads).min(bern_terms.max(1) / 8).max(1);
    for b in (0..bern_blocks).rev() {
        let first = 1 + bern_terms * b / bern_blocks;
        let last = bern_terms * (b + 1) / bern_blocks;
        if first <= last {
            blocks.push((true, first, last));
        }
    }
    if m >= n_direct {
        let direct_terms = m - n_direct + 1;
        let direct_blocks = (4 * threads).min(direct_terms / 8).max(1);
        for b in 0..direct_blocks {
            let first = n_direct + direct_terms * b / direct_blocks;
            let last = n_direct + direct_terms * (b + 1) / direct_blocks - 1;
            if first <= last {
                blocks.push((false, first, last));
            }
        }
    }

    let partials: Vec<Float> = blocks
        .into_par_iter()
        .map(|(bernoulli, first, last)| {
            if bernoulli {
                bernoulli_block(first, last, big_n, wp)
            } else {
                direct_block(first, last, big_n, wp)
            }
        })
        .collect();
    for p in partials {
        s += p;
    }

    let ln2 = Float::with_val(wp, Constant::Log2);
    (s / ln2).exp()
}

fn main() {
    let mut args = std::env::args().skip(1);
    let digits: usize = args
        .next()
        .map(|a| a.parse().expect("DIGITS must be a positive integer"))
        .unwrap_or(100);
    // rug's format precision counts significant digits; K0 has one
    // integer digit, so digits + 1 gives `digits` after the point.
    let text = format!("{:.*}", digits + 1, khinchin(digits));
    match args.next() {
        Some(path) => std::fs::write(&path, text + "\n").expect("cannot write OUTPUT_FILE"),
        None => println!("{}", text),
    }
}
