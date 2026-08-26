//! Khinchin's constant (OEIS A002210) via the accelerated zeta series.
//!
//! Same mathematics as ../../khinchin_fast.c:
//!
//! ```text
//! ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
//!                + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
//!
//! h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
//! ```
//!
//! Arbitrary precision comes from `rug` (GMP/MPFR).  The even zeta
//! values use the classical positive-term recurrence
//!
//! ```text
//! (n + 1/2) zeta(2n) = sum_{j=1}^{n-1} zeta(2j) zeta(2n-2j),  zeta(2) = pi^2/6,
//! ```
//!
//! which is numerically benign (all terms positive, error growth about
//! linear in n) and keeps the port free of Bernoulli numbers.  Guard
//! bits absorb the drift; unlike the C program, nothing here is
//! interval-certified.
//!
//! Run:  cargo run --release -- [digits]        (default 100)

use rug::float::Constant;
use rug::Float;

/// K0 with enough guard precision for `digits` correct decimals.
fn khinchin(digits: usize) -> Float {
    let prec = ((digits as f64 + 2.0) * std::f64::consts::LOG2_10).ceil() as u32 + 32;
    let wp = prec + (prec as f64).sqrt().ceil() as u32 + 64;
    let big_n = std::cmp::max(3, (wp as f64).powf(0.35) as u64);
    let m = (wp as f64 * std::f64::consts::LN_2 / (2.0 * (big_n as f64).ln())).ceil() as u64 + 1;

    // Finite logarithmic correction.
    let mut s = Float::new(wp);
    for k in 2..big_n {
        let a = (Float::with_val(wp, k - 1) / k).ln();
        let b = (Float::with_val(wp, k + 1) / k).ln();
        s -= a * b;
    }

    // zeta(2), zeta(4), ..., zeta(2M) by the positive recurrence.
    let mut z: Vec<Float> = Vec::with_capacity(m as usize);
    z.push(Float::with_val(wp, Constant::Pi).square() / 6u32);
    for n in 2..=m as usize {
        let mut acc = Float::new(wp);
        for j in 1..n {
            acc += &z[j - 1] * &z[n - j - 1];
        }
        z.push(acc * 2u32 / (2 * n as u64 + 1));
    }

    let mut powers: Vec<Float> = (2..big_n)
        .map(|k| Float::with_val(wp, 1u32) / (k * k))
        .collect();
    let mut h = Float::with_val(wp, 1u32);
    for n in 1..=m {
        let mut tail = z[n as usize - 1].clone() - 1u32;
        for p in &powers {
            tail -= p;
        }
        s += tail * &h / n;
        h -= Float::with_val(wp, 1u32) / (2 * n * (2 * n + 1));
        for (i, k) in (2..big_n).enumerate() {
            powers[i] /= k * k;
        }
    }

    let ln2 = Float::with_val(wp, Constant::Log2);
    (s / ln2).exp()
}

fn main() {
    let digits: usize = std::env::args()
        .nth(1)
        .map(|a| a.parse().expect("DIGITS must be a positive integer"))
        .unwrap_or(100);
    // rug's format precision counts significant digits; K0 has one
    // integer digit, so digits + 1 gives `digits` after the point.
    println!("{:.*}", digits + 1, khinchin(digits));
}
