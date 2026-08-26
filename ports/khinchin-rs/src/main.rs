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
//! Parallel with rayon: the recurrence convolutions (the dominant cost,
//! halved via their j <-> n-j symmetry) fan out across the thread pool,
//! and the main loop splits into blocks seeded like the C program's,
//! with h(first) summed directly and the power table started at
//! k^(-2 first).
//!
//! Run:  cargo run --release -- DIGITS [OUTPUT_FILE]     (default 100)
//!
//! With OUTPUT_FILE the decimal value plus a newline is written there
//! (the same contract as ../../khinchin_fast.c); otherwise stdout.

use rayon::prelude::*;
use rug::float::Constant;
use rug::ops::Pow;
use rug::Float;

/// K0 with enough guard precision for `digits` correct decimals.
fn khinchin(digits: usize) -> Float {
    let prec = ((digits as f64 + 2.0) * std::f64::consts::LOG2_10).ceil() as u32 + 32;
    let wp = prec + (prec as f64).sqrt().ceil() as u32 + 64;
    let big_n = std::cmp::max(3, (wp as f64).powf(0.35) as u64);
    let m = (wp as f64 * std::f64::consts::LN_2 / (2.0 * (big_n as f64).ln())).ceil() as usize + 1;

    // Finite logarithmic correction.
    let mut s = Float::new(wp);
    for k in 2..big_n {
        let a = (Float::with_val(wp, k - 1) / k).ln();
        let b = (Float::with_val(wp, k + 1) / k).ln();
        s -= a * b;
    }

    // zeta(2), zeta(4), ..., zeta(2M) by the positive recurrence, using
    // the symmetry of the convolution: only j <= (n-1)/2 is summed,
    // doubled, plus the middle square when n is even.
    let mut z: Vec<Float> = Vec::with_capacity(m);
    z.push(Float::with_val(wp, Constant::Pi).square() / 6u32);
    for n in 2..=m {
        let half = (n - 1) / 2;
        let mut acc = if half >= 64 {
            (1..=half)
                .into_par_iter()
                .fold(
                    || Float::new(wp),
                    |mut a, j| {
                        a += &z[j - 1] * &z[n - j - 1];
                        a
                    },
                )
                .reduce(
                    || Float::new(wp),
                    |mut a, b| {
                        a += b;
                        a
                    },
                )
        } else {
            let mut a = Float::new(wp);
            for j in 1..=half {
                a += &z[j - 1] * &z[n - j - 1];
            }
            a
        };
        acc *= 2u32;
        if n % 2 == 0 {
            acc += z[n / 2 - 1].clone().square();
        }
        z.push(acc * 2u32 / (2 * n as u64 + 1));
    }

    // Main loop in equal-count blocks, one rayon task per block.
    let w = std::cmp::max(1, std::cmp::min(rayon::current_num_threads(), m / 8));
    let partials: Vec<Float> = (0..w)
        .into_par_iter()
        .map(|b| {
            let first = 1 + m * b / w;
            let last = m * (b + 1) / w;
            let mut h = Float::new(wp);
            for j in 1..=(2 * first as u64 - 1) {
                let t = Float::with_val(wp, 1u32) / j;
                if j % 2 == 1 {
                    h += t;
                } else {
                    h -= t;
                }
            }
            let mut powers: Vec<Float> = (2..big_n)
                .map(|k| Float::with_val(wp, k * k).pow(first as u32).recip())
                .collect();
            let mut s = Float::new(wp);
            for n in first..=last {
                let mut tail = z[n - 1].clone() - 1u32;
                for p in &powers {
                    tail -= p;
                }
                s += tail * &h / n as u64;
                h -= Float::with_val(wp, 1u32) / (2 * n as u64 * (2 * n as u64 + 1));
                for (i, k) in (2..big_n).enumerate() {
                    powers[i] /= k * k;
                }
            }
            s
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
