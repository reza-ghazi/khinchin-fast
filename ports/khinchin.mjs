// Khinchin's constant (OEIS A002210) via the accelerated zeta series.
//
// Same mathematics as ../khinchin_fast.c:
//
//   ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
//                  + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
//
//   h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
//
// Runs entirely on native BigInt in fixed-point (scale 2^wp) — the
// "native bignum" data point for JavaScript, with no dependencies.
// Because BigInt has no transcendentals, the port carries its own
// fixed-point kit: pi by Machin's formula, logarithms via the
// atanh(1/q) series (ln((k±1)/k) = ±2 atanh(1/(2k±1)), ln 2 =
// 2 atanh(1/3)), and exp by argument-halving plus Taylor. The even
// zeta values come from the positive-term recurrence
// (n + 1/2) zeta(2n) = sum_j zeta(2j) zeta(2n-2j) — but only for
// n < n_direct, the C program's two-region split: above it each
// accelerated term IS the tail sum_{k=N..K} k^(-2n) with K <= ~8N,
// evaluated directly with no zeta values at all, which shortens the
// O(M^2) recurrence to O(n_direct^2). Guard bits absorb the drift;
// nothing here is interval-certified.
//
// Usage:  node khinchin.mjs DIGITS [OUTPUT_FILE]        (default 100)
// With OUTPUT_FILE the decimal value plus a newline is written there
// (the same contract as ../khinchin_fast.c); otherwise stdout.

import { writeFileSync } from "node:fs";

// atanh(1/q) = sum_{j>=0} 1 / ((2j+1) q^(2j+1)), fixed point 2^wp.
function atanhRecip(q, ONE) {
  const q2 = q * q;
  let power = ONE / q;
  let s = power;
  for (let j = 1n; power > 0n; j += 1n) {
    power /= q2;
    s += power / (2n * j + 1n);
  }
  return s;
}

// atan(1/q), same series with alternating signs (for Machin's pi).
function atanRecip(q, ONE) {
  const q2 = q * q;
  let power = ONE / q;
  let s = power;
  for (let j = 1n; power > 0n; j += 1n) {
    power /= q2;
    const term = power / (2n * j + 1n);
    s += j % 2n === 1n ? -term : term;
  }
  return s;
}

// exp(y) for 0 <= y < 1: halve the argument `shift` times, Taylor, square back.
function expFixed(y, wp, ONE) {
  const shift = 16n;
  const arg = y >> shift;
  let acc = ONE;
  let term = ONE;
  for (let j = 1n; term > 0n; j += 1n) {
    term = ((term * arg) >> wp) / j;
    acc += term;
  }
  for (let i = 0n; i < shift; i += 1n) {
    acc = (acc * acc) >> wp;
  }
  return acc;
}

// K0 as a fixed-point integer scaled by 2^wp, plus wp.
function khinchinFixed(digits) {
  const prec = BigInt(Math.ceil((digits + 2) * Math.log2(10)) + 32);
  const wp = prec + BigInt(Math.ceil(Math.sqrt(Number(prec)))) + 64n;
  const ONE = 1n << wp;
  const N = Math.max(3, Math.floor(Number(wp) ** 0.35));
  const log2N = Math.log2(N);
  const M = Math.ceil((Number(wp) * Math.LN2) / (2 * Math.log(N))) + 1;
  const nDirect = Math.min(Math.ceil(Number(wp) / (2 * (log2N + 3))), M + 1);

  // Finite logarithmic correction:
  // -ln((k-1)/k) ln((k+1)/k) = 4 atanh(1/(2k-1)) atanh(1/(2k+1)).
  let s = 0n;
  for (let k = 2; k < N; k++) {
    const a = atanhRecip(BigInt(2 * k - 1), ONE);
    const b = atanhRecip(BigInt(2 * k + 1), ONE);
    s += (4n * (a * b)) >> wp;
  }

  // pi = 16 atan(1/5) - 4 atan(1/239); zeta(2) = pi^2/6 seeds the
  // positive-term recurrence (n + 1/2) z_n = sum z_j z_{n-j}.
  const pi = 16n * atanRecip(5n, ONE) - 4n * atanRecip(239n, ONE);
  const z = new Array(nDirect);
  z[1] = ((pi * pi) >> wp) / 6n;
  for (let n = 2; n < nDirect; n++) {
    const half = (n - 1) >> 1;
    let raw = 0n;
    for (let j = 1; j <= half; j++) {
      raw += z[j] * z[n - j];
    }
    let acc = 2n * (raw >> wp);
    if (n % 2 === 0) {
      acc += (z[n / 2] * z[n / 2]) >> wp;
    }
    z[n] = (2n * acc) / BigInt(2 * n + 1);
  }

  // Main accelerated loop, ascending, incremental power table.
  const powers = [];
  for (let k = 2; k < N; k++) {
    powers.push(ONE / BigInt(k * k));
  }
  let h = ONE;
  for (let n = 1; n < nDirect; n++) {
    let tail = z[n] - ONE;
    for (const p of powers) {
      tail -= p;
    }
    s += (tail * h) / BigInt(n) >> wp;
    h += ONE / BigInt(2 * n + 1) - ONE / BigInt(2 * n);
    for (let i = 0, k = 2; k < N; i++, k++) {
      powers[i] /= BigInt(k * k);
    }
  }

  // Direct region: h continues ascending; the power table over
  // k in [N, K] is re-seeded per block of 32 terms.
  let n = nDirect;
  while (n <= M) {
    const last = Math.min(n + 31, M);
    let K = Math.min(
      Math.floor(2 ** ((Number(wp) + 32) / (2 * n))) + 1, 16 * N + 64);
    K = Math.max(K, N);
    const pw = [];
    for (let k = N; k <= K; k++) {
      pw.push(ONE / BigInt(k) ** BigInt(2 * n));
    }
    for (let m = n; m <= last; m++) {
      let tail = 0n;
      for (const p of pw) {
        tail += p;
      }
      s += (tail * h) / BigInt(m) >> wp;
      h += ONE / BigInt(2 * m + 1) - ONE / BigInt(2 * m);
      for (let i = 0, k = N; k <= K; i++, k++) {
        pw[i] /= BigInt(k * k);
      }
    }
    n = last + 1;
  }

  const ln2 = 2n * atanhRecip(3n, ONE);
  return [expFixed((s << wp) / ln2, wp, ONE), wp];
}

// Decimal string with `digits` digits after the point.
function khinchin(digits) {
  const [value, wp] = khinchinFixed(digits);
  const scaled = (value * 10n ** BigInt(digits) + (1n << (wp - 1n))) >> wp;
  const text = scaled.toString();
  return text[0] + "." + text.slice(1);
}

const digits = process.argv[2] ? parseInt(process.argv[2], 10) : 100;
const result = khinchin(digits);
if (process.argv[3]) {
  writeFileSync(process.argv[3], result + "\n");
} else {
  console.log(result);
}
