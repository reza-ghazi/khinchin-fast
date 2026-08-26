"""Khinchin's constant: two-region threaded variant of khinchin.py.

Same accelerated series as ../khinchin_fast.c and khinchin.py, but with
the C program's two-region split, which is what makes threads useful in
Python at all:

- Bernoulli region (n < n_direct): zeta(2n) is near 1 and needs mpmath's
  full-precision Bernoulli numbers, whose cached recurrence is strictly
  sequential.  Runs on the main thread.
- Direct region (n >= n_direct): the accelerated term IS the tail
  zeta(2n) - 1 - sum_{k<N} k^(-2n) = sum_{k>=N} k^(-2n), and to the
  needed accuracy only k <= K ~ 2^(wp/2n) <= ~8N contribute.  Each block
  of n is an independent, mpmath-free, pure-integer computation —
  submitted to a thread pool that runs concurrently with the Bernoulli
  region.

The split pays twice: the sequential Bernoulli recurrence covers only
n < n_direct (~2/3 of the range, and its cost is superlinear in the top
index), and the direct region parallelises.  On the GIL build the
workers serialise but the shorter recurrence still helps; on a
free-threaded build (see .venv-ft/) they genuinely overlap.  Worker code
touches no mpmath state, so thread safety is trivial on both.

Measured at 10,000 digits (24 cores; serial khinchin.py in
parentheses): pure-Python backend 39.9 s (64.5 s), gmpy2 backend 5.0 s
(9.0 s).  The attribution is humbling: a single-worker run on gmpy2
takes 4.85 s, so essentially the whole 1.8x is the restructure — the
direct region is cheap once written as short tail sums, and what
remains is the sequential Bernoulli region.  Amdahl's law, exactly as
khinchin.py's docstring predicts.

Usage:  python3 khinchin_mt.py DIGITS [OUTPUT_FILE]      (default 100)
"""

import math
import os
import sys
from concurrent.futures import ThreadPoolExecutor

from mpmath.libmp import (
    MPZ_ONE,
    MPZ_ZERO,
    from_man_exp,
    from_rational,
    mpf_abs,
    mpf_bernoulli,
    mpf_div,
    mpf_exp,
    mpf_log,
    mpf_mul,
    mpf_mul_int,
    mpf_pi,
    mpf_shift,
    to_fixed,
)
from mpmath.libmp.gammazeta import ln2_fixed

# mpmath's pure-Python backend converts huge ints to decimal strings
# internally; lift CPython's 4300-digit int/str safety limit or runs
# beyond ~1300 digits raise ValueError.
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


# N, ONE, K in the helpers mirror ../khinchin_fast.c and the series
# notation, as in khinchin.py.
# pylint: disable=invalid-name
def _alternating_harmonic(first, ONE):
    """h(first) = sum_{j=1}^{2*first-1} (-1)^(j+1)/j in fixed point."""
    h = MPZ_ZERO
    for j in range(1, 2 * first):
        h += ONE // j if j % 2 == 1 else -(ONE // j)
    return h


def _direct_block(first, last, N, wp, ONE):
    """Fixed-point sum of accelerated terms for n in [first, last].

    Pure integer arithmetic: no mpmath calls, so blocks are trivially
    thread-safe and, on a free-threaded interpreter, truly parallel.
    """
    K = min(int(2.0 ** ((wp + 32) / (2.0 * first))) + 1, 16 * N + 64)
    K = max(K, N)
    ks = range(N, K + 1)
    # k^(-2*first); floor underflow to 0 for large k is harmless (the
    # true value is below one ulp) and the ascending k*k divisions only
    # shrink the entries further.
    powers = [ONE // k ** (2 * first) for k in ks]
    h = _alternating_harmonic(first, ONE)
    s = MPZ_ZERO
    for n in range(first, last + 1):
        s += sum(powers) * h // n >> wp
        h += ONE // (2 * n + 1) - ONE // (2 * n)
        for i, k in enumerate(ks):
            powers[i] //= k * k
    return s


# pylint: disable=too-many-locals
def khinchin_fixed_two_region(prec, max_workers=None):
    """K0 as a fixed-point integer scaled by 2**prec."""
    wp = int(prec + prec**0.5 + 32)
    ONE = MPZ_ONE << wp
    N = max(3, int(wp**0.35))
    log2_N = math.log2(N)
    M = int(wp * math.log(2) / (2 * math.log(N))) + 1
    n_direct = min(int(math.ceil(wp / (2.0 * (log2_N + 3.0)))), M + 1)

    # Fan the direct region out first so it overlaps the sequential
    # Bernoulli region below.  Blocks of equal n-count, lowest (heaviest,
    # since K falls with n) first.
    workers = max_workers or os.cpu_count() or 1
    executor = ThreadPoolExecutor(max_workers=workers)
    futures = []
    # At least ~8 terms per block, or the per-block h(first) and power
    # table seeding outweighs the work.
    blocks = min(4 * workers, max(1, (M - n_direct + 1) // 8))
    for b in range(blocks):
        first = n_direct + (M + 1 - n_direct) * b // blocks
        last = n_direct + (M + 1 - n_direct) * (b + 1) // blocks - 1
        if first <= last:
            futures.append(executor.submit(_direct_block, first, last, N, wp, ONE))

    # Finite logarithmic correction (only ~N logs).
    s = MPZ_ZERO
    for k in range(2, N):
        a = mpf_log(from_rational(k - 1, k, wp), wp)
        b = mpf_log(from_rational(k + 1, k, wp), wp)
        s -= to_fixed(mpf_mul(a, b, wp), wp)

    # Bernoulli region on the main thread, exactly as in khinchin.py but
    # stopping at n_direct - 1.
    powers = [ONE // (k * k) for k in range(2, N)]
    h = ONE
    fac = from_man_exp(4, 0)  # 2 * (2n)!  at n = 1
    pi = mpf_pi(wp)
    pipow = twopi2 = mpf_shift(mpf_mul(pi, pi, wp), 2)  # (2 pi)^(2n)
    for n in range(1, n_direct):
        zeta2n = mpf_abs(mpf_bernoulli(2 * n, wp))
        zeta2n = mpf_mul(zeta2n, pipow, wp)
        zeta2n = mpf_div(zeta2n, fac, wp)
        tail = to_fixed(zeta2n, wp) - ONE - sum(powers)
        s += tail * h // n >> wp
        h += ONE // (2 * n + 1) - ONE // (2 * n)
        fac = mpf_mul_int(fac, (2 * n + 2) * (2 * n + 1), wp)
        pipow = mpf_mul(pipow, twopi2, wp)
        for i, k in enumerate(range(2, N)):
            powers[i] //= k * k

    for future in futures:
        s += future.result()
    executor.shutdown()

    s = (s << wp) // ln2_fixed(wp)
    K = mpf_exp(from_man_exp(s, -wp), wp)
    return to_fixed(K, prec)


def khinchin(digits, max_workers=None):
    """Decimal string of K0 with `digits` digits after the point."""
    prec = int((digits + 2) * math.log2(10)) + 32
    value = khinchin_fixed_two_region(prec, max_workers)
    scaled = (value * 10**digits + (MPZ_ONE << (prec - 1))) >> prec
    text = str(scaled)
    return text[0] + "." + text[1:]


if __name__ == "__main__":
    result = khinchin(int(sys.argv[1]) if len(sys.argv) > 1 else 100)
    if len(sys.argv) > 2:
        with open(sys.argv[2], "w", encoding="ascii") as output:
            output.write(result + "\n")
    else:
        print(result)
