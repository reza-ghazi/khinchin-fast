"""Khinchin's constant (OEIS A002210) via the accelerated zeta series.

Same mathematics as ../khinchin_fast.c: the Bailey-Borwein-Crandall series

    ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
                   + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),

    h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j,

with the small k pulled out so only M ~ P/(2 log2 N) terms are needed
instead of the ~P/2 of the plain series that mpmath's builtin
`mp.khinchin` uses (mpmath.libmp.gammazeta.khinchin_fixed).  The loop
runs in fixed-point integer arithmetic in the style of that routine, and
zeta(2n) comes from mpmath's cached Bernoulli numbers.

Measured on this machine (Python 3.14, mpmath 1.3.0, single core), time
to first evaluation in a fresh process:

    digits    mp.khinchin    this file    speedup
      1000        4.67 s       0.16 s        30x
      2000       36.2  s       0.98 s        37x

Unlike the C program this is float-guarded, not interval-certified: the
guard bits make the result reliable in practice, but only khinchin-fast
proves its digits.

This port is deliberately serial.  CPython's GIL rules out threads for
CPU-bound bignum work, and multiprocessing does not pay either: the
dominant cost is mpmath's sequential Bernoulli recurrence, which every
worker would have to repeat for its own block, so splitting the n-range
duplicates rather than divides the work.

Usage:  python3 khinchin.py DIGITS [OUTPUT_FILE]        (default 100)

With OUTPUT_FILE the decimal value plus a newline is written there (the
same contract as ../khinchin_fast.c); otherwise it prints to stdout.
"""

import math
import sys

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


# N, M, ONE, K mirror the naming of ../khinchin_fast.c and the series in
# the docstring; the local count is inherent to the numeric kernel.
# pylint: disable=invalid-name,too-many-locals
def khinchin_fixed_accel(prec: int) -> int:
    """K0 as a fixed-point integer scaled by 2**prec."""
    wp = int(prec + prec**0.5 + 32)
    ONE = MPZ_ONE << wp
    N = max(3, int(wp**0.35))
    M = int(wp * math.log(2) / (2 * math.log(N))) + 1

    # Finite logarithmic correction (only ~N logs).
    s = MPZ_ZERO
    for k in range(2, N):
        a = mpf_log(from_rational(k - 1, k, wp), wp)
        b = mpf_log(from_rational(k + 1, k, wp), wp)
        s -= to_fixed(mpf_mul(a, b, wp), wp)

    # Ascending fixed-point loop; absolute errors stay O(1) ulp because
    # the power table is updated by exact floor division.
    powers = [ONE // (k * k) for k in range(2, N)]
    h = ONE
    fac = from_man_exp(4, 0)  # 2 * (2n)!  at n = 1
    pi = mpf_pi(wp)
    pipow = twopi2 = mpf_shift(mpf_mul(pi, pi, wp), 2)  # (2 pi)^(2n)
    for n in range(1, M + 1):
        zeta2n = mpf_abs(mpf_bernoulli(2 * n, wp))
        zeta2n = mpf_mul(zeta2n, pipow, wp)
        zeta2n = mpf_div(zeta2n, fac, wp)
        tail = to_fixed(zeta2n, wp) - ONE - sum(powers)
        s += (tail * h // n) >> wp
        h += ONE // (2 * n + 1) - ONE // (2 * n)
        fac = mpf_mul_int(fac, (2 * n + 2) * (2 * n + 1), wp)
        pipow = mpf_mul(pipow, twopi2, wp)
        for i, k in enumerate(range(2, N)):
            powers[i] //= k * k

    # Truncation tail is < N^(-2M) < 2^(-wp); guard bits absorb it.
    s = (s << wp) // ln2_fixed(wp)
    K = mpf_exp(from_man_exp(s, -wp), wp)
    return to_fixed(K, prec)


def khinchin(digits: int) -> str:
    """Decimal string of K0 with `digits` digits after the point."""
    prec = int((digits + 2) * math.log2(10)) + 32
    value = khinchin_fixed_accel(prec)
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
