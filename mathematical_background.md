# Khinchin's Constant: Mathematical Background

This document explains the mathematics behind `khinchin-fast`: what Khinchin's
constant is, where it comes from, the history of its study and computation, the
formulas that make high-precision evaluation possible, how implementations
(including this one) work, and the state of digit-computation records.

## 1. What the constant says

Write a real number `x` as a regular continued fraction

```
x = a0 + 1/(a1 + 1/(a2 + 1/(a3 + ...))),      a0 integer, a1, a2, ... positive integers.
```

The integers `a1, a2, a3, ...` are the *partial quotients* (or coefficients) of
`x`. Khinchin's theorem (1934/35) states that for **almost every** real number
`x` — that is, for all `x` outside a set of Lebesgue measure zero — the
geometric mean of the first `n` partial quotients converges to one universal
constant, independent of `x`:

```
lim_{n->oo} (a1 a2 ... an)^(1/n)  =  K0  =  2.685452001065306445309714835481...
```

This is remarkable: the partial quotients of a "typical" number look wildly
irregular, yet their geometric mean is the same for almost every number. The
100 digits in `khinchin-100.txt` are the reference value used by this project's
tests.

The exceptional (measure-zero) set contains many familiar numbers:

- rationals (finite continued fraction);
- quadratic irrationals, whose expansions are eventually periodic — e.g. the
  golden ratio `phi = [1; 1, 1, 1, ...]` has geometric mean 1;
- Euler's number `e = [2; 1, 2, 1, 1, 4, 1, 1, 6, ...]`, whose regular pattern
  forces the geometric mean of its partial quotients to infinity along a
  different growth law.

Whether specific constants satisfy the theorem is unknown. Numerical evidence
is consistent with typicality for `pi`, `ln 2`, `ln 3`, the Euler–Mascheroni
constant `gamma`, Apéry's constant `zeta(3)`, the Feigenbaum constants, and —
self-referentially — `K0` itself (which would in particular make `K0`
irrational), but no proof exists for any of them. Numbers *constructed* to
satisfy the theorem are known (Wieting's "Khinchin sequence", 2008), so the
property is provably realizable; the open problem is establishing it for any
naturally occurring constant.

## 2. Mathematical basis: the Gauss map and its invariant measure

The engine behind the theorem is the *Gauss map* on the unit interval,

```
T(x) = 1/x - floor(1/x),        T(0) = 0,
```

which shifts continued fractions: if `x = [0; a1, a2, a3, ...]` then
`T(x) = [0; a2, a3, ...]`, and the first partial quotient is
`a1(x) = floor(1/x)`. Gauss discovered (in a notebook entry, and in an 1812
letter to Laplace) that the measure

```
dmu(x) = dx / ((1 + x) ln 2)
```

is invariant under `T`. Under this *Gauss measure*, the probability that a
partial quotient equals `k` is the **Gauss–Kuzmin distribution**:

```
P(a = k) = -log2( 1 - 1/(k+1)^2 )  =  log2( (k+1)^2 / (k(k+2)) ).
```

So `P(1) = log2(4/3) ~ 0.4150`, `P(2) = log2(9/8) ~ 0.1699`, and the tail
decays like `1/k^2` — slow enough that the *arithmetic* mean of partial
quotients diverges almost everywhere, which is why the geometric mean is the
natural statistic.

The Gauss map is ergodic with respect to `mu`. Applying the Birkhoff pointwise
ergodic theorem to the observable `f(x) = ln a1(x)` gives, for almost every
`x`,

```
(1/n) sum_{i=1}^{n} ln a_i  ->  Integral of ln a1 dmu  =  sum_{k>=1} ln(k) * P(a = k),
```

and exponentiating yields Khinchin's theorem with the closed form

```
K0 = prod_{k=1}^{oo} ( 1 + 1/(k(k+2)) )^(log2 k).
```

Khinchin's original 1935 proof did not use ergodic theory (he worked with
quasi-independence estimates for the digits); the clean ergodic-theoretic proof
is due to Ryll-Nardzewski (1951), after Knopp and others had developed the
measure-theoretic toolkit.

## 3. History

- **1812** — Gauss states, in a letter to Laplace, that he can prove the
  limiting distribution of continued-fraction remainders but has lost the
  proof; he asks about the error term. This is the *Gauss problem*.
- **1928** — R. O. Kuzmin solves Gauss's problem, proving convergence to the
  Gauss measure with an error bound `O(q^sqrt(n))`.
- **1929** — Paul Lévy independently proves it with the sharper geometric rate
  `O(q^n)`, `q = 0.7...`; the optimal constant (the Gauss–Kuzmin–Wirsing
  constant, ~0.3036) was found by Wirsing in 1974.
- **1934/1935** — Aleksandr Yakovlevich Khinchin publishes the metric theory of
  continued fractions ("Metrische Kettenbruchprobleme", Compositio
  Mathematica), proving the geometric-mean theorem and introducing what is now
  called Khinchin's constant. His short book *Continued Fractions* (1935; 3rd
  ed. 1961) remains the classic reference.
- **1951** — Ryll-Nardzewski recasts the proof via ergodic theory.
- **1959–1960s** — first serious numerical evaluations: Shanks and Wrench
  ("Khintchine's constant", *Amer. Math. Monthly* 66, 1959), extended by
  Wrench (*Math. Comp.* 14, 1960) — limited to modest precision because the
  defining product converges extremely slowly.
- **1997** — Bailey, Borwein and Crandall ("On the Khintchine constant",
  *Mathematics of Computation* 66) derive fast zeta-series accelerations and
  compute 7,350 digits.
- **1998** — Xavier Gourdon computes 110,000 digits with a refined version of
  the same family of series.
- **2000s–2020s** — incremental extensions by individuals (notably Remco
  Bloemen's `recmo/khinchin` project, whose reverse-Bernoulli idea this program
  adapts), reaching into the hundreds of thousands of digits. No quasi-linear
  algorithm has been found, so records remain tiny compared to constants like
  `pi`.

## 4. Formulas

### 4.1 Defining product (slow)

```
K0 = prod_{k=1}^{oo} ( 1 + 1/(k(k+2)) )^(log2 k)
```

The k-th factor is `1 + O(1/k^2)` with a `log k` exponent, so getting `d`
digits needs about `10^d` factors — hopeless beyond a few digits. Every
practical method starts from this product and accelerates it.

### 4.2 The Bailey–Borwein–Crandall zeta series

Taking logarithms and expanding `ln(1 + 1/(k(k+2)))` into series that can be
resummed against the zeta function gives

```
ln(2) * ln(K0) = sum_{n=1}^{oo}  (zeta(2n) - 1)/n  *  h(n),

h(n) = sum_{j=1}^{2n-1} (-1)^(j+1) / j       (alternating harmonic prefix),
```

Since `zeta(2n) - 1 ~ 2^(-2n)`, the series converges geometrically: about
`P/2` terms give `P` bits. This alone is what most library implementations
use.

### 4.3 Acceleration by pulling out small k (what this program computes)

Convergence is governed by the smallest `k` in the zeta tails. Subtracting the
first `N-2` powers from each zeta value and compensating with an explicit
logarithmic correction gives, for any `N >= 2`,

```
ln(2) * ln(K0) = - sum_{k=2}^{N-1} ln((k-1)/k) * ln((k+1)/k)
                 + sum_{n=1}^{oo} ( zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n) ) / n * h(n).
```

Now the n-th term has magnitude about `N^(-2n)`, so only
`M ~ P / (2 log2 N)` terms are needed for `P` bits, and the omitted tail is
rigorously bounded by `N^(-2M)`. This program picks `N ~ P^0.35` and tracks
every bound in Arb ball arithmetic, so the printed digits are certified.

Two further identities matter to the implementation:

- **Bernoulli representation of even zeta values** (used for small `n`, where
  `zeta(2n)` is near 1 and needs full precision):

  ```
  zeta(2n) = (-1)^(n+1) * (2 pi)^(2n) * B_{2n} / (2 * (2n)!)
  ```

  FLINT's reverse Bernoulli iterator produces the exact `B_{2n}` in descending
  order while sharing one power table — the "SCP" trick.

- **Direct tail identity** (used for large `n`, where this program departs
  from the classical implementations):

  ```
  zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)  =  sum_{k>=N} k^(-2n),
  ```

  and to absolute accuracy `2^(-P)` the right side needs only
  `k <= 2^(P/2n)` — barely more than `N` entries near the top of the range,
  all at low precision, with no Bernoulli numbers and no cancellation.

### 4.4 Hölder means and related constants

Khinchin's theorem generalizes to Hölder means of order `p < 1`:

```
K_p = ( sum_{k>=1} k^p * P(a = k) )^(1/p)     for almost every x,
```

with `K_0` (the limit `p -> 0`) the geometric mean above. The weights are the
Gauss–Kuzmin probabilities, so explicitly

```
K_p = ( sum_{k>=1} -k^p * log2( 1 - 1/(k+1)^2 ) )^(1/p),
```

finite exactly when `p < 1`: the harmonic-mean case is
`K_{-1} = 1.74540566240...`, while the arithmetic mean (`p = 1`) diverges
almost everywhere, as section 2 noted. A closely related result is the
**Lévy (Khinchin–Lévy) constant**: the denominators `q_n` of the convergents
of almost every `x` satisfy

```
q_n^(1/n)  ->  e^(pi^2 / (12 ln 2))  =  3.275822918721811...
```

### 4.5 Dilogarithm and integral forms

The Wikipedia article collects further closed forms. Two of them, both
re-verified numerically to 40 digits while preparing this document:

```
ln(K0/2) = (1/ln 2) * [ Li2(-1/2) + (1/2) sum_{k=2}^{oo} (-1)^k Li2(4/k^2) ]
```

with `Li2` the dilogarithm (the `k = 2` term contributes
`Li2(1) = pi^2/6`), and the Gauss-measure integral

```
Integral_{0}^{1} log2(floor(1/x)) / (1 + x) dx  =  ln K0,
```

which is just the Birkhoff integral of section 2 restated — its
`ln`-numerator variant without the `1/ln 2` normalization is OEIS A247038
`= ln 2 * ln K0` (section 7). Elegant, but neither is computationally
competitive: the dilogarithm series converges only like `1/k^2`, and the
integrand jumps at every `x = 1/k`, which also defeats naive numerical
quadrature.

## 5. Implementations

- **This program (`khinchin_fast.c`)** — the accelerated series of 4.3 with a
  two-region evaluation: FLINT's reverse Bernoulli iterator below
  `n_direct ~ P/(2(log2 N + 3))`, the literal low-precision tail sum above it,
  per-entry dropping precision in all power tables, cost-model load balancing
  over dynamically scheduled blocks, and proof-carrying decimal rounding
  (output digits are only printed when Arb certifies the scaled integer is
  unique). 200,000 digits take under a minute on a 24-thread laptop; the
  measured cost curve is `O(d^2.01)`.
- **FLINT/Arb** — `arb_const_khinchin()` implements the same series family
  with rigorous balls; this project's `KHINCHIN_BACKEND=arb` mode is an
  independent implementation in that style, kept as a cross-check oracle.
- **mpmath** (Python) — `mpmath.khinchin` uses the plain zeta series with
  no small-k acceleration (`khinchin_fixed` in `libmp/gammazeta.py`), so it
  needs ~P/2 full-precision Bernoulli terms; fine up to about a thousand
  digits. `ports/khinchin.py` in this repository applies the acceleration
  of 4.3 in the same fixed-point style and is a measured 30–37x faster at
  1000–2000 digits.
- **Mathematica** — `Khinchin` is a built-in constant with
  arbitrary-precision evaluation, but a slow one (508 s for 10,000 digits
  with Mathematica 14.2 on this machine); `ports/khinchin.wl` implements
  the accelerated series of 4.3 on top of Mathematica's symbolic
  `Zeta[2n]` and is 12.7–23x faster.
- **PARI/GP** — no built-in as of 2.18 (the OEIS entry's PARI program was
  removed in 2010); the `README.md` benchmark uses `ports/khinchin.gp`,
  which implements the same accelerated series (serial and `parvector`
  parallel). This program is 21–56x faster than the parallel GP script on the
  same machine, while producing identical digits.
- **y-cruncher** — does not support Khinchin's constant, because no
  binary-splitting-friendly series is known (custom-formula finite-product
  demos are approximations, not full-precision computations).

Beyond the Python and GP ports mentioned above, this repository's `ports/`
directory carries the accelerated series of 4.3 in five more systems:
Julia (threaded `BigFloat`), Rust (`rug` with rayon), Fortran (GNU MPFR
bound through `ISO_C_BINDING`, OpenMP-parallel), Maple, and Mathematica
(accelerated series over symbolic `Zeta[2n]`, 12.7–23x faster than the
built-in constant).
The Julia, Rust, and Fortran ports obtain the even zeta values from the
classical positive-term recurrence
`(n + 1/2) zeta(2n) = sum_{j=1}^{n-1} zeta(2j) zeta(2n-2j)` instead of
Bernoulli numbers — numerically benign, but `O(M^2)` full-precision
products, which bounds their speed. Every executable port takes
`DIGITS [OUTPUT_FILE]`, and all locally testable ones write output files
byte-identical to the C program's; `README.md` carries the measured
cross-language speed table.

## 6. Records

Digit-computation records for `K0` (published, widely cited):

| Year | Digits | Who / method |
|---:|---:|---|
| 1959–1960s | tens of digits | Shanks, Wrench and others, direct product/early series |
| 1997 | 7,350 | Bailey, Borwein, Crandall — zeta-series acceleration |
| 1998 | 110,000 | Gourdon — refined acceleration |
| 2016 | 1,000,000 | Carles Simó (University of Barcelona) — PARI, ~12 days, single core; announced through OEIS A002210, where a cached copy of the digits is kept with permission |

For perspective, this machine reproduces the 1998 record in about 15 seconds
and the 2016 million-digit computation in 57 minutes (measured). But because every known
algorithm costs `~d^2/polylog`, each 10x in digits costs ~100x in time:
million-digit computations are now commodity, while 10^9 digits would take
decades — versus minutes for `pi`, which has quasi-linear algorithms. The
fundamental open problem: find a single series with rational terms (amenable
to binary splitting) whose limit is `K0`, or prove none exists. `K0` is a sum
of infinitely many algebraically independent transcendentals
(`zeta(2), zeta(4), ...`), each needed to full precision, and nobody has
found a way around that. It is not even known whether `K0` is irrational.

## 7. The OEIS entries (A002210 and friends)

The decimal expansion of `K0` is entry **A002210** in the On-Line Encyclopedia
of Integer Sequences ("Decimal expansion of Khinchin's constant", keyword
`nice`), entered by N. J. A. Sloane; it predates the online OEIS, having
appeared in his 1973 *Handbook of Integer Sequences* and the 1995
*Encyclopedia of Integer Sequences* (ids N0609 and M1564). The constant is
named after the Soviet mathematician Aleksandr Yakovlevich Khinchin
(1894–1959); the OEIS normalized the spelling "Khinchin" (over the older
transliteration "Khintchine") in 2024.

Facts recorded on the entry, cross-checked against this program:

- The 104 digits in the entry's data field match `khinchin-1m.txt` digit for
  digit (verified as part of this project's review).
- **Digit provenance**: Harry J. Smith's b-file supplies 1,200 digits; Simon
  Plouffe hosts 1,024-digit and 110,000-digit tables; and Carles Simó's 2016
  computation of 10^6 digits (reported to the OEIS in October 2016) is cached
  on the entry with permission — the record this program reproduces in under
  an hour.
- **Formulas**: the defining product `prod_{k>=1} (1 + 1/(k(k+2)))^(log2 k)`,
  and `K0 = exp(A247038 / log 2)`, where **A247038** is the decimal expansion
  of `Integral_{x=0..1} ln(floor(1/x)) / (1+x) dx = 0.6847247885631571...`.
  That integral is precisely the unnormalized Birkhoff integral of section 2
  (the Gauss measure carries an extra `1/ln 2`), i.e. `ln(2) * ln(K0)` — the
  exact quantity this program's series computes before the final division by
  `ln 2` and exponentiation.
- **A002211** is the continued fraction expansion of `K0` itself:
  `[2; 1, 2, 5, 1, 1, 2, 1, 1, 3, 10, ...]`. Its digit statistics feed the
  self-referential numerical evidence of section 1 that `K0` obeys its own
  theorem.
- **Further references carried by the entry**: S. R. Finch, *Mathematical
  Constants* (Cambridge, 2003), pp. 59–65; I. Vardi, *Computational
  Recreations in Mathematica* (1991), p. 164; F. Le Lionnais, *Les Nombres
  Remarquables* (1983), p. 46; Ph. Flajolet and I. Vardi, "Zeta function
  expansions of some classical constants"; T. Wieting, "A Khinchin Sequence",
  *Proc. Amer. Math. Soc.* 136 (2008); Khinchin's original Compositio
  Mathematica papers (1935, 1936); and a 2013 Numberphile video ("Six
  Sequences") featuring the constant.

## 8. References

- A. Ya. Khinchin, *Continued Fractions*, 3rd ed., University of Chicago
  Press, 1964 (original Russian 1935).
- D. H. Bailey, J. M. Borwein, R. E. Crandall, "On the Khintchine constant",
  *Mathematics of Computation* 66 (1997), 417–431.
- C. Ryll-Nardzewski, "On the ergodic theorems II (Ergodic theory of continued
  fractions)", *Studia Mathematica* 12 (1951), 74–79.
- R. O. Kuzmin, "Sur un problème de Gauss", *Atti Congr. Intern. Bologne* 6
  (1928); P. Lévy, *Bull. Soc. Math. France* 57 (1929).
- X. Gourdon, P. Sebah, "Numbers, constants and computation" (web resource) —
  Khinchin's constant notes and the 110,000-digit computation.
- R. Bloemen, `https://github.com/recmo/khinchin` — reverse-Bernoulli
  acceleration this project adapts (see `NOTICE.md`).
- F. Johansson, FLINT/Arb documentation: `arb_const_khinchin`,
  `bernoulli_rev` — `https://flintlib.org`.
- S. R. Finch, *Mathematical Constants*, Encyclopedia of Mathematics and its
  Applications 94, Cambridge University Press, 2003, pp. 59–65.
- D. Shanks, J. W. Wrench, Jr., "Khintchine's constant", *Amer. Math.
  Monthly* 66 (1959), 276–279; J. W. Wrench, "Further evaluation of
  Khintchine's constant", *Math. Comp.* 14 (1960), 370–371.
- OEIS Foundation, entries A002210 (decimal expansion), A002211 (continued
  fraction), A247038 (`ln 2 * ln K0`) — `https://oeis.org/A002210`.
- Wikipedia, "Khinchin's constant" —
  `https://en.wikipedia.org/wiki/Khinchin%27s_constant` (dilogarithm and
  integral forms, Hölder-mean formulas, empirical-typicality lists).
