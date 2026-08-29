# Two Million Digits of Khinchin's Constant

**Reza K Ghazi** ([ORCID 0009-0009-3286-4675](https://orcid.org/0009-0009-3286-4675)) · August 27, 2026

*A record computation of Khinchin's constant to 2,000,000 decimal
digits — rigorously certified, independently cross-checked, and
completed in under nine hours on a single desktop machine. Code,
digits, and raw telemetry:
<https://github.com/reza-ghazi/khinchin-fast>.*

---

## 1. The constant

Write a real number `x` as a simple continued fraction:

```
x = a0 + 1/(a1 + 1/(a2 + 1/(a3 + ...)))
```

In 1935 Aleksandr Khinchin proved something remarkable: for almost
every real number (all but a set of Lebesgue measure zero), the
geometric mean of the partial quotients `a1, a2, a3, ...` converges to
one universal constant, independent of the number chosen:

```
lim_{n->oo} (a1 * a2 * ... * an)^(1/n)  =  K0  =  2.68545200106530644530...
```

The constant has a closed product form,

```
K0 = prod_{k=1}^{oo} ( 1 + 1/(k(k+2)) )^(log2 k),
```

and is entry [A002210](https://oeis.org/A002210) in the OEIS. It sits
at the center of the metric theory of continued fractions, yet almost
nothing is known about it as a number: it is not known whether `K0` is
irrational, let alone transcendental — and, in a pleasant irony, it is
not known whether `K0`'s own continued fraction obeys Khinchin's
theorem (numerically, it appears to).

## 2. Why the digits are hard

For constants like `pi`, hypergeometric series with rational terms
admit binary splitting, giving quasi-linear algorithms — which is why
`pi` is known to hundreds of trillions of digits. No such series is
known for `K0`. Every practical method goes through the even zeta
values: taking logarithms of the product and expanding gives the
Bailey–Borwein–Crandall series

```
ln(2) * ln(K0) = sum_{n=1}^{oo} (zeta(2n) - 1)/n * h(n),
      h(n) = sum_{j=1}^{2n-1} (-1)^(j+1) / j,
```

a sum over infinitely many algebraically independent transcendentals
`zeta(2), zeta(4), ...`, each needed essentially to full precision.
Nobody has found a way around that, and every known algorithm costs on
the order of `d^2` (up to logarithmic factors) for `d` digits. Each
10x in digits costs ~100x in time. That quadratic wall is why the
record history is so short:

| Year | Digits | Who / method |
|---:|---:|---|
| 1959–1960s | tens | Shanks, Wrench and others |
| 1997 | 7,350 | Bailey, Borwein, Crandall — accelerated zeta series |
| 1998 | 110,000 | Gourdon — refined acceleration |
| 2016 | 1,000,000 | Simó — PARI, ~12 days, single core |
| **2026** | **2,000,000** | **this computation — 8 h 47 m, one desktop machine** |

## 3. The algorithm

The program (`khinchin_fast.c`, a single C file on FLINT/Arb ball
arithmetic with OpenMP) evaluates the accelerated form of the BBC
series: for a cut-off `N`,

```
ln(2) * ln(K0) = - sum_{k=2}^{N-1} ln((k-1)/k) * ln((k+1)/k)
                 + sum_{n=1}^{oo} ( zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n) ) / n * h(n),
```

whose `n`-th term now shrinks like `N^(-2n)`, so about
`P / (2 log2 N)` terms suffice for `P` bits. Three further ideas carry
most of the speed:

- **Dropping precision.** A term of magnitude `N^(-2n)` contributes at
  most `P - 2n*log2(N)` bits to the final sum, so almost all
  arithmetic runs far below full precision. The shared power tables
  drop precision per entry the same way.
- **A two-region split.** For `n` above `~P/(2(log2 N + 3))` the term
  `zeta(2n) - 1 - sum k^(-2n)` is computed as the literal tail sum
  `sum_{k>=N} k^(-2n)` at low precision — no Bernoulli numbers at all.
  Below the split, where `zeta(2n)` is genuinely near 1 and needs full
  precision, FLINT's reverse Bernoulli iterator (`bernoulli_rev`)
  reconstructs it. Work is split into cost-balanced blocks scheduled
  dynamically across all threads.
- **Rigor by construction.** All arithmetic is Arb ball arithmetic:
  every intermediate is an interval guaranteed to contain the true
  value, every truncation is folded in as an explicit error bound, and
  the final decimal is produced by a proof-carrying rounding step that
  only prints digits when the ball certifies a unique correctly
  rounded value. There is no "estimated error" anywhere — a wrong
  digit is impossible; failure would be loud, not silent.

The zeta-acceleration structure is adapted from Remco Bloemen's
MIT-licensed `recmo/khinchin`; the repository's `MATHEMATICS.md` gives
the full derivation and history.

## 4. The record run

On August 27, 2026, the program computed 2,000,000 digits after the
decimal point in one run:

| | |
|---|---|
| Hardware | Intel Core Ultra 9 275HX (24 threads), 64 GB RAM |
| Software | Fedora Linux 44, GCC 16.2, FLINT 3.4.0, stock GMP/MPFR |
| Working precision | 6,643,988 bits |
| Wall time | 8 h 47 m 26 s (compute 31,644 s) |
| Average CPU | 931% (memory-bandwidth-bound at this size) |
| Peak memory | 55.8 GiB resident |

The run pushed the machine to its edge: the program's own conservative
memory preflight refuses 2,000,000 digits on a 64 GB machine
(estimate: 61.7 GiB) and had to be overridden with `--force`; measured
peak was 55.8 GiB, within 6 GiB of physical RAM, and a standby
swapfile went unused. Memory, not time, is the binding constraint —
the dominant allocation is the internal state of the reverse Bernoulli
iterators, which grows roughly quadratically with the digit count. The
local cost exponent, `~d^2.0` for small sizes, steepens to `~3.2`
between one and two million digits as that state overwhelms the
caches. The next doubling would need well over 128 GB of RAM and an
out-of-core redesign, and a billion digits remains out of reach for
everyone: the exponent, not the constant factor, is the obstruction.

## 5. Verification

Three independent legs support the result:

1. **Interval certification.** As described above, the computation
   carries a rigorous enclosure end to end; the printed digits are
   certified correctly rounded by construction, not checked after the
   fact.
2. **An independent backend.** The program ships a second,
   structurally different backend (per-term `arb_zeta_ui` at full
   precision) sharing only the outer series framing. The two agree
   byte-for-byte at every size tested.
3. **Agreement with the previous record.** The first 1,000,000 digits
   of the new result are identical, digit for digit, to Carles Simó's
   independent 2016 computation (different program, different
   algorithm details, different CPU era), as cached on OEIS A002210.

## 6. Availability

Everything is in the repository
<https://github.com/reza-ghazi/khinchin-fast> (MIT license): the C
program, the 2,000,000-digit result (`khinchin-2m.txt`), the raw run
telemetry (`khinchin-2m.log`), the mathematical background document,
and — as a byproduct of validating the algorithm — reference
implementations of the same accelerated series in seventeen further
languages, from FLINT-backed C++/Rust/Julia/Fortran to pure-bignum
Go/Java/Haskell/OCaml, every one verified byte-identical to the C
program's output. A single `make check` rebuilds and re-verifies the
whole set. The repository is permanently archived at Zenodo:
[doi:10.5281/zenodo.22134905](https://doi.org/10.5281/zenodo.22134905).

## Acknowledgments

The software and this write-up were developed with the assistance of
Claude (Anthropic), used as a programming and drafting tool under the
author's direction; all results are independently verifiable as
described in section 5.

## References

- A. Ya. Khinchin, *Metrische Kettenbruchprobleme*, Compositio
  Mathematica 1 (1935), 361–382; *Continued Fractions*, Dover reprint
  1997.
- D. H. Bailey, J. M. Borwein, R. E. Crandall, *On the Khintchine
  constant*, Mathematics of Computation 66 (1997), 417–431.
- X. Gourdon, P. Sebah, *The Khintchine constant* (Numbers, constants
  and computation), 1998–2002.
- C. Simó, *Computation of 10^6 digits of Khintchine's constant*,
  2016, announced via OEIS A002210.
- OEIS Foundation, entry A002210, <https://oeis.org/A002210>.
- F. Johansson, *Arb: efficient arbitrary-precision midpoint-radius
  interval arithmetic*, IEEE Transactions on Computers 66 (2017),
  1281–1292; FLINT, <https://flintlib.org>.
- R. Bloemen, `recmo/khinchin`, 2022,
  <https://github.com/recmo/khinchin>.
- R. K Ghazi, *khinchin-fast: Khinchin's constant to 2,000,000 certified
  decimal digits*, Zenodo, 2026,
  [doi:10.5281/zenodo.22134905](https://doi.org/10.5281/zenodo.22134905).
