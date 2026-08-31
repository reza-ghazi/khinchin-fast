<!--
Medium-ready rendition of ARTICLE.md. ARTICLE.md remains the canonical
write-up and is unchanged; this file exists only because Medium's editor
has no tables and no math rendering.

HOW TO USE
  1. https://medium.com/p/import — paste
     https://reza-ghazi.github.io/khinchin-fast/
     This is what sets rel=canonical back to the original; a plain paste
     cannot. Do NOT use the claude.ai artifact URL: it serves only a
     wrapper to outside fetchers, so the import fails. Let the import
     finish even though the draft will be messy.
  2. Delete the imported body, then paste this file's body over it.
     The canonical link stays attached to the draft.
  3. Handle the two [IMAGE] markers below (screenshot from the artifact
     page), then delete every HTML comment and marker line.
  4. Type ``` on an empty line in Medium to open a code block; paste the
     fenced content in, without the backticks.
-->

# Two Million Digits of Khinchin's Constant

<!-- Medium subtitle: paste the next line, select it, and use the "small T" -->
A record computation, rigorously certified, finished overnight on one desktop machine.

Almost every real number hides the same number inside it. Not most numbers in a loose sense — *almost every* one, in the precise measure-theoretic sense, with the exceptions forming a set of measure zero. Write down a real number at random and take the continued fraction expansion, and the geometric mean of its terms marches toward one universal constant that has nothing to do with the number you started with.

That constant is 2.68545200106530644530…, and on August 27, 2026 I computed two million of its digits — twice as many as anyone had before, in 8 hours 47 minutes on a desktop machine, with every printed digit mathematically certified rather than merely believed.

Here is how, and why it is harder than it sounds.

## 1. The constant

Write a real number `x` as a simple continued fraction:

```
x = a0 + 1/(a1 + 1/(a2 + 1/(a3 + ...)))
```

In 1935 Aleksandr Khinchin proved something remarkable: for almost every real number — all but a set of Lebesgue measure zero — the geometric mean of the partial quotients `a1, a2, a3, …` converges to one universal constant, independent of the number chosen:

```
lim (a1 · a2 · … · an)^(1/n)  =  K0  =  2.68545200106530644530…
n→∞
```

The constant has a closed product form,

```
K0 = ∏ ( 1 + 1/(k(k+2)) )^(log2 k),   k = 1, 2, 3, …
```

and is entry [A002210](https://oeis.org/A002210) in the On-Line Encyclopedia of Integer Sequences. It sits at the center of the metric theory of continued fractions, and yet almost nothing is known about it as a number. It is not known whether `K0` is irrational, let alone transcendental. And in a pleasant irony, it is not known whether `K0`'s own continued fraction obeys Khinchin's theorem — numerically, it appears to.

## 2. Why the digits are hard

For constants like π, hypergeometric series with rational terms admit binary splitting, giving quasi-linear algorithms. That is why π is known to hundreds of trillions of digits.

No such series is known for `K0`. Every practical method goes through the even zeta values. Taking logarithms of the product and expanding gives the Bailey–Borwein–Crandall series:

```
ln(2) · ln(K0) = sum_{n=1}^{oo} (zeta(2n) - 1)/n · h(n),

       h(n) = sum_{j=1}^{2n-1} (-1)^(j+1) / j
```

That is a sum over infinitely many algebraically independent transcendentals `zeta(2), zeta(4), …`, each needed essentially to full precision. Nobody has found a way around it, and every known algorithm costs on the order of `d²` (up to logarithmic factors) for `d` digits. Each 10× in digits costs roughly 100× in time.

That quadratic wall is why the record history is so short — five entries in seventy years:

<!-- [IMAGE] Upload paper/medium-assets/records-table.png here (drag it into
     the Medium editor, or click the + and choose the image). Medium has no
     tables, so this is the table as a picture. Caption suggestion:
     "Seventy years of record computations. Each 10x in digits costs ~100x
     in time." If you prefer text, delete this comment and keep the list
     below; otherwise delete the list once the image is in. -->

- **1959–1960s** — Shanks, Wrench and others: tens of digits
- **1997** — Bailey, Borwein, Crandall: 7,350 digits, via the accelerated zeta series
- **1998** — Gourdon: 110,000 digits, refined acceleration
- **2016** — Simó: 1,000,000 digits, PARI, about 12 days on a single core
- **2026** — this computation: **2,000,000 digits, 8 h 47 m, one desktop machine**

## 3. The algorithm

The program is a single C file built on FLINT/Arb ball arithmetic with OpenMP. It evaluates an accelerated form of the BBC series: for a cut-off `N`,

```
ln(2) · ln(K0) = - sum_{k=2}^{N-1} ln((k-1)/k) · ln((k+1)/k)

                 + sum_{n=1}^{oo} ( zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n) ) / n · h(n)
```

whose `n`-th term now shrinks like `N^(-2n)`, so about `P / (2 log2 N)` terms suffice for `P` bits. Three further ideas carry most of the speed.

**Dropping precision.** A term of magnitude `N^(-2n)` contributes at most `P - 2n·log2(N)` bits to the final sum, so almost all arithmetic in the program runs far below full precision. The shared power tables drop precision per entry the same way. This sounds like a micro-optimization and is in fact the difference between hours and days.

**A two-region split.** For `n` above roughly `P/(2(log2 N + 3))`, the quantity `zeta(2n) - 1 - sum k^(-2n)` is exactly the tail `sum_{k>=N} k^(-2n)`, and can be summed literally at low precision — no Bernoulli numbers at all. Below the split, where `zeta(2n)` is genuinely near 1 and needs full precision, FLINT's reverse Bernoulli iterator reconstructs it. Work is divided into cost-balanced blocks scheduled dynamically across all threads.

**Rigor by construction.** All arithmetic is Arb ball arithmetic: every intermediate value is an interval guaranteed to contain the true value, every truncation is folded in as an explicit error bound, and the final decimal expansion comes from a proof-carrying rounding step that prints digits only when the enclosing ball certifies a unique correctly rounded result.

There is no "estimated error" anywhere in the pipeline. A wrong digit is not unlikely — it is impossible. Failure would be loud rather than silent, which is the property you want when nobody on earth can check your output by other means.

The zeta-acceleration structure is adapted from Remco Bloemen's MIT-licensed [recmo/khinchin](https://github.com/recmo/khinchin).

## 4. The record run

On August 27, 2026, the program computed 2,000,000 digits after the decimal point in a single run:

<!-- [IMAGE] Optional: screenshot the run-parameters table from the artifact
     page. Otherwise this list carries the same information. -->

- **Hardware** — Intel Core Ultra 9 275HX (24 threads), 64 GB RAM
- **Software** — Fedora Linux 44, GCC 16.2, FLINT 3.4.0, stock GMP/MPFR
- **Working precision** — 6,643,988 bits
- **Wall time** — 8 h 47 m 26 s (compute: 31,644 s)
- **Average CPU** — 931% of 2400%, memory-bandwidth-bound at this size
- **Peak memory** — 55.8 GiB resident

The run pushed the machine to its edge. The program's own conservative memory preflight *refuses* 2,000,000 digits on a 64 GB machine — its estimate was 61.7 GiB — and had to be overridden with `--force`. Measured peak came in at 55.8 GiB, within 6 GiB of physical RAM, and a standby swapfile went untouched.

Memory, not time, is the binding constraint. The dominant allocation is the internal state of the reverse Bernoulli iterators, which grows roughly quadratically with the digit count. The local cost exponent, about `d^2.0` at moderate sizes, steepens to around `d^3.2` between one and two million digits as that state overwhelms the caches.

The next doubling would need well over 128 GB of RAM and an out-of-core redesign. A billion digits remains out of reach for everyone: the exponent, not the constant factor, is the obstruction.

## 5. Verification

A record you cannot check is a rumor. Three independent legs support this one.

**Interval certification.** The computation carries a rigorous enclosure end to end. The printed digits are certified correctly rounded by construction, not checked after the fact.

**An independent backend.** The program ships a second, structurally different backend — per-term `arb_zeta_ui` at full precision — that shares only the outer series framing with the fast one. The two agree byte for byte at every size tested.

**Agreement with the previous record.** The first 1,000,000 digits of the new result are identical, digit for digit, to Carles Simó's independent 2016 computation: different program, different algorithmic details, different CPU era, as cached on OEIS A002210.

## 6. Availability

Everything is in the repository at [github.com/reza-ghazi/khinchin-fast](https://github.com/reza-ghazi/khinchin-fast) under the MIT license: the C program, the 2,000,000-digit result, the raw run telemetry, the mathematical background document, and — as a byproduct of validating the algorithm — reference implementations of the same accelerated series in seventeen further languages, from FLINT-backed C++, Rust, Julia and Fortran to pure-bignum Go, Java, Haskell and OCaml. Every one is verified byte-identical to the C program's output, and a single `make check` rebuilds and re-verifies the whole set.

The repository is permanently archived at Zenodo, [doi:10.5281/zenodo.22134905](https://doi.org/10.5281/zenodo.22134905). Since August 31, 2026 the computation is also listed in the links of the constant's OEIS entry, [A002210](https://oeis.org/A002210), alongside Bailey–Borwein–Crandall and Simó.

## Acknowledgments

The software and this write-up were developed with the assistance of Claude (Anthropic), used as a programming and drafting tool under my direction. All results are independently verifiable as described in section 5.

## References

- A. Ya. Khinchin, *Metrische Kettenbruchprobleme*, Compositio Mathematica 1 (1935), 361–382; *Continued Fractions*, Dover reprint 1997.
- D. H. Bailey, J. M. Borwein, R. E. Crandall, *On the Khintchine constant*, Mathematics of Computation 66 (1997), 417–431.
- X. Gourdon, P. Sebah, *The Khintchine constant* (Numbers, constants and computation), 1998–2002.
- C. Simó, *Computation of 10⁶ digits of Khintchine's constant*, 2016, announced via OEIS A002210.
- OEIS Foundation, entry [A002210](https://oeis.org/A002210).
- F. Johansson, *Arb: efficient arbitrary-precision midpoint-radius interval arithmetic*, IEEE Transactions on Computers 66 (2017), 1281–1292; FLINT, [flintlib.org](https://flintlib.org).
- R. Bloemen, [recmo/khinchin](https://github.com/recmo/khinchin), 2022.
- R. K Ghazi, *khinchin-fast: Khinchin's constant to 2,000,000 certified decimal digits*, Zenodo, 2026, [doi:10.5281/zenodo.22134905](https://doi.org/10.5281/zenodo.22134905).

<!--
SUGGESTED MEDIUM TAGS (5 max, in this order):
  Mathematics, Number Theory, Programming, Computer Science, Open Source

PULL QUOTE candidate (Medium: select the line, click the " button):
  "A wrong digit is not unlikely — it is impossible."
-->
