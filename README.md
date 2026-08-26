# khinchin-fast

Computes Khinchin's constant using a parallel Bailey-Borwein-Crandall
acceleration and an SCP-style reverse Bernoulli iterator implemented with
FLINT/Arb. FLINT uses GMP for large integer arithmetic and Arb tracks a rigorous
interval; the program writes a correctly rounded decimal expansion.

The `ports/` directory carries implementations of the same accelerated
series for Python (30-37x faster than mpmath's built-in), Julia, Rust,
PARI/GP, Fortran, Maple, and Mathematica — see the reference-ports
section below.

## Build

```sh
make
```

Fedora packages required: `gcc`, `flint-devel`, `mpfr-devel`, and `gmp-devel`.

## Use

The number is the requested count of digits after the decimal point:

```sh
./khinchin-fast 1000000 khinchin-1m.txt
```

Benchmark calculation without converting or writing the digits:

```sh
./khinchin-fast --benchmark 1000000
```

The reverse Bernoulli backend is the default. The previous independent Arb
backend remains available for cross-checking:

```sh
KHINCHIN_BACKEND=arb ./khinchin-fast 10000 khinchin-arb-10k.txt
```

Timing statistics are printed to standard error. The result file contains only
the decimal value and a final newline. A conservative memory preflight rejects
jobs whose estimated peak exceeds 80% of physical RAM (see the memory section
below for the model and its calibration); `--force` overrides this check.

Tuning and diagnostic environment variables for the default backend:
`KHINCHIN_N` overrides the power-table cut-off, `KHINCHIN_DIRECT` overrides
the split between the direct and Bernoulli regions, and `KHINCHIN_VERBOSE=1`
prints per-block ranges and timings.

## Benchmark on this machine

Intel Core Ultra 9 275HX, 24 threads, FLINT 3.4.0, August 25, 2026:

| Digits after decimal | This version | Previous SCP | Earlier parallel Arb | Peak RAM |
|---:|---:|---:|---:|---:|
| 10,000 | 0.14 s | 0.34 s | 0.64 s | 37 MB |
| 20,000 | 0.43 s | 1.21 s | 3.47 s | 54 MB |
| 30,000 | 0.98 s | 2.64 s | 8.58 s | 77 MB |
| 50,000 | 2.67 s | 7.89 s | 34.61 s | 140 MB |
| 100,000 | 12.09 s | 37.40 s | not run | 342 MB |
| 200,000 | 55.75 s | 170.82 s | not run | 937 MB |
| 1,000,000 | 3413.29 s | not run | not run | 16.1 GB |

A log-log fit to the 10k-200k measurements is `time = O(digits^2.01)`, close
to the theoretical quadratic behaviour of the series. Beyond 200k the local
exponent steepens (about 2.5 between 200k and 1M) as the working set outgrows
the caches and the reverse Bernoulli iterators' quadratically growing state
starts to dominate: one million digits took 56.9 minutes, not the ~24 the
small-size fit projects. For context, the largest documented computation of
Khinchin's constant was 10^6 digits in about 12 days of single-core PARI
(University of Barcelona, 2016); this program reproduces it in under an hour.
Extrapolating to one billion digits still gives centuries: the exponent, not
the constant factor, is the obstruction (see the next sections), and memory
(16.1 GB peak at 1M, growing superlinearly) walls off the current in-memory
design well before time does — see the next section.

## Memory: the real limitation

Beyond roughly 200k digits, peak memory is dominated not by the shared power
tables (which grow about linearly with precision) but by the internal state of
FLINT's reverse Bernoulli iterators. An iterator started at `s = 2*n_direct`
stores about `2^(p/s)` scaled zeta powers whose per-entry precisions drop from
`p ~ s*(log2 s - 4.09)` bits — roughly quadratic growth in the digit count —
and about half the threads hold near-top iterators simultaneously, because the
heaviest blocks are scheduled first. Measured peak RSS (`/usr/bin/time -v`)
against the program's preflight estimate:

| Digits after decimal | Preflight estimate | Measured peak RSS |
|---:|---:|---:|
| 10,000 | 80 MiB | 37 MiB |
| 100,000 | 531 MiB | 337 MiB |
| 200,000 | 1,439 MiB | 947 MiB |
| 300,000 | 2,679 MiB | 1,957 MiB |
| 1,000,000 | 18.7 GiB | ~15.4 GiB |

The preflight models both allocations explicitly and is calibrated on the
measurements above to stay 1.2-1.6x over the actual peak — deliberately
conservative, since its job is to refuse jobs that would exhaust RAM. (An
earlier version modelled only the power tables and under-estimated the
million-digit run by more than 2x, which would have let a ~2.5M-digit job
sail into an out-of-memory crash.) On this 64 GB machine the preflight starts
refusing around 2 million digits, and the measured growth curve suggests such
a job would genuinely need ~55 GB: memory, not time, is the first hard wall
of the current all-in-memory design. Breaking past it would need out-of-core
Bernoulli generation or slice-wise recomputation of the zeta tails — both of
which give back part of the constant-factor speed this program exists to win.

## Comparison with PARI/GP

PARI/GP (2.18.1 development, multithreaded build, 24 threads) has no built-in
Khinchin constant, so the comparison uses a GP implementation of the same
accelerated series: a serial version with the incremental power table, and a
parallel version that splits the zeta range across threads with `parvector`.
The script lives at `ports/khinchin.gp`; the `psi`-based block seeding used
when this table was measured has since been replaced by direct summation
(~4x faster at 10k digits — current timings are in the by-language table
below). Wall-clock times as measured:

| Digits after decimal | GP serial | GP parallel (24 threads) | This program | Speedup vs GP parallel |
|---:|---:|---:|---:|---:|
| 10,000 | 4.62 s | 3.08 s | 0.14 s | 21x |
| 20,000 | not run | 16.72 s | 0.43 s | 39x |
| 50,000 | not run | 148.48 s | 2.67 s | 56x |

The GP results agree with this program's output digit for digit, which serves
as an additional independent cross-check (different zeta implementation,
different bignum stack above GMP).

## Reference ports in other languages

`ports/` holds implementations of the same accelerated series for other
systems. Every executable port takes `DIGITS [OUTPUT_FILE]` and, given a
file, writes the same decimal-plus-newline payload as the C program — all
seven produce byte-identical output files on this machine.

- `ports/khinchin.py` — fixed-point port on top of mpmath's internals.
  Measured 30x (1000 digits) to 37x (2000 digits) faster than mpmath's
  built-in `mp.khinchin` (the program on the OEIS entry), which uses the
  unaccelerated series. Deliberately serial: profiling shows 84% of the
  runtime inside mpmath's `mpf_bernoulli`, a strictly sequential cached
  recurrence, capping any multiprocess split at ~1.2x before IPC costs —
  and the GIL rules out threads. Installing gmpy2 switches mpmath's
  bignum backend to GMP automatically and is worth a further 3-7x (see
  the table below).
- `ports/khinchin.gp` — the serial and parallel (`parvector`) PARI/GP
  scripts behind the comparison table above. Each parallel block seeds its
  alternating-harmonic weight `h(first)` by direct summation, which
  measured ~4x faster at 10k digits than seeding with `psi`.
- `ports/khinchin.jl` — Julia on stock `BigFloat` (MPFR), no packages,
  multithreaded (run with `julia -t auto`): 6.3x over serial at 10k
  digits.
- `ports/khinchin-rs/` — Rust on `rug` (GMP/MPFR) with rayon: 13x over
  serial at 10k digits. `cargo run --release -- 1000 out.txt`.
- `ports/khinchin.f90` — Fortran bound directly to GNU MPFR through
  standard `ISO_C_BINDING` (Fortran itself has no arbitrary-precision
  arithmetic), OpenMP-parallel: 9.1x over serial at 10k digits.
- `ports/khinchin.mpl` — Maple version of the same series, serial.
  Verified with Maple 2024.2: byte-identical output to the C program at
  1000 digits, digit-exact against the reference file at 10,000.
- `ports/khinchin.wl` — Mathematica, serial. Implements the accelerated
  series directly (Mathematica evaluates `Zeta[2n]` symbolically through
  Bernoulli numbers), measured 12.7x faster than the built-in `Khinchin`
  constant at 1000 digits and 23x at 10,000 with Mathematica 14.2; the
  built-in — what the OEIS entry's program uses — is kept as
  `KhinchinBuiltin` for cross-checking. Byte-identical output to the C
  program at 1000 digits.

The Julia, Rust, and Fortran ports take their even zeta values from the
classical positive-term recurrence
`(n + 1/2) zeta(2n) = sum_j zeta(2j) zeta(2n-2j)` — numerically benign
(all terms positive, roughly linear error growth) and free of Bernoulli
numbers, at the price of `O(M^2)` full-precision products, which is what
dominates their runtimes beyond a few thousand digits. In all of them the
recurrence convolutions (halved via their `j <-> n-j` symmetry) and the
blocked main loop run in parallel.

### Measured speed by language

Same machine as the benchmarks above (24 threads), computing 1,000 and
10,000 digits after the decimal point; all outputs byte-identical:

| Implementation | Parallelism | 1,000 digits | 10,000 digits |
|---|---|---:|---:|
| C — `khinchin-fast` (this repo) | OpenMP | 0.007 s | 0.15 s |
| PARI/GP — `ports/khinchin.gp` | parvector | 0.05 s | 2.8 s |
| Rust — `ports/khinchin-rs` | rayon | 0.04 s | 9.8 s |
| Fortran + MPFR — `ports/khinchin.f90` | OpenMP | 0.05 s | 12.3 s |
| Julia — `ports/khinchin.jl` | threads | 0.22 s | 19.2 s |
| Python — `ports/khinchin.py` | serial | 0.16 s | 64.5 s |
| Python + gmpy2 — same file, GMP backend | serial | 0.05 s | 9.0 s |
| Mathematica — `ports/khinchin.wl` | serial | 0.12 s | 21.7 s |
| Maple — `ports/khinchin.mpl` | serial | 8.5 s | 298 s |
| mpmath built-in `mp.khinchin` | serial | 4.67 s | not run |
| Mathematica built-in `Khinchin` | serial | 0.57 s | 508 s |

The C program leads by roughly 20-440x at 10k digits: it is the only
implementation with per-term dropping precision, the two-region split,
and FLINT's reverse Bernoulli iterator. PARI/GP places second because its
`zeta(2n)` rides PARI's fast Bernoulli machinery, while the
Julia/Rust/Fortran ports pay the `O(M^2)` zeta recurrence, Python adds
mpmath's pure-Python bignum layer on a serial design (with gmpy2's GMP
backend the identical file jumps ahead of the Rust and Fortran ports),
and Maple pays its software-float `evalf` layer serially. Mathematica's symbolic `Zeta[2n]`
(exact Bernoulli rationals) places its port between Julia and Python
despite being serial — and 23x ahead of its own built-in at 10k digits.
Julia's time excludes JIT warmup.

## Why this is the fastest known approach

The claim has three layers — the formula, the zeta strategy inside it, and
the constant factors — and each was measured against its alternatives
rather than assumed.

**The formula.** Every known expression for `K0` was considered. The
defining product converges far too slowly for even tens of digits. The
dilogarithm and integral closed forms (`mathematical_background.md`,
section 4.5) converge only like `1/k^2` — verified numerically, and
useless for digit computation. The plain BBC zeta series needs `~P/2`
full-precision zeta terms; the accelerated form used here, with the small
`k` pulled out of every zeta tail, needs only `~P/(2 log2 N)`. That gap
alone is measurable: mpmath's built-in uses the plain series, and
`ports/khinchin.py`'s switch to the accelerated one — same language, same
bignum layer — is worth 30-37x. Gourdon/Sebah-style rearrangements are
the same family with different constants. The only thing that would beat
this family is a single series with rational terms amenable to binary
splitting — quasi-linear, like the formulas behind pi records — and no
such series for `K0` is known; finding one is an open problem (next
section). Within current mathematical knowledge, the accelerated zeta
series is the end of the line.

**The zeta values.** The accelerated series stands or falls on how fast
the `zeta(2n)` values arrive, and the ports in this repository ended up
being a controlled experiment across every practical strategy — same
machine, identical output digits:

| zeta(2n) strategy | Where | 10,000 digits |
|---|---|---:|
| Reverse Bernoulli iterator + direct low-precision tail region | C, this program | 0.15 s |
| `arb_zeta_ui` per term at full precision | `KHINCHIN_BACKEND=arb` | 0.64 s, and scaling worse (34.6 s vs 2.67 s at 50k) |
| PARI's Bernoulli machinery | `ports/khinchin.gp` | 2.8 s |
| Positive-term recurrence, `O(M^2)` | Julia/Rust/Fortran ports | 10-19 s |
| Symbolic exact `Zeta[2n]` (Bernoulli rationals) | Mathematica port | 21.7 s |

**The constant factors.** On top of the winning strategy, this program
adds per-term dropping precision (a term of magnitude `N^(-2n)` only
needs about `P - 2n log2 N` accurate bits, so most arithmetic runs far
below full precision), the two-region split that eliminates Bernoulli
numbers entirely where the exact `B_2n` would carry more bits than the
working precision itself, and cost-model load balancing across threads.
That layer is why it runs 28x faster than FLINT's own
`arb_const_khinchin()` and reproduced the 2016 million-digit record
roughly 300x faster (under an hour against ~12 single-core days).

What remains is constant-factor tuning at best — for example, whether
batch multimodular Bernoulli generation could edge out the reverse
iterators at very large sizes. The genuine walls are the quadratic
exponent, which is a property of the mathematics (next section), and the
memory growth documented above — neither yields to engineering.

## Known mathematical approaches (all ~quadratic)

Every published way to compute Khinchin's constant to high precision reduces
to summing many high-precision zeta-related terms:

1. **Bailey-Borwein-Crandall zeta acceleration (1997)** - what this program
   uses. Needs `M = P/(2 log2 N)` terms, each requiring the corresponding
   `zeta(2n)` to nearly full absolute precision. Cost about `d^2/polylog`.
2. **Gourdon/Sebah-style variants** - same family, different acceleration
   constants and series rearrangements. Same asymptotics, only constant-factor
   differences.
3. **Direct definition** - the product over continued-fraction densities,
   `prod (1 + 1/(k(k+2)))^(log2 k)` - converges far too slowly; the BBC series
   is the accelerated form of this.
4. **Bernoulli-number route** (the recmo/SCP trick this program runs) -
   replaces expensive independent zeta evaluations with a reverse Bernoulli
   recurrence. Big constant-factor win, same exponent.

The structural obstruction: for pi or log 2, the digits come from a single
hypergeometric series with rational terms, which binary splitting evaluates in
`O(M(d) log d)` - quasi-linear. Khinchin's constant is a sum of infinitely
many algebraically independent transcendentals (`zeta(2)`, `zeta(4)`, ...),
each needed to full precision. No one has found a single closed-form series
with rational terms for it. That is not a gap in engineering effort - it is an
open problem in mathematics. If you found such a series, that would be a
publishable result in its own right; it is the only thing that would make
10^9 digits feasible.

## Implementation notes

Term `n` of the accelerated sum has magnitude about `N^(-2n)`, so it only
needs about `P - 2n log2 N` accurate bits, and

```
zeta(2n) - 1 - sum_{2<=k<N} k^(-2n)  =  sum_{k>=N} k^(-2n).
```

The SCP backend splits the zeta range at `n_direct ~ P/(2(log2 N + 3))`:

- **Direct region** (`n >= n_direct`, about the top third of the range): the
  term is computed as the literal tail sum over `k` in `[N, K]` with
  `K ~ 2^(P/2n) <= 8N`, entirely at low precision - no Bernoulli numbers, no
  cancellation, no full-precision arithmetic. A rigorous integral bound covers
  the omitted `k > K`. Before this split, profiling showed the reverse
  Bernoulli iterator dominating exactly here, because the exact `B_2n` carry
  about `2n log2 n` bits - more than the working precision itself near the top
  of the range.
- **Bernoulli region** (`n < n_direct`): FLINT's reverse Bernoulli iterator
  reconstructs `zeta(2n)`, which is near 1 and genuinely needs full precision
  there. Power-table entries and the tail assembly still run at their own
  dropping precisions.

Work is divided into about four blocks per thread, sized by an analytic cost
model, sorted heaviest first, and scheduled dynamically. Four blocks per
thread is the measured optimum: two per thread performs the same within noise,
while eight per thread is ~20% slower because every extra Bernoulli block pays
its own iterator initialisation. Each worker traverses
its block with a shared power table whose entries start at the precision their
final, most demanding use requires. The final interval includes the proven
truncation bounds, and decimal rounding succeeds only when Arb proves the
scaled integer is unique.

The backend was adapted from the algorithmic approach in Remco Bloemen's MIT
licensed `recmo/khinchin` work and FLINT's reverse Bernoulli implementation. See
`NOTICE.md` for source links and licenses.
