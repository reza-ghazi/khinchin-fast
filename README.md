# khinchin-fast

Computes Khinchin's constant using a parallel Bailey-Borwein-Crandall
acceleration and an SCP-style reverse Bernoulli iterator implemented with
FLINT/Arb. FLINT uses GMP for large integer arithmetic and Arb tracks a rigorous
interval; the program writes a correctly rounded decimal expansion.

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
parallel version that splits the zeta range across threads with `parvector`,
seeding each block's alternating harmonic weight with `psi`. Wall-clock times:

| Digits after decimal | GP serial | GP parallel (24 threads) | This program | Speedup vs GP parallel |
|---:|---:|---:|---:|---:|
| 10,000 | 4.62 s | 3.08 s | 0.14 s | 21x |
| 20,000 | not run | 16.72 s | 0.43 s | 39x |
| 50,000 | not run | 148.48 s | 2.67 s | 56x |

The GP results agree with this program's output digit for digit, which serves
as an additional independent cross-check (different zeta implementation,
different bignum stack above GMP).

## Reference ports in other languages

`ports/` holds single-file implementations of the same accelerated series for
other systems:

- `ports/khinchin.py` — fixed-point port on top of mpmath's internals.
  Measured 30x (1000 digits) to 37x (2000 digits) faster than mpmath's
  built-in `mp.khinchin` (the program on the OEIS entry), which uses the
  unaccelerated series; output verified against this repo's digit files.
- `ports/khinchin.gp` — the serial and parallel (`parvector`) PARI/GP
  scripts behind the comparison table above; both verified to 1000 digits.
- `ports/khinchin.mpl` — Maple version of the same series. No Maple
  installation on this machine, so reviewed but not machine-tested.
- `ports/khinchin.wl` — Mathematica. Packages the built-in `N[Khinchin, d]`
  exactly as the OEIS entry's program does: with no local Mathematica to
  benchmark against, there is no evidence a hand-rolled series would beat
  the built-in.

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
