# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Single-file C program (`khinchin_fast.c`) that computes Khinchin's constant to a
requested number of decimal digits using FLINT/Arb ball arithmetic, OpenMP, and
GMP/MPFR. There is no test framework, no CI, and no VCS here — the project is
the one source file plus a Makefile and a reference digit file.

## Build and run

```sh
make                       # -O3 -march=native -flto -fopenmp, links flint/mpfr/gmp
make clean
./khinchin-fast 10000 out.txt          # digits AFTER the decimal point
./khinchin-fast --benchmark 10000      # compute only, no conversion/write
./khinchin-fast --force 10000 out.txt  # bypass the physical-memory preflight
KHINCHIN_BACKEND=arb ./khinchin-fast 10000 out-arb.txt   # cross-check backend
```

Requires `flint-devel`, `mpfr-devel`, `gmp-devel` (Fedora names). `-march=native`
means the binary is not portable off the build machine.

## Verification (the closest thing to a test suite)

`khinchin-100.txt` holds 100 correct digits after the decimal point (no trailing
newline). Regression check after any change to the math:

```sh
./khinchin-fast 100 /tmp/k100.txt
diff <(tr -d '\n' < /tmp/k100.txt) <(tr -d '\n' < khinchin-100.txt)
```

Then cross-check a larger size against the independent backend — the two
backends share only the outer BBC framing, so agreement is meaningful:

```sh
./khinchin-fast 10000 a.txt && KHINCHIN_BACKEND=arb ./khinchin-fast 10000 b.txt && diff a.txt b.txt
```

Timing/telemetry (`backend=`, `threads=`, `compute_seconds=`, …) goes to stderr
as `key=value` lines; the output file contains only the decimal value plus a
newline. Keep that split — stdout/file purity is what makes `diff` checks work.

## Architecture

Two interchangeable backends compute the same Bailey–Borwein–Crandall series for
`log2(K)`:

- `khinchin_parallel()` (`KHINCHIN_BACKEND=arb`) — calls `arb_zeta_ui()` per
  term at full precision, ascending `n`, equal-count thread blocks. Slow but
  independent; kept untouched as the cross-check oracle.
- `khinchin_scp()` (default) — the optimized backend described below.

### SCP backend structure

Key fact driving everything: term `n` of the accelerated sum has magnitude
about `N^(-2n)`, so it only needs about `P - 2n*log2(N)` accurate bits, and
`zeta(2n) - 1 - Σ_{k<N} k^(-2n) = Σ_{k>=N} k^(-2n)`.

The range `n = 1..M` splits at `n_direct ≈ P/(2(log2N+3))`:

- **Direct region** (`n >= n_direct`, `scp_direct_range()`): the term is
  computed as the literal tail sum over `k in [N, K]`, `K ≈ 2^(P/2n) <= ~8N`,
  all at low precision — no Bernoulli numbers, no cancellation, no
  full-precision arithmetic. A rigorous integral bound for the omitted `k > K`
  is added per term via `arb_add_error()`.
- **Bernoulli region** (`n < n_direct`, `scp_bernoulli_range()`): FLINT's
  `bernoulli_rev` iterator reconstructs `zeta(2n)`, which is near 1 and hence
  genuinely needs full precision there. Power-table entries and the tail
  assembly still use dropping per-entry/per-term precision.

Per-entry table precision is `P - 2*first*log2(k) + guard` where `first` is the
block's lowest `n` — relative errors persist through the incremental `*k^2`
updates, so an entry must start accurate enough for its final, most demanding
use. Blocks (4× threads) are sized for equal estimated cost (`term_cost()`),
sorted heaviest-first, and scheduled `omp for schedule(dynamic)`.

**Trap that already bit once:** each block must accumulate into its own fresh
`arb_t` (`block_sum`), then be added to the thread partial at full precision.
The range functions round their accumulator at the block's own (low) precision;
letting them add directly into a partial that already holds a larger earlier
block destroys it and certification fails.

Tuning/debug knobs (env vars): `KHINCHIN_N` overrides the table cut-off `N`,
`KHINCHIN_DIRECT` overrides `n_direct`, `KHINCHIN_VERBOSE=1` prints per-block
`n` ranges and timings. Timings at 50k digits were flat for `N` in ~67–100, so
the `precision^0.35` default stands.

The truncated `n > M` tail is bounded by `N^(-2M)` and folded in with
`arb_add_error()`, so the final ball is rigorous. Rounding is proof-carrying:
`rounded_scaled_integer()` multiplies by `10^digits`, adds 1/2, floors, and
requires `arb_get_unique_fmpz()` to succeed — never replace this with an
unchecked conversion. Guard bits come from `working_precision()`
(`log2(10)*(digits+1) + 128`); widen that rather than weakening the uniqueness
check if certification ever fails.

`estimated_peak_bytes()` predates the variable-precision tables and now
overestimates (conservative, safe direction) — keep it an upper bound if the
table sizes change.

## Conventions

C99-ish with `_POSIX_C_SOURCE 200809L`, 4-space indent, Allman braces, FLINT
types (`ulong`, `slong`, `arb_t`, `fmpz_t`) throughout. Every `arb_init`/
`fmpz_init` has a matching clear on all paths, including the error returns in
`main()`, which also call `flint_cleanup()`. Builds are warning-clean under
`-Wall -Wextra -Wpedantic`; keep them that way.

`NOTICE.md` records that the SCP approach is adapted from Remco Bloemen's
MIT-licensed `recmo/khinchin` and that FLINT is LGPL — update it if a new
third-party algorithm or dependency is introduced.

`Khinchin demo - 20260825-185607.txt` is an unrelated y-cruncher validation
artifact (a 100-digit finite-product demo), not an input or output of this
program.
