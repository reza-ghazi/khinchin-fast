# Third-party notices

The SCP-style zeta acceleration was adapted from the algorithmic approach in
Remco Bloemen's `khinchin` project:

<https://github.com/recmo/khinchin>

That project is Copyright (c) 2022 Remco Bloemen and licensed under the MIT
License.

This program links to FLINT, Arb, MPFR, GMP, and the OpenMP runtime. The reverse
Bernoulli iterator and rigorous ball arithmetic are provided by FLINT:

<https://github.com/flintlib/flint>

FLINT is distributed under the GNU Lesser General Public License, version 3 or
later. MPFR and GMP are distributed under their respective GNU licenses.

The reference ports under `ports/` use further third-party software as
external libraries or interpreters (nothing from them is redistributed
here): `khinchin.py` builds on mpmath (BSD 3-clause); `khinchin.gp` runs
under PARI/GP (GPL v2+); `khinchin-rs` depends on the Rust crates rug
(LGPL v3+) and rayon (MIT/Apache-2.0); `khinchin.f90` and `khinchin.jl`
bind or wrap GNU MPFR (LGPL v3+).
