# Khinchin's constant (OEIS A002210) via the accelerated zeta series.
#
# Same mathematics as ../khinchin_fast.c:
#
#   ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
#                  + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
#
#   h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
#
# The heavy lifting is done by FLINT/Arb — the same library behind the C
# program — reached with plain `ccall` against the system libflint, using
# only Arb's pointer-based API (_arb_vec_init) so no C struct layouts are
# needed and no Julia packages are required.  Each zeta(2n) comes from
# arb_zeta_ui at full precision (the strategy of the C program's
# KHINCHIN_BACKEND=arb cross-check backend); the n-range splits across
# Julia threads in equal-count blocks, and the final ball's midpoint is
# pulled into a BigFloat through arb_get_interval_mpfr.
#
# Requires libflint (Fedora: flint or flint-devel).  Multithreaded: run
# with `julia -t auto`.
#
# Usage:  julia -t auto khinchin.jl DIGITS [OUTPUT_FILE]   (default 100)
# With OUTPUT_FILE the decimal value plus a newline is written there
# (the same contract as ../khinchin_fast.c); otherwise stdout.

using Base.Threads

const libflint = "libflint"

arb_new() = ccall((:_arb_vec_init, libflint), Ptr{Cvoid}, (Clong,), 1)
arb_free(x) = ccall((:_arb_vec_clear, libflint), Cvoid, (Ptr{Cvoid}, Clong), x, 1)
arb_zero!(x) = ccall((:arb_zero, libflint), Cvoid, (Ptr{Cvoid},), x)
arb_one!(x) = ccall((:arb_one, libflint), Cvoid, (Ptr{Cvoid},), x)
arb_set_ui!(x, u) = ccall((:arb_set_ui, libflint), Cvoid, (Ptr{Cvoid}, Culong), x, u)
for (jl, c) in ((:arb_add!, :arb_add), (:arb_sub!, :arb_sub),
                (:arb_mul!, :arb_mul), (:arb_div!, :arb_div))
    @eval $jl(r, a, b, prec) =
        ccall(($(QuoteNode(c)), libflint), Cvoid,
              (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Clong), r, a, b, prec)
end
for (jl, c) in ((:arb_sub_ui!, :arb_sub_ui), (:arb_div_ui!, :arb_div_ui),
                (:arb_mul_ui!, :arb_mul_ui), (:arb_pow_ui!, :arb_pow_ui))
    @eval $jl(r, a, u, prec) =
        ccall(($(QuoteNode(c)), libflint), Cvoid,
              (Ptr{Cvoid}, Ptr{Cvoid}, Culong, Clong), r, a, u, prec)
end
arb_inv!(r, a, prec) = ccall((:arb_inv, libflint), Cvoid,
    (Ptr{Cvoid}, Ptr{Cvoid}, Clong), r, a, prec)
arb_log!(r, a, prec) = ccall((:arb_log, libflint), Cvoid,
    (Ptr{Cvoid}, Ptr{Cvoid}, Clong), r, a, prec)
arb_exp!(r, a, prec) = ccall((:arb_exp, libflint), Cvoid,
    (Ptr{Cvoid}, Ptr{Cvoid}, Clong), r, a, prec)
arb_log_ui!(r, u, prec) = ccall((:arb_log_ui, libflint), Cvoid,
    (Ptr{Cvoid}, Culong, Clong), r, u, prec)
arb_zeta_ui!(r, u, prec) = ccall((:arb_zeta_ui, libflint), Cvoid,
    (Ptr{Cvoid}, Culong, Clong), r, u, prec)
arb_get_interval_mpfr!(lo, hi, x) =
    ccall((:arb_get_interval_mpfr, libflint), Cvoid,
          (Ref{BigFloat}, Ref{BigFloat}, Ptr{Cvoid}), lo, hi, x)

# Accelerated terms for n in [first, last]: zeta(2n) via arb_zeta_ui,
# shared incremental power table, ascending harmonic weight.  Writes the
# block sum into sb.  Pure ccall work — safe on any Julia thread.
function zeta_block!(sb, first, last, N, wp)
    zt, h, term, recip = arb_new(), arb_new(), arb_new(), arb_new()
    arb_zero!(sb)

    # h(first) = 1 - sum_{n<first} 1/(2n(2n+1)).
    arb_one!(h)
    for n in 1:first-1
        arb_set_ui!(recip, 2n)
        arb_mul_ui!(recip, recip, 2n + 1, wp)
        arb_inv!(recip, recip, wp)
        arb_sub!(h, h, recip, wp)
    end

    # Power table entry k starts at k^(-2(first-1)); advanced once per n.
    pw = [arb_new() for _ in 2:N-1]
    for (i, k) in enumerate(2:N-1)
        arb_set_ui!(pw[i], k)
        arb_pow_ui!(pw[i], pw[i], 2 * (first - 1), wp)
        arb_inv!(pw[i], pw[i], wp)
    end

    for n in first:last
        arb_zeta_ui!(zt, 2n, wp)
        arb_sub_ui!(zt, zt, 1, wp)
        for (i, k) in enumerate(2:N-1)
            arb_div_ui!(pw[i], pw[i], k * k, wp)
            arb_sub!(zt, zt, pw[i], wp)
        end
        arb_div_ui!(term, zt, n, wp)
        arb_mul!(term, term, h, wp)
        arb_add!(sb, sb, term, wp)
        arb_set_ui!(recip, 2n)
        arb_mul_ui!(recip, recip, 2n + 1, wp)
        arb_inv!(recip, recip, wp)
        arb_sub!(h, h, recip, wp)
    end

    foreach(arb_free, pw)
    arb_free(recip); arb_free(term); arb_free(h); arb_free(zt)
    return nothing
end

function khinchin(digits::Integer)
    prec = ceil(Int, (digits + 2) * log2(10)) + 32
    wp = prec + ceil(Int, sqrt(prec)) + 64
    N = max(3, floor(Int, wp^0.35))
    M = ceil(Int, wp * log(2) / (2 * log(N))) + 1

    s, t1, t2 = arb_new(), arb_new(), arb_new()
    arb_zero!(s)

    # Finite logarithmic correction.
    for k in 2:N-1
        arb_set_ui!(t1, k - 1)
        arb_div_ui!(t1, t1, k, wp)
        arb_log!(t1, t1, wp)
        arb_set_ui!(t2, k + 1)
        arb_div_ui!(t2, t2, k, wp)
        arb_log!(t2, t2, wp)
        arb_mul!(t1, t1, t2, wp)
        arb_sub!(s, s, t1, wp)
    end

    # Zeta range in equal-count blocks, one task per thread.
    W = max(1, min(nthreads(), M))
    partials = [arb_new() for _ in 1:W]
    @threads :static for b in 1:W
        zeta_block!(partials[b], 1 + (M * (b - 1)) ÷ W, (M * b) ÷ W, N, wp)
    end
    for p in partials
        arb_add!(s, s, p, wp)
        arb_free(p)
    end

    arb_log_ui!(t1, 2, wp)
    arb_div!(s, s, t1, wp)
    arb_exp!(s, s, wp)

    # Midpoint out through MPFR into a BigFloat.
    result = setprecision(BigFloat, wp) do
        lo = BigFloat()
        hi = BigFloat()
        arb_get_interval_mpfr!(lo, hi, s)
        lo
    end
    arb_free(t2); arb_free(t1); arb_free(s)
    result
end

"""Decimal string of K0 with `digits` digits after the point."""
function khinchin_string(digits::Integer)
    prec = ceil(Int, (digits + 2) * log2(10)) + 32
    scaled = setprecision(() -> round(BigInt, khinchin(digits) * big(10)^digits),
                          BigFloat, prec)
    text = string(scaled)
    text[1:1] * "." * text[2:end]
end

if abspath(PROGRAM_FILE) == @__FILE__
    text = khinchin_string(isempty(ARGS) ? 100 : parse(Int, ARGS[1]))
    if length(ARGS) >= 2
        open(f -> println(f, text), ARGS[2], "w")
    else
        println(text)
    end
end
