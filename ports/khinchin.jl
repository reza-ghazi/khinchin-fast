# Khinchin's constant (OEIS A002210) via the accelerated zeta series.
#
# Same mathematics as ../khinchin_fast.c:
#
#   ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
#                  + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
#
#   h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
#
# Runs on BigFloat (MPFR) with no packages.  The even zeta values come
# from the classical positive-term recurrence
#
#   (n + 1/2) zeta(2n) = sum_{j=1}^{n-1} zeta(2j) zeta(2n-2j),   zeta(2) = pi^2/6,
#
# which is numerically benign (every term is positive, so accumulated
# relative error grows only about linearly in n) and keeps the port free
# of Bernoulli numbers and external zeta implementations.  The guard
# bits absorb that drift; unlike the C program, nothing here is
# interval-certified.
#
# Usage:  julia khinchin.jl [digits]        (default 100)

function khinchin(digits::Integer)
    prec = ceil(Int, (digits + 2) * log2(10)) + 32
    wp = prec + ceil(Int, sqrt(prec)) + 64
    setprecision(BigFloat, wp) do
        N = max(3, floor(Int, wp^0.35))
        M = ceil(Int, wp * log(2) / (2 * log(N))) + 1

        S = -sum(log(big(k - 1) / k) * log(big(k + 1) / k) for k in 2:N-1)

        # zeta(2), zeta(4), ..., zeta(2M) by the positive recurrence.
        z = Vector{BigFloat}(undef, M)
        z[1] = big(pi)^2 / 6
        for n in 2:M
            acc = zero(BigFloat)
            for j in 1:n-1
                acc += z[j] * z[n-j]
            end
            z[n] = 2acc / (2n + 1)
        end

        powers = [inv(big(k)^2) for k in 2:N-1]
        h = one(BigFloat)
        for n in 1:M
            tail = z[n] - 1 - sum(powers)
            S += tail * h / n
            h -= inv(big(2n) * (2n + 1))
            for (i, k) in enumerate(2:N-1)
                powers[i] /= k^2
            end
        end
        exp(S / log(big(2)))
    end
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
    println(khinchin_string(isempty(ARGS) ? 100 : parse(Int, ARGS[1])))
end
