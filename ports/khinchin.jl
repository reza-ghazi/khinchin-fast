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
# Multithreaded: run with `julia -t auto khinchin.jl [digits]`.  The
# recurrence convolutions (the dominant cost, using their j <-> n-j
# symmetry) fan out across threads, and the main loop splits into
# per-thread blocks seeded like the C program's, with h(first) summed
# directly and the power table started at k^(-2 first).
#
# Usage:  julia -t auto khinchin.jl DIGITS [OUTPUT_FILE]   (default 100)
# With OUTPUT_FILE the decimal value plus a newline is written there
# (the same contract as ../khinchin_fast.c); otherwise stdout.

using Base.Threads

function khinchin(digits::Integer)
    prec = ceil(Int, (digits + 2) * log2(10)) + 32
    wp = prec + ceil(Int, sqrt(prec)) + 64
    setprecision(BigFloat, wp) do
        N = max(3, floor(Int, wp^0.35))
        M = ceil(Int, wp * log(2) / (2 * log(N))) + 1

        S = -sum(log(big(k - 1) / k) * log(big(k + 1) / k) for k in 2:N-1)

        # zeta(2), zeta(4), ..., zeta(2M) by the positive recurrence,
        # using the symmetry of the convolution: only j <= (n-1)/2 is
        # summed, doubled, plus the middle square when n is even.
        z = Vector{BigFloat}(undef, M)
        z[1] = big(pi)^2 / 6
        for n in 2:M
            half = (n - 1) ÷ 2
            acc = if half >= 64 && nthreads() > 1
                nch = nthreads()
                parts = Vector{BigFloat}(undef, nch)
                @threads :static for c in 1:nch
                    lo = 1 + ((c - 1) * half) ÷ nch
                    hi = (c * half) ÷ nch
                    a = zero(BigFloat)
                    for j in lo:hi
                        a += z[j] * z[n-j]
                    end
                    parts[c] = a
                end
                sum(parts)
            else
                a = zero(BigFloat)
                for j in 1:half
                    a += z[j] * z[n-j]
                end
                a
            end
            acc *= 2
            iseven(n) && (acc += z[n÷2]^2)
            z[n] = 2acc / (2n + 1)
        end

        # Main loop in equal-count blocks, one task per block.
        W = max(1, min(nthreads(), M ÷ 8))
        partials = Vector{BigFloat}(undef, W)
        @threads :static for b in 1:W
            first = 1 + (M * (b - 1)) ÷ W
            last = (M * b) ÷ W
            h = zero(BigFloat)
            for j in 1:2first-1
                h += (isodd(j) ? 1 : -1) / big(j)
            end
            powers = [inv(big(k)^(2first)) for k in 2:N-1]
            s = zero(BigFloat)
            for n in first:last
                tail = z[n] - 1 - sum(powers)
                s += tail * h / n
                h -= inv(big(2n) * (2n + 1))
                for (i, k) in enumerate(2:N-1)
                    powers[i] /= k^2
                end
            end
            partials[b] = s
        end
        S += sum(partials)

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
    text = khinchin_string(isempty(ARGS) ? 100 : parse(Int, ARGS[1]))
    if length(ARGS) >= 2
        open(f -> println(f, text), ARGS[2], "w")
    else
        println(text)
    end
end
