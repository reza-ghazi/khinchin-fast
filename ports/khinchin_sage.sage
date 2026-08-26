# Khinchin's constant (OEIS A002210) via the accelerated zeta series.
#
# Same mathematics as ../khinchin_fast.c:
#
#   ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
#                  + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
#
#   h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
#
# Sage ships FLINT/Arb, so this port computes in RealBallField — every
# zeta(2n) is arb_zeta at full precision (the strategy of the C
# program's KHINCHIN_BACKEND=arb backend) and the whole loop is rigorous
# ball arithmetic until the final midpoint extraction.  Serial.
#
# Tested with Sage 10.9: byte-identical output to the C program at 1000
# digits, digit-exact at 10,000 (4.2 s compute, serial; sage startup
# adds ~3.5 s of wall time). Testing required repairing the machine's
# editable Sage build first - its planarity module was grafted from
# upstream to survive Fedora 44's planarity 5.0 API break.
#
# Usage:  sage khinchin_sage.sage DIGITS [OUTPUT_FILE]   (default 100)
# With OUTPUT_FILE the decimal value plus a newline is written there
# (the same contract as ../khinchin_fast.c); otherwise stdout.

import sys


def khinchin_string(digits):
    prec = ceil((digits + 2) * log(10.0, 2)) + 32
    wp = prec + ceil(sqrt(prec)) + 64
    R = RealBallField(wp)
    N = max(3, floor(RDF(wp) ^ 0.35))
    M = ceil(wp * log(2.0) / (2 * log(RDF(N)))) + 1

    S = -sum((R(k - 1) / k).log() * (R(k + 1) / k).log()
             for k in range(2, N))

    ks = list(range(2, N))
    powers = [~(R(k) ^ 2) for k in ks]
    h = R(1)
    for n in range(1, M + 1):
        tail = R(2 * n).zeta() - 1 - sum(powers)
        S += tail / n * h
        h -= ~(R(2 * n) * (2 * n + 1))
        powers = [p / (k * k) for p, k in zip(powers, ks)]

    K = (S / R(2).log()).exp()
    scaled = (RealField(prec)(K.mid()) * 10 ^ digits).round()
    text = str(scaled)
    return text[0] + "." + text[1:]


args = sys.argv[1:]
digits = int(args[0]) if args else 100
result = khinchin_string(digits)
if len(args) > 1:
    with open(args[1], "w") as output:
        output.write(result + "\n")
else:
    print(result)
