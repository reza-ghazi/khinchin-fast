# Khinchin's constant (OEIS A002210) via the accelerated zeta series.
#
# Same mathematics as ../khinchin_fast.c:
#
#   ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
#                  + sum_{n>=1} (Zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
#
#   h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
#
# Maple evaluates Zeta at even integer arguments through Bernoulli
# numbers, so pulling the small k out of the zeta tails cuts the term
# count from ~bits/2 to ~bits/(2 log2 N), the same acceleration the C
# program uses.  There is no Maple program on the OEIS entry to compare
# against.
#
# NOTE: written without access to a Maple installation - reviewed
# against the tested Python and PARI/GP ports in this directory, but
# not machine-tested.
#
# Usage:  read "khinchin.mpl": khinchin(1000);

khinchin := proc(d::posint)
    local bits, N, M, S, pw, h, k, n, oldDigits, result;
    oldDigits := Digits;
    Digits := d + 20;
    bits := ceil((d + 20) * ln(10.0) / ln(2.0));
    N := max(3, floor(evalf(bits^0.35)));
    M := ceil(bits * ln(2.0) / (2 * ln(N))) + 1;
    S := -add(evalf(ln((k-1)/k) * ln((k+1)/k)), k = 2 .. N-1);
    pw := Array(2 .. N-1, [seq(evalf(1/k^2), k = 2 .. N-1)]);
    h := 1.0;
    for n from 1 to M do
        S := S + (evalf(Zeta(2*n)) - 1 - add(pw[k], k = 2 .. N-1)) / n * h;
        h := h - evalf(1/(2*n*(2*n+1)));
        for k from 2 to N-1 do pw[k] := pw[k]/k^2 end do;
    end do;
    result := exp(S / ln(2.0));
    Digits := oldDigits;
    result;
end proc:
