(* Khinchin's constant (OEIS A002210) via the accelerated zeta series.

   Same mathematics as ../khinchin_fast.c:

     ln(2) ln(K0) = -Sum[Log[(k-1)/k] Log[(k+1)/k], {k, 2, N-1}]
                    + Sum[(Zeta[2n] - 1 - Sum[k^(-2n), {k, 2, N-1}])/n h[n], {n, 1, M}],

     h[n] = Sum[(-1)^(j+1)/j, {j, 1, 2n-1}].

   Mathematica evaluates Zeta[2n] symbolically through Bernoulli numbers,
   which makes the accelerated series far faster than the built-in
   Khinchin constant that the OEIS entry's program uses: measured with
   Mathematica 14.2 on this machine, 0.045 s vs 0.57 s at 1000 digits and
   the built-in needs 508 s at 10,000. The built-in is kept below as a
   cross-check.

   Usage:  Get["ports/khinchin.wl"];
           KhinchinString[1000]
           KhinchinToFile[1000, "khinchin-wl-1k.txt"]
   KhinchinToFile stores the decimal value plus a newline, the same
   contract as ../khinchin_fast.c. *)

KhinchinAccel[d_Integer?Positive] := Module[{bits, nn, m, wp, s, pw, h},
  bits = Ceiling[(d + 20) Log[2, 10]];
  wp = d + 25;
  nn = Max[3, Floor[bits^0.35]];
  m = Ceiling[bits Log[2]/(2 Log[nn])] + 1;
  s = -Sum[N[Log[(k - 1)/k] Log[(k + 1)/k], wp], {k, 2, nn - 1}];
  pw = Table[N[1/k^2, wp], {k, 2, nn - 1}];
  h = N[1, wp];
  Do[
    s += (N[Zeta[2 n], wp] - 1 - Total[pw]) h/n;
    h -= N[1/(2 n (2 n + 1)), wp];
    pw = pw/Table[k^2, {k, 2, nn - 1}],
    {n, 1, m}];
  Exp[s/Log[N[2, wp]]]];

(* d digits after the decimal point, as a decimal string. *)
KhinchinString[d_Integer?Positive] :=
  ToString @ NumberForm[N[KhinchinAccel[d], d + 1], d + 1,
    ExponentFunction -> (Null &)]

(* Digit list, as on the OEIS entry. *)
KhinchinDigits[d_Integer?Positive] :=
  First @ RealDigits @ N[KhinchinAccel[d], d + 1]

(* Decimal value plus a newline into a text file. *)
KhinchinToFile[d_Integer?Positive, path_String] :=
  Export[path, KhinchinString[d] <> "\n", "Text"]

(* The built-in constant, as the OEIS entry's program uses it
   (RealDigits[N[Khinchin, 100]][[1]], V. J. S. Orlovsky, 2009).
   Much slower; kept as an independent cross-check. *)
KhinchinBuiltin[d_Integer?Positive] := N[Khinchin, d + 1]
