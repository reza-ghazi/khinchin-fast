(* Khinchin's constant (OEIS A002210) in Mathematica.

   Mathematica ships an arbitrary-precision built-in for this constant,
   and the program on the OEIS entry simply uses it:

       RealDigits[N[Khinchin, 100]][[1]]
       (* Vladimir Joseph Stephan Orlovsky, Jun 18 2009 *)

   Mathematica is not installed on this machine, so the accelerated
   series from ../khinchin_fast.c could not be benchmarked against the
   built-in here; absent evidence it would be faster, this file just
   packages the OEIS/built-in approach. *)

(* d digits after the decimal point, as a decimal string. *)
KhinchinString[d_Integer?Positive] :=
  ToString @ NumberForm[N[Khinchin, d + 1], d + 1, ExponentFunction -> (Null &)]

(* Digit list, as on the OEIS entry. *)
KhinchinDigits[d_Integer?Positive] := First @ RealDigits @ N[Khinchin, d + 1]

(* Decimal value plus a newline into a text file, the same contract as
   ../khinchin_fast.c. *)
KhinchinToFile[d_Integer?Positive, path_String] :=
  Export[path, KhinchinString[d] <> "\n", "Text"]

(* Examples: KhinchinString[100]
             KhinchinToFile[1000, "khinchin-wl-1k.txt"] *)
