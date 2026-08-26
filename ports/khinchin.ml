(* Khinchin's constant (OEIS A002210) via the accelerated zeta series.

   Same mathematics as ../khinchin_fast.c:

     ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
                    + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),

     h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.

   Runs on Zarith (OCaml's GMP-backed integers) in fixed point (scale
   2^wp).  The port carries its own fixed-point kit: pi by Machin's
   formula, logarithms via the atanh(1/q) series, exp by
   argument-halving plus Taylor.  The even zeta values come from the
   positive-term recurrence (n + 1/2) zeta(2n) = sum_j zeta(2j)
   zeta(2n-2j) — but only below n_direct, the C program's two-region
   split: above it each term is a literal tail sum with no zeta values,
   shortening the O(M^2) recurrence to O(n_direct^2).  Guard bits
   absorb the drift; nothing here is interval-certified.

   Tested with OCaml 5.2.1 (opam) and Zarith 1.14: byte-identical
   output to the C program at 1000 digits (0.07 s), digit-exact at
   10,000 (56 s serial).

   Build and run:
     ocamlfind ocamlopt -package zarith -linkpkg khinchin.ml -o khinchin-ml
     ./khinchin-ml DIGITS [OUTPUT_FILE]

   With OUTPUT_FILE the decimal value plus a newline is written there
   (the same contract as ../khinchin_fast.c); otherwise stdout. *)

let wp = ref 0
let one = ref Z.zero

let mul_shift a b = Z.shift_right (Z.mul a b) !wp

(* atanh(1/q) = sum_{j>=0} 1/((2j+1) q^(2j+1)); atan(1/q) alternates. *)
let arc_recip ?(alternating = false) q =
  let q2 = Z.of_int (q * q) in
  let power = ref (Z.div !one (Z.of_int q)) in
  let s = ref !power in
  let j = ref 1 in
  while Z.sign !power > 0 do
    power := Z.div !power q2;
    let term = Z.div !power (Z.of_int ((2 * !j) + 1)) in
    if alternating && !j mod 2 = 1 then s := Z.sub !s term
    else s := Z.add !s term;
    incr j
  done;
  !s

(* exp(y) for 0 <= y < 1: halve 16 times, Taylor, square back. *)
let exp_fixed y =
  let arg = Z.shift_right y 16 in
  let acc = ref !one in
  let term = ref !one in
  let j = ref 1 in
  while Z.sign !term > 0 do
    term := Z.div (mul_shift !term arg) (Z.of_int !j);
    acc := Z.add !acc !term;
    incr j
  done;
  for _ = 1 to 16 do
    acc := mul_shift !acc !acc
  done;
  !acc

(* K0 scaled by 2^wp. *)
let khinchin_fixed digits =
  let prec =
    int_of_float (ceil (float_of_int (digits + 2) *. log 10. /. log 2.)) + 32
  in
  wp := prec + int_of_float (ceil (sqrt (float_of_int prec))) + 64;
  one := Z.shift_left Z.one !wp;
  let big_n = max 3 (int_of_float (float_of_int !wp ** 0.35)) in
  let m =
    int_of_float
      (ceil (float_of_int !wp *. log 2. /. (2. *. log (float_of_int big_n))))
    + 1
  in
  let n_direct =
    min
      (int_of_float
         (ceil
            (float_of_int !wp
            /. (2. *. ((log (float_of_int big_n) /. log 2.) +. 3.)))))
      (m + 1)
  in

  (* -ln((k-1)/k) ln((k+1)/k) = 4 atanh(1/(2k-1)) atanh(1/(2k+1)). *)
  let s = ref Z.zero in
  for k = 2 to big_n - 1 do
    let a = arc_recip ((2 * k) - 1) and b = arc_recip ((2 * k) + 1) in
    s := Z.add !s (Z.mul (Z.of_int 4) (mul_shift a b))
  done;

  (* pi = 16 atan(1/5) - 4 atan(1/239); zeta(2) = pi^2/6 seeds the
     positive-term recurrence (n + 1/2) z_n = sum z_j z_{n-j}. *)
  let pi =
    Z.sub
      (Z.mul (Z.of_int 16) (arc_recip ~alternating:true 5))
      (Z.mul (Z.of_int 4) (arc_recip ~alternating:true 239))
  in
  let z = Array.make (m + 1) Z.zero in
  z.(1) <- Z.div (mul_shift pi pi) (Z.of_int 6);
  for n = 2 to n_direct - 1 do
    let half = (n - 1) / 2 in
    let raw = ref Z.zero in
    for j = 1 to half do
      raw := Z.add !raw (Z.mul z.(j) z.(n - j))
    done;
    let acc = ref (Z.mul (Z.of_int 2) (Z.shift_right !raw !wp)) in
    if n mod 2 = 0 then acc := Z.add !acc (mul_shift z.(n / 2) z.(n / 2));
    z.(n) <- Z.div (Z.mul (Z.of_int 2) !acc) (Z.of_int ((2 * n) + 1))
  done;

  (* Main accelerated loop, ascending, incremental power table. *)
  let ks = Array.init (big_n - 2) (fun i -> i + 2) in
  let powers =
    Array.map (fun k -> Z.div !one (Z.of_int (k * k))) ks
  in
  let h = ref !one in
  for n = 1 to n_direct - 1 do
    let tail = ref (Z.sub z.(n) !one) in
    Array.iter (fun p -> tail := Z.sub !tail p) powers;
    s :=
      Z.add !s
        (Z.shift_right (Z.div (Z.mul !tail !h) (Z.of_int n)) !wp);
    h :=
      Z.sub
        (Z.add !h (Z.div !one (Z.of_int ((2 * n) + 1))))
        (Z.div !one (Z.of_int (2 * n)));
    Array.iteri
      (fun i k -> powers.(i) <- Z.div powers.(i) (Z.of_int (k * k)))
      ks
  done;

  (* Direct region (the C program's split): literal tail sums over
     k in [N, K], no zeta values; h continues ascending, power table
     re-seeded per block of 32 terms. *)
  let n = ref n_direct in
  while !n <= m do
    let last = min (!n + 31) m in
    let kk =
      max big_n
        (min
           (int_of_float
              (2. ** ((float_of_int !wp +. 32.) /. (2. *. float_of_int !n)))
           + 1)
           ((16 * big_n) + 64))
    in
    let dks = Array.init (kk - big_n + 1) (fun i -> i + big_n) in
    let pw =
      Array.map
        (fun k -> Z.div !one (Z.pow (Z.of_int k) (2 * !n)))
        dks
    in
    for mm = !n to last do
      let tail = ref Z.zero in
      Array.iter (fun p -> tail := Z.add !tail p) pw;
      s :=
        Z.add !s
          (Z.shift_right (Z.div (Z.mul !tail !h) (Z.of_int mm)) !wp);
      h :=
        Z.sub
          (Z.add !h (Z.div !one (Z.of_int ((2 * mm) + 1))))
          (Z.div !one (Z.of_int (2 * mm)));
      Array.iteri
        (fun i k -> pw.(i) <- Z.div pw.(i) (Z.of_int (k * k)))
        dks
    done;
    n := last + 1
  done;

  let ln2 = Z.mul (Z.of_int 2) (arc_recip 3) in
  exp_fixed (Z.div (Z.shift_left !s !wp) ln2)

let () =
  let digits =
    if Array.length Sys.argv > 1 then int_of_string Sys.argv.(1) else 100
  in
  let value = khinchin_fixed digits in
  let scaled =
    Z.shift_right
      (Z.add
         (Z.mul value (Z.pow (Z.of_int 10) digits))
         (Z.shift_left Z.one (!wp - 1)))
      !wp
  in
  let text = Z.to_string scaled in
  let result = String.sub text 0 1 ^ "." ^ String.sub text 1 (String.length text - 1) in
  if Array.length Sys.argv > 2 then (
    let oc = open_out Sys.argv.(2) in
    output_string oc (result ^ "\n");
    close_out oc)
  else print_endline result
