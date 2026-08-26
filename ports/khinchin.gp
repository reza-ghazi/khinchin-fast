\\ Khinchin's constant (OEIS A002210) via the accelerated zeta series.
\\
\\ Same mathematics as ../khinchin_fast.c:
\\
\\   ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
\\                  + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
\\
\\   h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j = psi(2n) - psi(n).
\\
\\ PARI has no built-in for this constant (the OEIS entry's PARI program
\\ was removed in 2010).  khinchin(d) is the serial version;
\\ khinchin_par(d) splits the zeta range across threads with parvector,
\\ seeding each block's alternating-harmonic weight with psi.  These are
\\ the scripts behind the PARI/GP comparison table in ../README.md
\\ (21-56x slower than the C program, digit-for-digit identical).
\\
\\ Usage:  gp -q ports/khinchin.gp
\\         ? \p 1000
\\         ? khinchin(1000)
\\         ? khinchin_par(10000)
\\ (set \p to the digits you want displayed; the functions compute at
\\ d + 20 internally via localprec)

khinchin(d) =
{
  my(bits, N, M, S, pw, h);
  localprec(d + 20);
  bits = ceil((d + 20) * log(10) / log(2));
  N = max(3, floor(bits^0.35));
  M = ceil(bits * log(2) / (2 * log(N))) + 1;
  S = -sum(k = 2, N - 1, log((k-1)/k) * log((k+1)/k) * 1.);
  pw = vector(N - 2, i, 1. / (i + 1)^2);   \\ pw[i] = k^(-2n), k = i+1
  h = 1.;
  for(n = 1, M,
    S += (zeta(2*n) - 1 - vecsum(pw)) / n * h;
    h -= 1. / (2*n * (2*n + 1));
    pw = vector(N - 2, i, pw[i] / (i + 1)^2));
  exp(S / log(2));
}

khinchin_par(d, blocks = 4 * default(nbthreads)) =
{
  my(bits, N, M, S, parts);
  localprec(d + 20);
  bits = ceil((d + 20) * log(10) / log(2));
  N = max(3, floor(bits^0.35));
  M = ceil(bits * log(2) / (2 * log(N))) + 1;
  blocks = min(blocks, M);
  S = -sum(k = 2, N - 1, log((k-1)/k) * log((k+1)/k) * 1.);
  parts = parvector(blocks, b,
    my(first = 1 + (M * (b - 1)) \ blocks,
       last  = (M * b) \ blocks,
       h  = psi(2. * first) - psi(1. * first),
       pw = vector(N - 2, i, (i + 1)^(-2. * first)),
       s  = 0.);
    for(n = first, last,
      s += (zeta(2*n) - 1 - vecsum(pw)) / n * h;
      h -= 1. / (2*n * (2*n + 1));
      pw = vector(N - 2, i, pw[i] / (i + 1)^2));
    s);
  exp((S + vecsum(parts)) / log(2));
}
