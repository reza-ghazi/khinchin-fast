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
\\ each block seeding its alternating-harmonic weight h(first) by direct
\\ summation (much cheaper than psi at high precision).  These are
\\ the scripts behind the PARI/GP comparison table in ../README.md
\\ (21-56x slower than the C program, digit-for-digit identical).
\\
\\ Usage:  gp -q ports/khinchin.gp
\\         ? \p 1000
\\         ? khinchin(1000)
\\         ? khinchin_par(10000)
\\         ? khinchin_write(10000, "khinchin-gp-10k.txt")
\\ (set \p to the digits you want displayed; the functions compute at
\\ d + 20 internally via localprec.  khinchin_write stores the decimal
\\ value plus a newline, the same contract as ../khinchin_fast.c)

\\ parvector threads get their own PARI stacks; the default 8 MB cap
\\ overflows near 10k digits, so pre-size them and allow growth to 1 GB.
default(threadsize, max(default(threadsize), 1 << 26));
default(threadsizemax, max(default(threadsizemax), 1 << 30));

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
       h  = sum(j = 1, 2 * first - 1, (-1)^(j + 1) / (1. * j)),
       pw = vector(N - 2, i, (i + 1)^(-2. * first)),
       s  = 0.);
    for(n = first, last,
      s += (zeta(2*n) - 1 - vecsum(pw)) / n * h;
      h -= 1. / (2*n * (2*n + 1));
      pw = vector(N - 2, i, pw[i] / (i + 1)^2));
    s);
  exp((S + vecsum(parts)) / log(2));
}

khinchin_write(d, path) =
{
  my(K  = khinchin_par(d),
     sc = floor(K * 10^d + 1/2),
     fr = Strchr(Vecsmall(Str(10^d + sc % 10^d))[2..-1]),
     f  = fileopen(path, "w"));
  filewrite(f, Str(sc \ 10^d, ".", fr));
  fileclose(f);
}
