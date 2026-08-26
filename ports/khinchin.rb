# Khinchin's constant (OEIS A002210) via the accelerated zeta series.
#
# Same mathematics as ../khinchin_fast.c:
#
#   ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
#                  + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
#
#   h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
#
# Runs entirely on Ruby's native arbitrary-precision Integer in fixed
# point (scale 2^wp) — the "native bignum" data point for Ruby, no gems.
# Integers have no transcendentals, so the port carries its own
# fixed-point kit: pi by Machin's formula, logarithms via the atanh(1/q)
# series, exp by argument-halving plus Taylor. The even zeta values come
# from the positive-term recurrence (n + 1/2) zeta(2n) = sum_j zeta(2j)
# zeta(2n-2j) — but only below n_direct, the C program's two-region
# split: above it each term is a literal tail sum with no zeta values,
# shortening the O(M^2) recurrence to O(n_direct^2). Serial (the GVL
# rules out threads for CPU-bound work). Guard bits absorb the drift; nothing here
# is interval-certified.
#
# Usage:  ruby khinchin.rb DIGITS [OUTPUT_FILE]        (default 100)
# With OUTPUT_FILE the decimal value plus a newline is written there
# (the same contract as ../khinchin_fast.c); otherwise stdout.

# atanh(1/q) = sum_{j>=0} 1/((2j+1) q^(2j+1)); atan(1/q) alternates.
def arc_recip(q, one, alternating: false)
  q2 = q * q
  power = one / q
  s = power
  j = 1
  while power.positive?
    power /= q2
    term = power / (2 * j + 1)
    s += alternating && j.odd? ? -term : term
    j += 1
  end
  s
end

# exp(y) for 0 <= y < 1: halve 16 times, Taylor, square back.
def exp_fixed(y, wp, one)
  arg = y >> 16
  acc = one
  term = one
  j = 1
  while term.positive?
    term = ((term * arg) >> wp) / j
    acc += term
    j += 1
  end
  16.times { acc = (acc * acc) >> wp }
  acc
end

# K0 scaled by 2^wp, plus wp.
def khinchin_fixed(digits)
  prec = ((digits + 2) * Math.log2(10)).ceil + 32
  wp = prec + Math.sqrt(prec).ceil + 64
  one = 1 << wp
  n_cut = [3, (wp**0.35).floor].max
  m = (wp * Math.log(2) / (2 * Math.log(n_cut))).ceil + 1
  n_direct = [(wp / (2 * (Math.log2(n_cut) + 3))).ceil, m + 1].min

  # -ln((k-1)/k) ln((k+1)/k) = 4 atanh(1/(2k-1)) atanh(1/(2k+1)).
  s = 0
  (2...n_cut).each do |k|
    a = arc_recip(2 * k - 1, one)
    b = arc_recip(2 * k + 1, one)
    s += (4 * a * b) >> wp
  end

  # pi = 16 atan(1/5) - 4 atan(1/239); zeta(2) = pi^2/6 seeds the
  # positive-term recurrence (n + 1/2) z_n = sum z_j z_{n-j}.
  pi = 16 * arc_recip(5, one, alternating: true) -
       4 * arc_recip(239, one, alternating: true)
  z = Array.new(n_direct)
  z[1] = ((pi * pi) >> wp) / 6
  (2...n_direct).each do |n|
    half = (n - 1) / 2
    raw = 0
    (1..half).each { |j| raw += z[j] * z[n - j] }
    acc = 2 * (raw >> wp)
    acc += (z[n / 2] * z[n / 2]) >> wp if n.even?
    z[n] = (2 * acc) / (2 * n + 1)
  end

  # Main accelerated loop, ascending, incremental power table.
  ks = (2...n_cut).to_a
  powers = ks.map { |k| one / (k * k) }
  h = one
  (1...n_direct).each do |n|
    tail = z[n] - one
    powers.each { |p| tail -= p }
    s += (tail * h / n) >> wp
    h += one / (2 * n + 1) - one / (2 * n)
    ks.each_with_index { |k, i| powers[i] /= k * k }
  end

  # Direct region (the C program's split): literal tail sums over
  # k in [N, K], no zeta values; h continues ascending, power table
  # re-seeded per block of 32 terms.
  n = n_direct
  while n <= m
    last = [n + 31, m].min
    kk = [(2.0**((wp + 32.0) / (2 * n))).floor + 1, 16 * n_cut + 64].min
    kk = [kk, n_cut].max
    dks = (n_cut..kk).to_a
    pw = dks.map { |k| one / (k**(2 * n)) }
    (n..last).each do |mm|
      tail = 0
      pw.each { |p| tail += p }
      s += (tail * h / mm) >> wp
      h += one / (2 * mm + 1) - one / (2 * mm)
      dks.each_with_index { |k, i| pw[i] /= k * k }
    end
    n = last + 1
  end

  ln2 = 2 * arc_recip(3, one)
  [exp_fixed((s << wp) / ln2, wp, one), wp]
end

digits = ARGV[0] ? Integer(ARGV[0]) : 100
value, wp = khinchin_fixed(digits)
scaled = (value * 10**digits + (1 << (wp - 1))) >> wp
text = scaled.to_s
result = "#{text[0]}.#{text[1..]}"
if ARGV[1]
  File.write(ARGV[1], result + "\n")
else
  puts result
end
