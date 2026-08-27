// Khinchin's constant (OEIS A002210) via the accelerated zeta series.
//
// Same mathematics as ../../khinchin_fast.c:
//
//   ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
//                  + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
//
//   h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
//
// Runs entirely on System.Numerics.BigInteger in fixed point (scale
// 2^wp) — the "native bignum" data point for C#, no packages. Like its
// siblings the port carries its own fixed-point kit (pi by Machin's
// formula, logarithms via the atanh(1/q) series, exp by
// argument-halving plus Taylor) and takes the even zeta values from the
// positive-term recurrence (n + 1/2) zeta(2n) = sum_j zeta(2j)
// zeta(2n-2j) — but only below n_direct, the C program's two-region
// split: above it each term is a literal tail sum with no zeta values.
// The recurrence convolutions (halved via their j <-> n-j symmetry) fan
// out over Parallel.For. Guard bits absorb the drift; nothing here is
// interval-certified.
//
// Run:  dotnet run -c Release --project khinchin-cs -- DIGITS [OUTPUT_FILE]
//
// With OUTPUT_FILE the decimal value plus a newline is written there
// (the same contract as ../../khinchin_fast.c); otherwise stdout.

using System.Numerics;

static class Khinchin
{
    static int wp;
    static BigInteger one;

    // atanh(1/q) = sum_{j>=0} 1/((2j+1) q^(2j+1)); atan(1/q) alternates.
    static BigInteger ArcRecip(long q, bool alternating)
    {
        BigInteger q2 = q * q;
        BigInteger power = one / q;
        BigInteger s = power;
        for (long j = 1; power > 0; j++)
        {
            power /= q2;
            BigInteger term = power / (2 * j + 1);
            s += alternating && j % 2 == 1 ? -term : term;
        }
        return s;
    }

    // exp(y) for 0 <= y < 1: halve 16 times, Taylor, square back.
    static BigInteger ExpFixed(BigInteger y)
    {
        BigInteger arg = y >> 16;
        BigInteger acc = one, term = one;
        for (long j = 1; term > 0; j++)
        {
            term = ((term * arg) >> wp) / j;
            acc += term;
        }
        for (int i = 0; i < 16; i++)
            acc = (acc * acc) >> wp;
        return acc;
    }

    // K0 scaled by 2^wp.
    static BigInteger KhinchinFixed(int digits)
    {
        int prec = (int)Math.Ceiling((digits + 2) * Math.Log2(10)) + 32;
        wp = prec + (int)Math.Ceiling(Math.Sqrt(prec)) + 64;
        one = BigInteger.One << wp;
        int bigN = Math.Max(3, (int)Math.Pow(wp, 0.35));
        int m = (int)Math.Ceiling(wp * Math.Log(2) / (2 * Math.Log(bigN))) + 1;
        int nDirect = Math.Min(
            (int)Math.Ceiling(wp / (2 * (Math.Log2(bigN) + 3))), m + 1);

        // -ln((k-1)/k) ln((k+1)/k) = 4 atanh(1/(2k-1)) atanh(1/(2k+1)).
        BigInteger s = 0;
        for (long k = 2; k < bigN; k++)
            s += (4 * (ArcRecip(2 * k - 1, false) * ArcRecip(2 * k + 1, false))) >> wp;

        // pi = 16 atan(1/5) - 4 atan(1/239); zeta(2) = pi^2/6 seeds the
        // positive-term recurrence, convolutions via Parallel.For.
        BigInteger pi = 16 * ArcRecip(5, true) - 4 * ArcRecip(239, true);
        var z = new BigInteger[nDirect];
        z[1] = ((pi * pi) >> wp) / 6;
        for (int n = 2; n < nDirect; n++)
        {
            int half = (n - 1) / 2;
            BigInteger raw = 0;
            if (half >= 64)
            {
                int nn = n;
                object gate = new();
                Parallel.For(1, half + 1,
                    () => BigInteger.Zero,
                    (j, _, local) => local + z[j] * z[nn - j],
                    local => { lock (gate) raw += local; });
            }
            else
            {
                for (int j = 1; j <= half; j++)
                    raw += z[j] * z[n - j];
            }
            BigInteger acc = (raw >> wp) << 1;
            if (n % 2 == 0)
                acc += (z[n / 2] * z[n / 2]) >> wp;
            z[n] = (acc << 1) / (2 * n + 1);
        }

        // Main accelerated loop, ascending, incremental power table.
        var powers = new BigInteger[bigN - 2];
        for (int k = 2; k < bigN; k++)
            powers[k - 2] = one / (k * k);
        BigInteger h = one;
        for (long n = 1; n < nDirect; n++)
        {
            BigInteger tail = z[n] - one;
            foreach (var pk in powers)
                tail -= pk;
            s += (tail * h / n) >> wp;
            h += one / (2 * n + 1) - one / (2 * n);
            for (int k = 2; k < bigN; k++)
                powers[k - 2] /= (long)k * k;
        }

        // Direct region (the C program's split): literal tail sums over
        // k in [N, K], no zeta values; h continues ascending, power
        // table re-seeded per block of 32 terms.
        for (int n = nDirect; n <= m;)
        {
            int last = Math.Min(n + 31, m);
            int bigK = Math.Min(
                (int)Math.Pow(2, (wp + 32.0) / (2.0 * n)) + 1, 16 * bigN + 64);
            bigK = Math.Max(bigK, bigN);
            var pw = new BigInteger[bigK - bigN + 1];
            for (int k = bigN; k <= bigK; k++)
                pw[k - bigN] = one / BigInteger.Pow(k, 2 * n);
            for (long mm = n; mm <= last; mm++)
            {
                BigInteger tail = 0;
                foreach (var pk in pw)
                    tail += pk;
                s += (tail * h / mm) >> wp;
                h += one / (2 * mm + 1) - one / (2 * mm);
                for (int k = bigN; k <= bigK; k++)
                    pw[k - bigN] /= (long)k * k;
            }
            n = last + 1;
        }

        BigInteger ln2 = ArcRecip(3, false) << 1;
        return ExpFixed((s << wp) / ln2);
    }

    static void Main(string[] args)
    {
        int digits = args.Length > 0 ? int.Parse(args[0]) : 100;
        BigInteger value = KhinchinFixed(digits);
        BigInteger scaled =
            (value * BigInteger.Pow(10, digits) + (BigInteger.One << (wp - 1))) >> wp;
        string text = scaled.ToString();
        string result = text[..1] + "." + text[1..];
        if (args.Length > 1)
            File.WriteAllText(args[1], result + "\n");
        else
            Console.WriteLine(result);
    }
}
