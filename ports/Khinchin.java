// Khinchin's constant (OEIS A002210) via the accelerated zeta series.
//
// Same mathematics as ../khinchin_fast.c:
//
//   ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
//                  + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
//
//   h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
//
// Runs entirely on java.math.BigInteger in fixed point (scale 2^wp) —
// the "native bignum" data point for Java, no dependencies. BigInteger
// has no transcendentals, so the port carries its own fixed-point kit:
// pi by Machin's formula, logarithms via the atanh(1/q) series, exp by
// argument-halving plus Taylor. The even zeta values come from the
// positive-term recurrence (n + 1/2) zeta(2n) = sum_j zeta(2j)
// zeta(2n-2j) — O(M^2) bignum products, with the convolutions (halved
// via their j <-> n-j symmetry) run through parallel streams. Guard
// bits absorb the drift; nothing here is interval-certified.
//
// Run:  java Khinchin.java DIGITS [OUTPUT_FILE]        (default 100)
//
// With OUTPUT_FILE the decimal value plus a newline is written there
// (the same contract as ../khinchin_fast.c); otherwise stdout.

import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.stream.IntStream;

public final class Khinchin {
    static int wp;
    static BigInteger one;

    static BigInteger mulShift(BigInteger a, BigInteger b) {
        return a.multiply(b).shiftRight(wp);
    }

    // atanh(1/q) = sum_{j>=0} 1/((2j+1) q^(2j+1)); atan(1/q) alternates.
    static BigInteger arcRecip(long q, boolean alternating) {
        BigInteger q2 = BigInteger.valueOf(q * q);
        BigInteger power = one.divide(BigInteger.valueOf(q));
        BigInteger s = power;
        for (long j = 1; power.signum() > 0; j++) {
            power = power.divide(q2);
            BigInteger term = power.divide(BigInteger.valueOf(2 * j + 1));
            s = (alternating && j % 2 == 1) ? s.subtract(term) : s.add(term);
        }
        return s;
    }

    // exp(y) for 0 <= y < 1: halve 16 times, Taylor, square back.
    static BigInteger expFixed(BigInteger y) {
        BigInteger arg = y.shiftRight(16);
        BigInteger acc = one, term = one;
        for (long j = 1; term.signum() > 0; j++) {
            term = mulShift(term, arg).divide(BigInteger.valueOf(j));
            acc = acc.add(term);
        }
        for (int i = 0; i < 16; i++) {
            acc = mulShift(acc, acc);
        }
        return acc;
    }

    // K0 scaled by 2^wp.
    static BigInteger khinchinFixed(int digits) {
        int prec = (int) Math.ceil((digits + 2) * (Math.log(10) / Math.log(2))) + 32;
        wp = prec + (int) Math.ceil(Math.sqrt(prec)) + 64;
        one = BigInteger.ONE.shiftLeft(wp);
        int bigN = Math.max(3, (int) Math.pow(wp, 0.35));
        int m = (int) Math.ceil(wp * Math.log(2) / (2 * Math.log(bigN))) + 1;

        // -ln((k-1)/k) ln((k+1)/k) = 4 atanh(1/(2k-1)) atanh(1/(2k+1)).
        BigInteger s = BigInteger.ZERO;
        for (long k = 2; k < bigN; k++) {
            s = s.add(mulShift(arcRecip(2 * k - 1, false),
                               arcRecip(2 * k + 1, false)).shiftLeft(2));
        }

        // pi = 16 atan(1/5) - 4 atan(1/239); zeta(2) = pi^2/6 seeds the
        // positive-term recurrence, convolutions via parallel streams.
        BigInteger pi = arcRecip(5, true).shiftLeft(4)
                .subtract(arcRecip(239, true).shiftLeft(2));
        BigInteger[] z = new BigInteger[m + 1];
        z[1] = mulShift(pi, pi).divide(BigInteger.valueOf(6));
        for (int n = 2; n <= m; n++) {
            int half = (n - 1) / 2;
            final int nn = n;
            BigInteger raw;
            if (half >= 64) {
                raw = IntStream.rangeClosed(1, half).parallel()
                        .mapToObj(j -> z[j].multiply(z[nn - j]))
                        .reduce(BigInteger.ZERO, BigInteger::add);
            } else {
                raw = BigInteger.ZERO;
                for (int j = 1; j <= half; j++) {
                    raw = raw.add(z[j].multiply(z[n - j]));
                }
            }
            BigInteger acc = raw.shiftRight(wp).shiftLeft(1);
            if (n % 2 == 0) {
                acc = acc.add(mulShift(z[n / 2], z[n / 2]));
            }
            z[n] = acc.shiftLeft(1).divide(BigInteger.valueOf(2L * n + 1));
        }

        // Main accelerated loop, ascending, incremental power table.
        BigInteger[] powers = new BigInteger[bigN - 2];
        for (int k = 2; k < bigN; k++) {
            powers[k - 2] = one.divide(BigInteger.valueOf((long) k * k));
        }
        BigInteger h = one;
        for (long n = 1; n <= m; n++) {
            BigInteger tail = z[(int) n].subtract(one);
            for (BigInteger p : powers) {
                tail = tail.subtract(p);
            }
            s = s.add(tail.multiply(h).divide(BigInteger.valueOf(n)).shiftRight(wp));
            h = h.add(one.divide(BigInteger.valueOf(2 * n + 1)))
                 .subtract(one.divide(BigInteger.valueOf(2 * n)));
            for (int k = 2; k < bigN; k++) {
                powers[k - 2] = powers[k - 2].divide(BigInteger.valueOf((long) k * k));
            }
        }

        BigInteger ln2 = arcRecip(3, false).shiftLeft(1);
        return expFixed(s.shiftLeft(wp).divide(ln2));
    }

    public static void main(String[] args) throws Exception {
        int digits = args.length > 0 ? Integer.parseInt(args[0]) : 100;
        BigInteger value = khinchinFixed(digits);
        BigInteger scaled = value.multiply(BigInteger.TEN.pow(digits))
                .add(BigInteger.ONE.shiftLeft(wp - 1)).shiftRight(wp);
        String text = scaled.toString();
        String result = text.charAt(0) + "." + text.substring(1);
        if (args.length > 1) {
            Files.writeString(Path.of(args[1]), result + "\n");
        } else {
            System.out.println(result);
        }
    }
}
