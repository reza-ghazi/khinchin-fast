// Khinchin's constant (OEIS A002210) via the accelerated zeta series.
//
// Same mathematics as ../khinchin_fast.c:
//
//	ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
//	               + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
//
//	h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
//
// Runs entirely on math/big.Int in fixed point (scale 2^wp) — the
// "native bignum" data point for Go, no dependencies. math/big has no
// transcendentals, so the port carries its own fixed-point kit: pi by
// Machin's formula, logarithms via the atanh(1/q) series, exp by
// argument-halving plus Taylor. The even zeta values come from the
// positive-term recurrence (n + 1/2) zeta(2n) = sum_j zeta(2j)
// zeta(2n-2j) — O(M^2) bignum products, with the convolutions (halved
// via their j <-> n-j symmetry) fanned out over goroutines. Guard bits
// absorb the drift; nothing here is interval-certified.
//
// Run:  go run khinchin.go DIGITS [OUTPUT_FILE]        (default 100)
//
// With OUTPUT_FILE the decimal value plus a newline is written there
// (the same contract as ../khinchin_fast.c); otherwise stdout.
package main

import (
	"fmt"
	"math"
	"math/big"
	"os"
	"runtime"
	"strconv"
	"sync"
)

var (
	wp  uint
	one *big.Int
)

func mulShift(a, b *big.Int) *big.Int {
	return new(big.Int).Rsh(new(big.Int).Mul(a, b), wp)
}

// atanh(1/q) = sum_{j>=0} 1 / ((2j+1) q^(2j+1)); atan(1/q) alternates.
func arcRecip(q int64, alternating bool) *big.Int {
	q2 := big.NewInt(q * q)
	power := new(big.Int).Div(one, big.NewInt(q))
	s := new(big.Int).Set(power)
	term := new(big.Int)
	for j := int64(1); power.Sign() > 0; j++ {
		power.Div(power, q2)
		term.Div(power, big.NewInt(2*j+1))
		if alternating && j%2 == 1 {
			s.Sub(s, term)
		} else {
			s.Add(s, term)
		}
	}
	return s
}

// exp(y) for 0 <= y < 1: halve `shift` times, Taylor, square back.
func expFixed(y *big.Int) *big.Int {
	const shift = 16
	arg := new(big.Int).Rsh(y, shift)
	acc := new(big.Int).Set(one)
	term := new(big.Int).Set(one)
	for j := int64(1); term.Sign() > 0; j++ {
		term = mulShift(term, arg)
		term.Div(term, big.NewInt(j))
		acc.Add(acc, term)
	}
	for i := 0; i < shift; i++ {
		acc = mulShift(acc, acc)
	}
	return acc
}

// khinchinFixed returns K0 scaled by 2^wp.
func khinchinFixed(digits int) *big.Int {
	prec := uint(math.Ceil(float64(digits+2)*math.Log2(10))) + 32
	wp = prec + uint(math.Ceil(math.Sqrt(float64(prec)))) + 64
	one = new(big.Int).Lsh(big.NewInt(1), wp)
	bigN := int(math.Max(3, math.Floor(math.Pow(float64(wp), 0.35))))
	m := int(math.Ceil(float64(wp)*math.Ln2/(2*math.Log(float64(bigN))))) + 1

	// Finite logarithmic correction:
	// -ln((k-1)/k) ln((k+1)/k) = 4 atanh(1/(2k-1)) atanh(1/(2k+1)).
	s := new(big.Int)
	for k := int64(2); k < int64(bigN); k++ {
		p := mulShift(arcRecip(2*k-1, false), arcRecip(2*k+1, false))
		s.Add(s, p.Lsh(p, 2))
	}

	// pi = 16 atan(1/5) - 4 atan(1/239); zeta(2) = pi^2/6 seeds the
	// positive-term recurrence (n + 1/2) z_n = sum z_j z_{n-j}, whose
	// convolutions run in parallel over goroutines.
	pi := new(big.Int).Sub(
		new(big.Int).Lsh(arcRecip(5, true), 4),
		new(big.Int).Lsh(arcRecip(239, true), 2))
	workers := runtime.NumCPU()
	z := make([]*big.Int, m+1)
	z[1] = new(big.Int).Div(mulShift(pi, pi), big.NewInt(6))
	for n := 2; n <= m; n++ {
		half := (n - 1) / 2
		raw := new(big.Int)
		if half >= 64 {
			partials := make([]*big.Int, workers)
			var wg sync.WaitGroup
			for w := 0; w < workers; w++ {
				wg.Add(1)
				go func(w int) {
					defer wg.Done()
					local := new(big.Int)
					product := new(big.Int)
					for j := 1 + half*w/workers; j <= half*(w+1)/workers; j++ {
						local.Add(local, product.Mul(z[j], z[n-j]))
					}
					partials[w] = local
				}(w)
			}
			wg.Wait()
			for _, p := range partials {
				raw.Add(raw, p)
			}
		} else {
			product := new(big.Int)
			for j := 1; j <= half; j++ {
				raw.Add(raw, product.Mul(z[j], z[n-j]))
			}
		}
		acc := new(big.Int).Lsh(new(big.Int).Rsh(raw, wp), 1)
		if n%2 == 0 {
			acc.Add(acc, mulShift(z[n/2], z[n/2]))
		}
		z[n] = acc.Div(acc.Lsh(acc, 1), big.NewInt(int64(2*n+1)))
	}

	// Main accelerated loop, ascending, incremental power table.
	powers := make([]*big.Int, 0, bigN-2)
	for k := int64(2); k < int64(bigN); k++ {
		powers = append(powers, new(big.Int).Div(one, big.NewInt(k*k)))
	}
	h := new(big.Int).Set(one)
	tail := new(big.Int)
	for n := int64(1); n <= int64(m); n++ {
		tail.Sub(z[n], one)
		for _, p := range powers {
			tail.Sub(tail, p)
		}
		term := new(big.Int).Mul(tail, h)
		term.Div(term, big.NewInt(n))
		s.Add(s, term.Rsh(term, wp))
		h.Add(h, new(big.Int).Div(one, big.NewInt(2*n+1)))
		h.Sub(h, new(big.Int).Div(one, big.NewInt(2*n)))
		for i, k := 0, int64(2); k < int64(bigN); i, k = i+1, k+1 {
			powers[i].Div(powers[i], big.NewInt(k*k))
		}
	}

	ln2 := new(big.Int).Lsh(arcRecip(3, false), 1)
	y := new(big.Int).Div(new(big.Int).Lsh(s, wp), ln2)
	return expFixed(y)
}

func main() {
	digits := 100
	if len(os.Args) > 1 {
		digits, _ = strconv.Atoi(os.Args[1])
	}
	value := khinchinFixed(digits)
	scaled := new(big.Int).Mul(value, new(big.Int).Exp(
		big.NewInt(10), big.NewInt(int64(digits)), nil))
	scaled.Add(scaled, new(big.Int).Lsh(big.NewInt(1), wp-1))
	scaled.Rsh(scaled, wp)
	text := scaled.String()
	result := text[:1] + "." + text[1:]
	if len(os.Args) > 2 {
		if err := os.WriteFile(os.Args[2], []byte(result+"\n"), 0o644); err != nil {
			panic(err)
		}
	} else {
		fmt.Println(result)
	}
}
