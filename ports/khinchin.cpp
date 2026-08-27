// Khinchin's constant (OEIS A002210) via the accelerated zeta series.
//
// Same mathematics and structure as ../khinchin_fast.c, including its
// two-region split and dropping precision:
//
// - Bernoulli region (n < n_direct): zeta(2n) reconstructed from exact
//   Bernoulli numbers streamed by FLINT's reverse Bernoulli iterator,
//   one iterator per block (C++ can use bernoulli_rev_t natively - no
//   shim needed, unlike the Rust port).
// - Direct region (n >= n_direct): literal tail sums over k in [N, K],
//   K ~ 2^(wp/2n) <= ~8N, no Bernoulli numbers at all.
// - A term of magnitude N^(-2n) only needs ~wp - 2n log2 N accurate
//   bits, and table entries start at the precision of their final, most
//   demanding use.
//
// Arithmetic is Arb through a small RAII wrapper; blocks fan out over
// std::thread, heaviest first via an atomic work index. The final ball's
// radius is tracked by Arb throughout, but unlike the C program the
// decimal output is midpoint-based, not certified.
//
// Build and run:
//   g++ -O3 -march=native -std=c++17 -o khinchin-cpp khinchin.cpp -lflint -lmpfr -lgmp -pthread
//   ./khinchin-cpp DIGITS [OUTPUT_FILE]     (default 100)
//
// With OUTPUT_FILE the decimal value plus a newline is written there
// (the same contract as ../khinchin_fast.c); otherwise stdout.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include <mpfr.h>

#include <flint/arb.h>
#include <flint/bernoulli.h>
#include <flint/fmpz.h>

namespace {

// Minimal RAII wrapper over arb_t.
class Arb {
  public:
    Arb() { arb_init(value_); }
    explicit Arb(ulong n) {
        arb_init(value_);
        arb_set_ui(value_, n);
    }
    ~Arb() { arb_clear(value_); }
    Arb(const Arb&) = delete;
    Arb& operator=(const Arb&) = delete;
    Arb(Arb&&) = delete;
    arb_ptr get() { return value_; }
    arb_srcptr get() const { return value_; }

  private:
    arb_t value_;
};

slong clamp_prec(double bits, slong wp) {
    return static_cast<slong>(std::max(64.0, std::min(bits, double(wp))));
}

// h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j, accumulated at prec.
void alternating_harmonic(Arb& h, ulong n, slong prec) {
    Arb t;
    arb_zero(h.get());
    for (ulong j = 1; j <= 2 * n - 1; ++j) {
        arb_one(t.get());
        arb_div_ui(t.get(), t.get(), j, prec);
        if (j % 2 == 1)
            arb_add(h.get(), h.get(), t.get(), prec);
        else
            arb_sub(h.get(), h.get(), t.get(), prec);
    }
}

struct Params {
    slong wp;
    ulong N;
    double log2_N;
    ulong M;
    ulong n_direct;
};

// Accelerated terms for n in [first, last] via the reverse Bernoulli
// iterator, descending, with dropping precision as in the C program.
void bernoulli_block(Arb& sum, ulong first, ulong last, const Params& p) {
    const slong acc_prec =
        clamp_prec(p.wp - 2.0 * first * p.log2_N + 64.0, p.wp);
    bernoulli_rev_t iter;
    bernoulli_rev_init(iter, 2 * last);
    fmpz_t num, den;
    fmpz_init(num);
    fmpz_init(den);

    Arb factor, tau2, inv_tau2, zeta, tail, h, t;
    arb_const_pi(tau2.get(), p.wp);
    arb_mul_2exp_si(tau2.get(), tau2.get(), 1);
    arb_mul(tau2.get(), tau2.get(), tau2.get(), p.wp);
    arb_inv(inv_tau2.get(), tau2.get(), p.wp);
    arb_pow_ui(factor.get(), tau2.get(), last, p.wp);
    arb_fac_ui(t.get(), 2 * last, p.wp);
    arb_div(factor.get(), factor.get(), t.get(), p.wp);
    arb_mul_2exp_si(factor.get(), factor.get(), -1);

    std::vector<Arb> powers(p.N - 2);
    std::vector<slong> table_prec(p.N - 2);
    for (ulong k = 2; k < p.N; ++k) {
        table_prec[k - 2] = clamp_prec(
            p.wp - 2.0 * first * std::log2(double(k)) + 32.0, p.wp);
        arb_set_ui(powers[k - 2].get(), k);
        arb_pow_ui(powers[k - 2].get(), powers[k - 2].get(), 2 * last,
                   table_prec[k - 2]);
        arb_inv(powers[k - 2].get(), powers[k - 2].get(), table_prec[k - 2]);
    }
    alternating_harmonic(h, last, acc_prec);
    arb_zero(sum.get());

    for (ulong n = last;; --n) {
        const slong tail_prec =
            clamp_prec(double(p.wp) - 2.0 * n + 32.0, p.wp);
        bernoulli_rev_next(num, den, iter);
        arb_set_round_fmpz(zeta.get(), num, p.wp);
        arb_div_fmpz(zeta.get(), zeta.get(), den, p.wp);
        arb_mul(zeta.get(), zeta.get(), factor.get(), p.wp);
        arb_abs(zeta.get(), zeta.get());
        arb_sub_ui(tail.get(), zeta.get(), 1, tail_prec);
        for (ulong k = 2; k < p.N; ++k)
            arb_sub(tail.get(), tail.get(), powers[k - 2].get(), tail_prec);
        arb_mul(tail.get(), tail.get(), h.get(), tail_prec);
        arb_div_ui(tail.get(), tail.get(), n, tail_prec);
        arb_add(sum.get(), sum.get(), tail.get(), acc_prec);
        if (n == first) break;
        arb_mul_ui(factor.get(), factor.get(), 2 * n * (2 * n - 1), p.wp);
        arb_mul(factor.get(), factor.get(), inv_tau2.get(), p.wp);
        for (ulong k = 2; k < p.N; ++k)
            arb_mul_ui(powers[k - 2].get(), powers[k - 2].get(), k * k,
                       table_prec[k - 2]);
        arb_one(t.get());
        arb_div_ui(t.get(), t.get(), (2 * n - 2) * (2 * n - 1), acc_prec);
        arb_add(h.get(), h.get(), t.get(), acc_prec);
    }
    fmpz_clear(den);
    fmpz_clear(num);
    bernoulli_rev_clear(iter);
}

// Accelerated terms for n in [first, last] as literal tail sums over
// k in [N, K], descending, with dropping precision.
void direct_block(Arb& sum, ulong first, ulong last, const Params& p) {
    ulong K = std::min<ulong>(
        ulong(std::pow(2.0, (p.wp + 32.0) / (2.0 * first))) + 1,
        16 * p.N + 64);
    K = std::max(K, p.N);
    const slong acc_prec =
        clamp_prec(p.wp - 2.0 * first * p.log2_N + 64.0, p.wp);

    std::vector<Arb> powers(K - p.N + 1);
    std::vector<slong> table_prec(K - p.N + 1);
    for (ulong k = p.N; k <= K; ++k) {
        table_prec[k - p.N] = clamp_prec(
            p.wp - 2.0 * first * std::log2(double(k)) + 32.0, p.wp);
        arb_set_ui(powers[k - p.N].get(), k);
        arb_pow_ui(powers[k - p.N].get(), powers[k - p.N].get(), 2 * last,
                   table_prec[k - p.N]);
        arb_inv(powers[k - p.N].get(), powers[k - p.N].get(),
                table_prec[k - p.N]);
    }
    Arb h, tail, t;
    alternating_harmonic(h, last, acc_prec);
    arb_zero(sum.get());

    for (ulong n = last;; --n) {
        const slong sum_prec =
            clamp_prec(p.wp - 2.0 * n * p.log2_N + 32.0, p.wp);
        arb_zero(tail.get());
        for (ulong k = p.N; k <= K; ++k)
            arb_add(tail.get(), tail.get(), powers[k - p.N].get(), sum_prec);
        arb_mul(tail.get(), tail.get(), h.get(), sum_prec);
        arb_div_ui(tail.get(), tail.get(), n, sum_prec);
        arb_add(sum.get(), sum.get(), tail.get(), acc_prec);
        if (n == first) break;
        for (ulong k = p.N; k <= K; ++k)
            arb_mul_ui(powers[k - p.N].get(), powers[k - p.N].get(), k * k,
                       table_prec[k - p.N]);
        arb_one(t.get());
        arb_div_ui(t.get(), t.get(), (2 * n - 2) * (2 * n - 1), acc_prec);
        arb_add(h.get(), h.get(), t.get(), acc_prec);
    }
}

std::string khinchin(ulong digits) {
    Params p;
    const slong prec =
        slong(std::ceil((digits + 2) * std::log2(10.0))) + 32;
    p.wp = prec + slong(std::ceil(std::sqrt(double(prec)))) + 64;
    p.N = std::max<ulong>(3, ulong(std::pow(double(p.wp), 0.35)));
    p.log2_N = std::log2(double(p.N));
    p.M = ulong(std::ceil(p.wp * std::log(2.0) /
                          (2.0 * std::log(double(p.N))))) + 1;
    p.n_direct = std::min<ulong>(
        ulong(std::ceil(p.wp / (2.0 * (p.log2_N + 3.0)))), p.M + 1);

    Arb S, t1, t2;
    for (ulong k = 2; k < p.N; ++k) {
        arb_set_ui(t1.get(), k - 1);
        arb_div_ui(t1.get(), t1.get(), k, p.wp);
        arb_log(t1.get(), t1.get(), p.wp);
        arb_set_ui(t2.get(), k + 1);
        arb_div_ui(t2.get(), t2.get(), k, p.wp);
        arb_log(t2.get(), t2.get(), p.wp);
        arb_mul(t1.get(), t1.get(), t2.get(), p.wp);
        arb_sub(S.get(), S.get(), t1.get(), p.wp);
    }

    // Equal-count blocks per region, heaviest (high-n Bernoulli) first,
    // pulled by worker threads from an atomic index.
    struct Block { bool bernoulli; ulong first, last; };
    std::vector<Block> blocks;
    const ulong threads =
        std::max(1u, std::thread::hardware_concurrency());
    const ulong bt = p.n_direct - 1;
    const ulong bb = std::max<ulong>(1, std::min(4 * threads, bt / 8));
    for (ulong b = bb; b-- > 0;) {
        const ulong first = 1 + bt * b / bb, last = bt * (b + 1) / bb;
        if (first <= last) blocks.push_back({true, first, last});
    }
    if (p.M >= p.n_direct) {
        const ulong dt = p.M - p.n_direct + 1;
        const ulong db = std::max<ulong>(1, std::min(4 * threads, dt / 8));
        for (ulong b = 0; b < db; ++b) {
            const ulong first = p.n_direct + dt * b / db;
            const ulong last = p.n_direct + dt * (b + 1) / db - 1;
            if (first <= last) blocks.push_back({false, first, last});
        }
    }
    std::vector<Arb> partials(blocks.size());
    std::atomic<size_t> next{0};
    std::vector<std::thread> pool;
    for (ulong w = 0; w < threads; ++w)
        pool.emplace_back([&] {
            for (size_t i; (i = next.fetch_add(1)) < blocks.size();) {
                const Block& blk = blocks[i];
                if (blk.bernoulli)
                    bernoulli_block(partials[i], blk.first, blk.last, p);
                else
                    direct_block(partials[i], blk.first, blk.last, p);
            }
            flint_cleanup();
        });
    for (auto& th : pool) th.join();
    for (auto& part : partials)
        arb_add(S.get(), S.get(), part.get(), p.wp);

    arb_log_ui(t1.get(), 2, p.wp);
    arb_div(S.get(), S.get(), t1.get(), p.wp);
    arb_exp(S.get(), S.get(), p.wp);

    // Midpoint out through MPFR, correctly rounded decimal string.
    mpfr_t lo, hi;
    mpfr_init2(lo, p.wp);
    mpfr_init2(hi, p.wp);
    arb_get_interval_mpfr(lo, hi, S.get());
    mpfr_exp_t expo;
    char* str = mpfr_get_str(nullptr, &expo, 10, digits + 1, lo, MPFR_RNDN);
    std::string text(str, str + 1);
    text += ".";
    text.append(str + 1, str + digits + 1);
    mpfr_free_str(str);
    mpfr_clear(hi);
    mpfr_clear(lo);
    return text;
}

}  // namespace

int main(int argc, char** argv) {
    const ulong digits = argc > 1 ? std::stoul(argv[1]) : 100;
    const std::string result = khinchin(digits);
    if (argc > 2) {
        std::ofstream out(argv[2]);
        out << result << "\n";
    } else {
        std::cout << result << "\n";
    }
    flint_cleanup();
    return 0;
}
