#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>

#include <flint/arb.h>
#include <flint/bernoulli.h>
#include <flint/flint.h>
#include <flint/fmpz.h>

static double seconds_since(const struct timespec *start)
{
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (double) (now.tv_sec - start->tv_sec)
        + 1e-9 * (double) (now.tv_nsec - start->tv_nsec);
}

static void usage(FILE *stream, const char *program)
{
    fprintf(stream,
        "Usage: %s DIGITS OUTPUT_FILE\n"
        "       %s --force DIGITS OUTPUT_FILE\n"
        "       %s --benchmark DIGITS\n\n"
        "DIGITS is the number of digits after the decimal point.\n"
        "The normal mode writes a correctly rounded decimal value to OUTPUT_FILE.\n"
        "--force overrides the physical-memory preflight check.\n",
        program, program, program);
}

static int parse_digits(const char *text, ulong *digits)
{
    char *end = NULL;
    unsigned long long value;

    if (text[0] == '\0' || text[0] == '-')
        return 0;

    errno = 0;
    value = strtoull(text, &end, 10);
    if (errno != 0 || *end != '\0' || value > ULONG_MAX)
        return 0;

    *digits = (ulong) value;
    return 1;
}

static slong working_precision(ulong digits)
{
    /* ceil(log2(10) * (digits + 1)) plus ample rounding guard bits. */
    const long double log2_10 = 3.3219280948873623478703194294893902L;
    long double bits = log2_10 * ((long double) digits + 1.0L) + 128.0L;

    if (bits > (long double) WORD_MAX)
    {
        fprintf(stderr, "Requested precision is too large for this build.\n");
        exit(EXIT_FAILURE);
    }

    return (slong) bits + 1;
}

static long double estimated_peak_bytes(ulong digits, slong precision,
    int threads)
{
    ulong N = (ulong) pow((double) precision, 0.35);
    long double number_bytes = (long double) precision / 8.0L + 64.0L;
    long double parallel_tables;
    long double decimal_conversion;

    if (N < 2)
        N = 2;
    parallel_tables = 2.0L * threads * (N - 2) * number_bytes;
    decimal_conversion = (long double) digits
        + 3.0L * (long double) precision / 8.0L;
    return 2.0L * (parallel_tables + decimal_conversion)
        + 64.0L * 1024.0L * 1024.0L;
}

static long double physical_memory_bytes(void)
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGESIZE);
    if (pages <= 0 || page_size <= 0)
        return 0.0L;
    return (long double) pages * (long double) page_size;
}

/*
   Bailey-Borwein-Crandall acceleration, using the same rigorous tail bound as
   FLINT's arb_const_khinchin(), with independent zeta terms split by thread.
*/
static void khinchin_parallel(arb_t result, slong precision, int threads)
{
    ulong N = (ulong) pow((double) precision, 0.35);
    ulong M;
    arb_ptr partials;
    arb_t t, u;

    if (N < 2)
        N = 2;
    M = (ulong) ceil(((double) precision * 0.69314718055994530942)
        / (2.0 * log((double) N)));
    if (M < 1)
        M = 1;
    if ((ulong) threads > M)
        threads = (int) M;

    arb_init(t);
    arb_init(u);
    arb_zero(result);

    /* Finite logarithmic correction. */
    for (ulong k = 2; k < N; ++k)
    {
        arb_set_ui(t, k - 1);
        arb_div_ui(t, t, k, precision);
        arb_log(t, t, precision);

        arb_set_ui(u, k + 1);
        arb_div_ui(u, u, k, precision);
        arb_log(u, u, precision);

        arb_mul(t, t, u, precision);
        arb_sub(result, result, t, precision);
    }

    partials = _arb_vec_init(threads);
    omp_set_dynamic(0);

#pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        ulong first = 1 + (M * (ulong) tid) / (ulong) threads;
        ulong last = (M * (ulong) (tid + 1)) / (ulong) threads;
        arb_ptr powers = _arb_vec_init(N - 2);
        arb_t zeta_tail, h, term, reciprocal;

        arb_init(zeta_tail);
        arb_init(h);
        arb_init(term);
        arb_init(reciprocal);
        arb_one(h);

        /* h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j. */
        for (ulong n = 1; n < first; ++n)
        {
            arb_set_ui(reciprocal, 2 * n);
            arb_mul_ui(reciprocal, reciprocal, 2 * n + 1, precision);
            arb_inv(reciprocal, reciprocal, precision);
            arb_sub(h, h, reciprocal, precision);
        }

        /* Start each power at k^(-2(first-1)); the loop advances it once. */
        for (ulong k = 2; k < N; ++k)
        {
            arb_set_ui(powers + (k - 2), k);
            arb_mul(powers + (k - 2), powers + (k - 2),
                powers + (k - 2), precision);
            if (first > 1)
            {
                arb_pow_ui(powers + (k - 2), powers + (k - 2), first - 1,
                    precision);
                arb_inv(powers + (k - 2), powers + (k - 2), precision);
            }
            else
            {
                arb_one(powers + (k - 2));
            }
        }

        arb_zero(partials + tid);
        for (ulong n = first; n <= last; ++n)
        {
            arb_zeta_ui(zeta_tail, 2 * n, precision);
            arb_sub_ui(zeta_tail, zeta_tail, 1, precision);

            for (ulong k = 2; k < N; ++k)
            {
                arb_div_ui(powers + (k - 2), powers + (k - 2), k * k,
                    precision);
                arb_sub(zeta_tail, zeta_tail, powers + (k - 2), precision);
            }

            arb_div_ui(term, zeta_tail, n, precision);
            arb_mul(term, term, h, precision);
            arb_add(partials + tid, partials + tid, term, precision);

            arb_set_ui(reciprocal, 2 * n);
            arb_mul_ui(reciprocal, reciprocal, 2 * n + 1, precision);
            arb_inv(reciprocal, reciprocal, precision);
            arb_sub(h, h, reciprocal, precision);
        }

        arb_clear(reciprocal);
        arb_clear(term);
        arb_clear(h);
        arb_clear(zeta_tail);
        _arb_vec_clear(powers, N - 2);
    }

    for (int i = 0; i < threads; ++i)
        arb_add(result, result, partials + i, precision);

    /* The omitted positive tail is bounded by N^(-2M). */
    arb_set_ui(t, N);
    arb_pow_ui(t, t, 2 * M, MAG_BITS);
    arb_inv(t, t, MAG_BITS);
    arb_add_error(result, t);

    arb_log_ui(t, 2, precision);
    arb_div(result, result, t, precision);
    arb_exp(result, result, precision);

    _arb_vec_clear(partials, threads);
    arb_clear(u);
    arb_clear(t);
}

static slong clamp_precision(double bits, slong maximum)
{
    if (bits < 64.0)
        return 64;
    if (bits >= (double) maximum)
        return maximum;
    return (slong) bits;
}

/* Estimated relative cost of accelerated term n, in bits of linear
   bignum work. Used only for load balancing. */
static double term_cost(ulong n, ulong n_direct, ulong N, slong precision,
    double log2_N)
{
    double table_bits = fmax(64.0,
        (double) precision - 2.0 * (double) n * log2_N);

    if (n >= n_direct)
    {
        /* Direct tail sum over k in [N, 2^(precision/2n)]. */
        double K = pow(2.0,
            ((double) precision + 32.0) / (2.0 * (double) n));
        double entries = fmax(1.0, K - (double) N + 1.0);

        return entries * (1.25 * table_bits + 160.0) + 15.0 * table_bits;
    }
    else
    {
        /* Power-table update plus the reverse Bernoulli iterator's own
           Dirichlet sum (working precision of the exact B_2n) plus the
           full-precision Bernoulli-to-zeta conversion. */
        double s = 2.0 * (double) n;
        double iterator_precision = s * fmax(1.0, log2(s) - 4.09);
        double zeta_terms = fmax(1.0,
            pow(2.0, iterator_precision / (s - 1.0)));
        double bernoulli_bits = zeta_terms * fmax(64.0,
            iterator_precision - s * (log2(zeta_terms) - 1.44));

        return (double) N * table_bits + bernoulli_bits
            + 24.0 * (double) precision;
    }
}

/* h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j = psi(2n) - psi(n). */
static void alternating_harmonic(arb_t h, ulong n, slong precision)
{
    arb_t x;

    arb_init(x);
    arb_set_ui(x, 2 * n);
    arb_digamma(h, x, precision);
    arb_set_ui(x, n);
    arb_digamma(x, x, precision);
    arb_sub(h, h, x, precision);
    arb_clear(x);
}

/*
   Direct evaluation of the accelerated terms for a block of large n:
   zeta(2n) - 1 - sum_{2<=k<N} k^(-2n) = sum_{k>=N} k^(-2n), which to
   absolute accuracy 2^(-precision) needs only k up to about
   2^(precision/2n). The term is a short positive sum at low precision:
   no Bernoulli numbers, no cancellation, no full-precision arithmetic.

   Term n of the accelerated sum has magnitude about N^(-2n), so it only
   needs about precision - 2n*log2(N) accurate bits. Power table entry k
   is kept at the precision demanded by the lowest n in the block, since
   relative errors persist through the incremental updates.
*/
static void scp_direct_range(arb_t accumulator, ulong first, ulong last,
    ulong N, slong precision, double log2_N)
{
    ulong K = (ulong) pow(2.0,
        ((double) precision + 32.0) / (2.0 * (double) first)) + 1;
    ulong count;
    /* Safety clamp: with the default n_direct, K <= about 8N. If a
       misconfigured split ever asks for more, the truncation bound below
       stays rigorous and certification fails loudly instead. */
    ulong K_limit = 16 * N + 64;
    ulong exponent = 2 * last;
    slong accumulator_precision = clamp_precision(
        (double) precision - 2.0 * (double) first * log2_N + 64.0,
        precision);
    arb_ptr powers;
    slong *table_precision;
    arb_t h, zeta_tail, term, reciprocal, bound;

    if (K < N)
        K = N;
    if (K > K_limit)
        K = K_limit;
    count = K - N + 1;
    powers = _arb_vec_init(count);
    table_precision = flint_malloc(count * sizeof(slong));

    arb_init(h);
    arb_init(zeta_tail);
    arb_init(term);
    arb_init(reciprocal);
    arb_init(bound);

    alternating_harmonic(h, last, accumulator_precision);

    for (ulong k = N; k <= K; ++k)
    {
        table_precision[k - N] = clamp_precision(
            (double) precision - 2.0 * (double) first * log2((double) k)
            + 32.0, precision);
        arb_set_ui(powers + (k - N), k);
        arb_pow_ui(powers + (k - N), powers + (k - N), exponent,
            table_precision[k - N]);
        arb_inv(powers + (k - N), powers + (k - N),
            table_precision[k - N]);
    }

    for (ulong n = last; n >= first; --n)
    {
        slong sum_precision = clamp_precision(
            (double) precision - 2.0 * (double) n * log2_N + 32.0,
            precision);

        exponent = 2 * n;
        arb_zero(zeta_tail);
        for (ulong i = 0; i < count; ++i)
            arb_add(zeta_tail, zeta_tail, powers + i, sum_precision);

        /* Rigorous bound for the omitted k > K:
           sum_{k>K} k^(-2n) <= (K+1)^(-2n) (1 + (K+1)/(2n-1)). */
        arb_set_ui(bound, K + 1);
        arb_pow_ui(bound, bound, exponent, 32);
        arb_inv(bound, bound, 32);
        arb_set_ui(reciprocal, K + 1);
        arb_div_ui(reciprocal, reciprocal, exponent - 1, 32);
        arb_add_ui(reciprocal, reciprocal, 1, 32);
        arb_mul(bound, bound, reciprocal, 32);
        arb_add_error(zeta_tail, bound);

        arb_div_ui(term, zeta_tail, n, sum_precision);
        arb_mul(term, term, h, sum_precision);
        arb_add(accumulator, accumulator, term, accumulator_precision);

        if (n > first)
        {
            for (ulong k = N; k <= K; ++k)
                arb_mul_ui(powers + (k - N), powers + (k - N), k * k,
                    table_precision[k - N]);

            arb_set_ui(reciprocal, 2 * (n - 1));
            arb_mul_ui(reciprocal, reciprocal, 2 * n - 1,
                accumulator_precision);
            arb_inv(reciprocal, reciprocal, accumulator_precision);
            arb_add(h, h, reciprocal, accumulator_precision);
        }

        if (n == first)
            break;
    }

    arb_clear(bound);
    arb_clear(reciprocal);
    arb_clear(term);
    arb_clear(zeta_tail);
    arb_clear(h);
    flint_free(table_precision);
    _arb_vec_clear(powers, count);
}

/*
   Bernoulli-number evaluation of the accelerated terms for a block of
   small n, where the direct tail sum would need too many k. FLINT's
   reverse Bernoulli iterator shares one power table across the descending
   block of even zeta values. zeta(2n) is near 1, so reconstructing it
   from B_2n requires full precision even though the term itself does not;
   the power table uses the same per-entry dropping precision as above.
*/
static void scp_bernoulli_range(arb_t accumulator, ulong first, ulong last,
    ulong N, slong precision, double log2_N)
{
    ulong exponent = 2 * last;
    slong accumulator_precision = clamp_precision(
        (double) precision - 2.0 * (double) first * log2_N + 64.0,
        precision);
    bernoulli_rev_t iterator;
    arb_ptr powers = _arb_vec_init(N - 2);
    slong *table_precision = flint_malloc((N > 2 ? N - 2 : 1)
        * sizeof(slong));
    arb_t factor, tau2, inv_tau2, zeta_tail, h, term, reciprocal, x;
    fmpz_t numerator, denominator;

    arb_init(factor);
    arb_init(tau2);
    arb_init(inv_tau2);
    arb_init(zeta_tail);
    arb_init(h);
    arb_init(term);
    arb_init(reciprocal);
    arb_init(x);
    fmpz_init(numerator);
    fmpz_init(denominator);

    bernoulli_rev_init(iterator, exponent);

    arb_const_pi(tau2, precision);
    arb_mul_2exp_si(tau2, tau2, 1);
    arb_mul(tau2, tau2, tau2, precision);
    arb_inv(inv_tau2, tau2, precision);

    arb_pow_ui(factor, tau2, last, precision);
    arb_fac_ui(x, exponent, precision);
    arb_div(factor, factor, x, precision);
    arb_mul_2exp_si(factor, factor, -1);

    alternating_harmonic(h, last, accumulator_precision);

    for (ulong k = 2; k < N; ++k)
    {
        table_precision[k - 2] = clamp_precision(
            (double) precision - 2.0 * (double) first * log2((double) k)
            + 32.0, precision);
        arb_set_ui(powers + (k - 2), k);
        arb_pow_ui(powers + (k - 2), powers + (k - 2), exponent,
            table_precision[k - 2]);
        arb_inv(powers + (k - 2), powers + (k - 2), table_precision[k - 2]);
    }

    for (ulong n = last; n >= first; --n)
    {
        slong tail_precision = clamp_precision(
            (double) (precision - 2 * (slong) n) + 32.0, precision);

        exponent = 2 * n;
        bernoulli_rev_next(numerator, denominator, iterator);

        arb_set_round_fmpz(zeta_tail, numerator, precision);
        arb_div_fmpz(zeta_tail, zeta_tail, denominator, precision);
        arb_mul(zeta_tail, zeta_tail, factor, precision);
        arb_abs(zeta_tail, zeta_tail);
        arb_sub_ui(zeta_tail, zeta_tail, 1, tail_precision);

        for (ulong k = 2; k < N; ++k)
            arb_sub(zeta_tail, zeta_tail, powers + (k - 2), tail_precision);

        arb_div_ui(term, zeta_tail, n, tail_precision);
        arb_mul(term, term, h, tail_precision);
        arb_add(accumulator, accumulator, term, accumulator_precision);

        if (n > first)
        {
            arb_mul_ui(factor, factor, exponent, precision);
            arb_mul_ui(factor, factor, exponent - 1, precision);
            arb_mul(factor, factor, inv_tau2, precision);

            for (ulong k = 2; k < N; ++k)
                arb_mul_ui(powers + (k - 2), powers + (k - 2), k * k,
                    table_precision[k - 2]);

            arb_set_ui(reciprocal, 2 * (n - 1));
            arb_mul_ui(reciprocal, reciprocal, 2 * n - 1,
                accumulator_precision);
            arb_inv(reciprocal, reciprocal, accumulator_precision);
            arb_add(h, h, reciprocal, accumulator_precision);
        }

        if (n == first)
            break;
    }

    fmpz_clear(denominator);
    fmpz_clear(numerator);
    arb_clear(x);
    arb_clear(reciprocal);
    arb_clear(term);
    arb_clear(h);
    arb_clear(zeta_tail);
    arb_clear(inv_tau2);
    arb_clear(tau2);
    arb_clear(factor);
    flint_free(table_precision);
    _arb_vec_clear(powers, N - 2);
    bernoulli_rev_clear(iterator);
}

/*
   SCP-style backend. The accelerated zeta range splits at n_direct: below
   it workers run the reverse Bernoulli iterator, above it the direct tail
   sum. Blocks are sized for equal estimated cost, not equal term counts.
*/
static void khinchin_scp(arb_t result, slong precision, int threads)
{
    ulong N = (ulong) pow((double) precision, 0.35);
    ulong M, n_direct;
    double log2_N;
    int blocks;
    ulong *block_first, *block_last;
    int *block_order;
    double *block_weight, *block_seconds;
    arb_ptr partials;
    arb_t t;
    const char *n_override = getenv("KHINCHIN_N");

    if (n_override != NULL && n_override[0] != '\0')
    {
        ulong parsed;
        if (parse_digits(n_override, &parsed) && parsed >= 2)
            N = parsed;
    }
    if (N < 2)
        N = 2;
    M = (ulong) ceil(((double) precision * 0.69314718055994530942)
        / (2.0 * log((double) N)));
    if (M < 1)
        M = 1;
    if ((ulong) threads > M)
        threads = (int) M;
    log2_N = log2((double) N);

    /* Above n_direct the direct tail sum needs at most about 8N entries. */
    n_direct = (ulong) ceil((double) precision / (2.0 * (log2_N + 3.0)));
    if (n_direct < 1)
        n_direct = 1;
    n_override = getenv("KHINCHIN_DIRECT");
    if (n_override != NULL && n_override[0] != '\0')
    {
        ulong parsed;
        if (parse_digits(n_override, &parsed) && parsed >= 1)
            n_direct = parsed;
    }

    arb_init(t);
    arb_zero(result);

    partials = _arb_vec_init(threads);
    omp_set_dynamic(0);

    /* Finite logarithmic correction, strided across threads. */
#pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        arb_t a, b;

        arb_init(a);
        arb_init(b);
        arb_zero(partials + tid);
        for (ulong k = 2 + (ulong) tid; k < N; k += (ulong) threads)
        {
            arb_set_ui(a, k - 1);
            arb_div_ui(a, a, k, precision);
            arb_log(a, a, precision);

            arb_set_ui(b, k + 1);
            arb_div_ui(b, b, k, precision);
            arb_log(b, b, precision);

            arb_mul(a, a, b, precision);
            arb_sub(partials + tid, partials + tid, a, precision);
        }
        arb_clear(b);
        arb_clear(a);
    }
    for (int i = 0; i < threads; ++i)
        arb_add(result, result, partials + i, precision);

    /* Equal-estimated-cost contiguous blocks over n = 1..M. Blocks
       outnumber threads and are scheduled dynamically, heaviest first, so
       errors in the cost model only misplace a fraction of one thread's
       share of the work. */
    blocks = 4 * threads;
    if ((ulong) blocks > M)
        blocks = (int) M;
    block_first = flint_malloc((ulong) blocks * sizeof(ulong));
    block_last = flint_malloc((ulong) blocks * sizeof(ulong));
    block_weight = flint_malloc((ulong) blocks * sizeof(double));
    block_order = flint_malloc((ulong) blocks * sizeof(int));
    block_seconds = flint_malloc((ulong) blocks * sizeof(double));
    {
        double total = 0.0, accumulated = 0.0, closed = 0.0;
        int index = 0;

        for (ulong n = 1; n <= M; ++n)
            total += term_cost(n, n_direct, N, precision, log2_N);
        block_first[0] = 1;
        for (ulong n = 1; n <= M && index < blocks - 1; ++n)
        {
            accumulated += term_cost(n, n_direct, N, precision, log2_N);
            if (n >= block_first[index]
                && (accumulated * blocks >= total * (index + 1)
                    || M - n == (ulong) (blocks - 1 - index)))
            {
                block_last[index] = n;
                block_weight[index] = accumulated - closed;
                closed = accumulated;
                ++index;
                block_first[index] = n + 1;
            }
        }
        block_last[blocks - 1] = M;
        block_weight[blocks - 1] = total - closed;
    }
    for (int i = 0; i < blocks; ++i)
        block_order[i] = i;
    for (int i = 1; i < blocks; ++i)
    {
        int moved = block_order[i];
        int j = i;
        while (j > 0 && block_weight[block_order[j - 1]]
            < block_weight[moved])
        {
            block_order[j] = block_order[j - 1];
            --j;
        }
        block_order[j] = moved;
    }

#pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        arb_t block_sum;

        arb_init(block_sum);
        arb_zero(partials + tid);
        /* Each block accumulates alone: the range functions round their
           accumulator at the precision matched to the block's own term
           magnitudes, which would destroy an accumulator already holding
           larger earlier blocks. */
#pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < blocks; ++i)
        {
            int block = block_order[i];
            ulong first = block_first[block];
            ulong last = block_last[block];
            double block_start = omp_get_wtime();

            arb_zero(block_sum);
            if (last >= n_direct)
                scp_direct_range(block_sum,
                    first > n_direct ? first : n_direct, last,
                    N, precision, log2_N);
            if (first < n_direct)
                scp_bernoulli_range(block_sum, first,
                    last < n_direct - 1 ? last : n_direct - 1,
                    N, precision, log2_N);
            arb_add(partials + tid, partials + tid, block_sum, precision);
            block_seconds[block] = omp_get_wtime() - block_start;
        }
        arb_clear(block_sum);
    }

    if (getenv("KHINCHIN_VERBOSE") != NULL)
        for (int i = 0; i < blocks; ++i)
            fprintf(stderr, "block %d: n in [%lu, %lu], %.3f s\n",
                i, block_first[i], block_last[i], block_seconds[i]);
    flint_free(block_order);
    flint_free(block_weight);
    flint_free(block_seconds);
    flint_free(block_last);
    flint_free(block_first);

    for (int i = 0; i < threads; ++i)
        arb_add(result, result, partials + i, precision);

    arb_set_ui(t, N);
    arb_pow_ui(t, t, 2 * M, MAG_BITS);
    arb_inv(t, t, MAG_BITS);
    arb_add_error(result, t);

    arb_log_ui(t, 2, precision);
    arb_div(result, result, t, precision);
    arb_exp(result, result, precision);

    _arb_vec_clear(partials, threads);
    arb_clear(t);
}

static int rounded_scaled_integer(fmpz_t rounded, const arb_t value,
    ulong digits, slong precision)
{
    fmpz_t power10;
    arb_t scaled, half, floored;
    int unique;

    fmpz_init(power10);
    arb_init(scaled);
    arb_init(half);
    arb_init(floored);

    fmpz_ui_pow_ui(power10, 10, digits);
    arb_set_fmpz(scaled, power10);
    arb_mul(scaled, scaled, value, precision);

    arb_one(half);
    arb_mul_2exp_si(half, half, -1);
    arb_add(scaled, scaled, half, precision);
    arb_floor(floored, scaled, precision);
    unique = arb_get_unique_fmpz(rounded, floored);

    arb_clear(floored);
    arb_clear(half);
    arb_clear(scaled);
    fmpz_clear(power10);
    return unique;
}

static int write_decimal(FILE *output, const fmpz_t scaled, ulong digits)
{
    char *decimal = fmpz_get_str(NULL, 10, scaled);
    size_t length;
    int ok = 1;

    if (decimal == NULL)
        return 0;

    length = strlen(decimal);
    if (digits == 0)
    {
        ok = fprintf(output, "%s\n", decimal) >= 0;
    }
    else if (length > digits)
    {
        size_t integer_digits = length - (size_t) digits;
        ok = fwrite(decimal, 1, integer_digits, output) == integer_digits;
        ok = ok && fputc('.', output) != EOF;
        ok = ok && fwrite(decimal + integer_digits, 1, (size_t) digits, output)
            == (size_t) digits;
        ok = ok && fputc('\n', output) != EOF;
    }
    else
    {
        ok = fputs("0.", output) >= 0;
        for (size_t i = length; ok && i < digits; ++i)
            ok = fputc('0', output) != EOF;
        ok = ok && fwrite(decimal, 1, length, output) == length;
        ok = ok && fputc('\n', output) != EOF;
    }

    flint_free(decimal);
    return ok;
}

int main(int argc, char **argv)
{
    int benchmark_only = 0;
    int force = 0;
    const char *digits_text;
    const char *output_path = NULL;
    ulong digits;
    slong precision;
    arb_t khinchin;
    fmpz_t rounded;
    struct timespec total_start, compute_start, conversion_start;
    double compute_seconds, conversion_seconds = 0.0, total_seconds;
    int threads;
    const char *backend_name;
    long double estimated_bytes, physical_bytes;

    if (argc == 3 && strcmp(argv[1], "--benchmark") == 0)
    {
        benchmark_only = 1;
        digits_text = argv[2];
    }
    else if (argc == 4 && strcmp(argv[1], "--force") == 0)
    {
        force = 1;
        digits_text = argv[2];
        output_path = argv[3];
    }
    else if (argc == 3)
    {
        digits_text = argv[1];
        output_path = argv[2];
    }
    else
    {
        usage(stderr, argv[0]);
        return EXIT_FAILURE;
    }

    if (!parse_digits(digits_text, &digits))
    {
        fprintf(stderr, "Invalid DIGITS value: %s\n", digits_text);
        return EXIT_FAILURE;
    }

    precision = working_precision(digits);
    threads = omp_get_max_threads();
    backend_name = getenv("KHINCHIN_BACKEND");
    if (backend_name == NULL || backend_name[0] == '\0')
        backend_name = "scp";
    if (strcmp(backend_name, "scp") != 0 && strcmp(backend_name, "arb") != 0)
    {
        fprintf(stderr, "KHINCHIN_BACKEND must be 'scp' or 'arb'.\n");
        return EXIT_FAILURE;
    }
    estimated_bytes = estimated_peak_bytes(digits, precision, threads);
    physical_bytes = physical_memory_bytes();
    if (!force && physical_bytes > 0.0L
        && estimated_bytes > physical_bytes * 0.80L)
    {
        fprintf(stderr,
            "Refusing job: estimated peak memory %.2Lf GiB exceeds 80%% of "
            "physical RAM (%.2Lf GiB).\nUse --force to override.\n",
            estimated_bytes / 1073741824.0L,
            physical_bytes / 1073741824.0L);
        return EXIT_FAILURE;
    }
    arb_init(khinchin);
    fmpz_init(rounded);

    clock_gettime(CLOCK_MONOTONIC, &total_start);
    compute_start = total_start;
    if (strcmp(backend_name, "scp") == 0)
        khinchin_scp(khinchin, precision, threads);
    else
        khinchin_parallel(khinchin, precision, threads);
    compute_seconds = seconds_since(&compute_start);

    if (!benchmark_only)
    {
        FILE *output;

        clock_gettime(CLOCK_MONOTONIC, &conversion_start);
        if (!rounded_scaled_integer(rounded, khinchin, digits, precision))
        {
            fprintf(stderr, "Internal error: unable to certify decimal rounding.\n");
            if (getenv("KHINCHIN_VERBOSE") != NULL)
            {
                fprintf(stderr, "computed ball: ");
                arb_printd(khinchin, 30);
                fprintf(stderr, "\n");
            }
            fmpz_clear(rounded);
            arb_clear(khinchin);
            flint_cleanup();
            return EXIT_FAILURE;
        }

        output = fopen(output_path, "wb");
        if (output == NULL)
        {
            fprintf(stderr, "Cannot open %s: %s\n", output_path, strerror(errno));
            fmpz_clear(rounded);
            arb_clear(khinchin);
            flint_cleanup();
            return EXIT_FAILURE;
        }

        if (!write_decimal(output, rounded, digits) || fclose(output) != 0)
        {
            fprintf(stderr, "Failed while writing %s: %s\n", output_path,
                strerror(errno));
            fmpz_clear(rounded);
            arb_clear(khinchin);
            flint_cleanup();
            return EXIT_FAILURE;
        }
        conversion_seconds = seconds_since(&conversion_start);
    }

    total_seconds = seconds_since(&total_start);
    fprintf(stderr,
        "digits_after_decimal=%lu\n"
        "working_precision_bits=%ld\n"
        "backend=%s\n"
        "threads=%d\n"
        "estimated_peak_memory_mib=%.2Lf\n"
        "compute_seconds=%.6f\n"
        "conversion_and_write_seconds=%.6f\n"
        "total_seconds=%.6f\n",
        digits, precision, backend_name, threads, estimated_bytes / 1048576.0L,
        compute_seconds, conversion_seconds, total_seconds);

    fmpz_clear(rounded);
    arb_clear(khinchin);
    flint_cleanup();
    return EXIT_SUCCESS;
}
