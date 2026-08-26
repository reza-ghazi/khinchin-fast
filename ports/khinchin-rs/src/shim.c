/* Minimal C shim around FLINT's reverse Bernoulli iterator.

   bernoulli_rev_t's struct layout (it embeds arb_t and fmpz arrays) is
   not something Rust FFI can declare portably, so the shim owns the
   iterator behind an opaque pointer and hands each B_n out as GMP mpz
   numerator/denominator pairs, which rug::Integer can receive directly. */

#include <stdlib.h>
#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/bernoulli.h>

typedef struct
{
    bernoulli_rev_t iterator;
    fmpz_t numerator;
    fmpz_t denominator;
} k_rev;

void *k_rev_alloc(ulong n)
{
    k_rev *r = malloc(sizeof(k_rev));
    bernoulli_rev_init(r->iterator, n);
    fmpz_init(r->numerator);
    fmpz_init(r->denominator);
    return r;
}

/* Writes the next B_n (descending even n) into num/den. */
void k_rev_next(mpz_ptr num, mpz_ptr den, void *p)
{
    k_rev *r = p;
    bernoulli_rev_next(r->numerator, r->denominator, r->iterator);
    fmpz_get_mpz(num, r->numerator);
    fmpz_get_mpz(den, r->denominator);
}

void k_rev_free(void *p)
{
    k_rev *r = p;
    bernoulli_rev_clear(r->iterator);
    fmpz_clear(r->denominator);
    fmpz_clear(r->numerator);
    free(r);
}
