#!/usr/bin/perl
# Khinchin's constant (OEIS A002210) via the accelerated zeta series.
#
# Same mathematics as ../khinchin_fast.c:
#
#   ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
#                  + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
#
#   h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
#
# Runs on Math::BigInt in fixed point (scale 2^wp) — the "native bignum"
# data point for Perl, core modules only (the GMP backend is picked up
# automatically when installed, FastCalc or pure-Perl Calc otherwise;
# expect pure-Perl backends to be slow).  The port carries its own
# fixed-point kit: pi by Machin's formula, logarithms via the atanh(1/q)
# series, exp by argument-halving plus Taylor; even zeta values come
# from the positive-term recurrence (n + 1/2) zeta(2n) = sum_j zeta(2j)
# zeta(2n-2j) - but only below n_direct, the C program's two-region
# split: above it each term is a literal tail sum with no zeta values,
# shortening the O(M^2) recurrence to O(n_direct^2).  Guard bits absorb
# the drift; nothing here is interval-certified.
#
# Usage:  perl khinchin.pl DIGITS [OUTPUT_FILE]        (default 100)
# With OUTPUT_FILE the decimal value plus a newline is written there
# (the same contract as ../khinchin_fast.c); otherwise stdout.

use strict;
use warnings;
use Math::BigInt try => 'GMP,FastCalc,Calc';
use POSIX qw(ceil floor);

my ($wp, $one);

# atanh(1/q) = sum_{j>=0} 1/((2j+1) q^(2j+1)); atan(1/q) alternates.
sub arc_recip {
    my ($q, $alternating) = @_;
    my $q2    = $q * $q;
    my $power = $one->copy->bdiv($q);
    my $s     = $power->copy;
    for (my $j = 1 ; $power->is_positive ; $j++) {
        $power->bdiv($q2);
        my $term = $power->copy->bdiv(2 * $j + 1);
        if ($alternating && $j % 2 == 1) { $s->bsub($term) }
        else                             { $s->badd($term) }
    }
    return $s;
}

# exp(y) for 0 <= y < 1: halve 16 times, Taylor, square back.
sub exp_fixed {
    my ($y) = @_;
    my $arg  = $y->copy->brsft(16);
    my $acc  = $one->copy;
    my $term = $one->copy;
    for (my $j = 1 ; $term->is_positive ; $j++) {
        $term->bmul($arg)->brsft($wp)->bdiv($j);
        $acc->badd($term);
    }
    $acc->bmul($acc)->brsft($wp) for 1 .. 16;
    return $acc;
}

sub khinchin_fixed {
    my ($digits) = @_;
    my $prec = ceil(($digits + 2) * log(10) / log(2)) + 32;
    $wp  = $prec + ceil(sqrt($prec)) + 64;
    $one = Math::BigInt->new(1)->blsft($wp);
    my $N = 3 > floor($wp**0.35) ? 3 : floor($wp**0.35);
    my $M = ceil($wp * log(2) / (2 * log($N))) + 1;
    my $nd = ceil($wp / (2 * (log($N) / log(2) + 3)));
    $nd = $M + 1 if $nd > $M + 1;

    # -ln((k-1)/k) ln((k+1)/k) = 4 atanh(1/(2k-1)) atanh(1/(2k+1)).
    # NOTE: never call b* mutators inside argument lists - several of
    # them (bdiv, brsft) return (quotient, remainder) in list context,
    # and the stray remainder is then taken as an accuracy parameter
    # that silently rounds the invocant and poisons later results.
    my $s = Math::BigInt->new(0);
    for my $k (2 .. $N - 1) {
        my $p = arc_recip(2 * $k - 1, 0)->bmul(arc_recip(2 * $k + 1, 0));
        $p->brsft($wp)->bmul(4);
        $s->badd($p);
    }

    # pi = 16 atan(1/5) - 4 atan(1/239); zeta(2) = pi^2/6 seeds the
    # positive-term recurrence (n + 1/2) z_n = sum z_j z_{n-j}.
    my $pi = arc_recip(5, 1)->bmul(16)->bsub(arc_recip(239, 1)->bmul(4));
    my @z;
    $z[1] = $pi->copy->bmul($pi)->brsft($wp)->bdiv(6);
    for my $n (2 .. $nd - 1) {
        my $half = int(($n - 1) / 2);
        my $raw  = Math::BigInt->new(0);
        for my $j (1 .. $half) {
            my $prod = $z[$j]->copy->bmul($z[ $n - $j ]);
            $raw->badd($prod);
        }
        my $acc = $raw->brsft($wp)->bmul(2);
        if ($n % 2 == 0) {
            my $sq = $z[ $n / 2 ]->copy->bmul($z[ $n / 2 ]);
            $sq->brsft($wp);
            $acc->badd($sq);
        }
        $z[$n] = $acc->bmul(2)->bdiv(2 * $n + 1);
    }

    # Main accelerated loop, ascending, incremental power table.
    my @ks     = (2 .. $N - 1);
    # scalar() matters: map blocks run in list context, where bdiv
    # returns (quotient, remainder) and would double the array.
    my @powers = map { scalar $one->copy->bdiv($_ * $_) } @ks;
    my $h      = $one->copy;
    for my $n (1 .. $nd - 1) {
        my $tail = $z[$n]->copy->bsub($one);
        $tail->bsub($_) for @powers;
        $tail->bmul($h)->bdiv($n)->brsft($wp);
        $s->badd($tail);
        my $plus  = $one->copy->bdiv(2 * $n + 1);
        my $minus = $one->copy->bdiv(2 * $n);
        $h->badd($plus)->bsub($minus);
        $powers[$_]->bdiv($ks[$_] * $ks[$_]) for 0 .. $#ks;
    }

    # Direct region (the C program's split): literal tail sums over
    # k in [N, K], no zeta values; h continues ascending, power table
    # re-seeded per block of 32 terms.
    my $n = $nd;
    while ($n <= $M) {
        my $lastn = $n + 31 < $M ? $n + 31 : $M;
        my $K = floor(2**(($wp + 32) / (2 * $n))) + 1;
        $K = 16 * $N + 64 if $K > 16 * $N + 64;
        $K = $N if $K < $N;
        my @dks = ($N .. $K);
        my @pw;
        for my $k (@dks) {
            my $e = Math::BigInt->new($k)->bpow(2 * $n);
            my $p = $one->copy->bdiv($e);
            push @pw, $p;
        }
        for my $mm ($n .. $lastn) {
            my $tail = Math::BigInt->new(0);
            $tail->badd($_) for @pw;
            $tail->bmul($h)->bdiv($mm)->brsft($wp);
            $s->badd($tail);
            my $plus  = $one->copy->bdiv(2 * $mm + 1);
            my $minus = $one->copy->bdiv(2 * $mm);
            $h->badd($plus)->bsub($minus);
            $pw[$_]->bdiv($dks[$_] * $dks[$_]) for 0 .. $#dks;
        }
        $n = $lastn + 1;
    }

    my $ln2 = arc_recip(3, 0)->bmul(2);
    my $y = $s->blsft($wp)->bdiv($ln2);
    return exp_fixed($y);
}

my $digits = @ARGV ? int($ARGV[0]) : 100;
my $value  = khinchin_fixed($digits);
my $half_ulp = Math::BigInt->new(1)->blsft($wp - 1);
my $scaled = $value->bmul(Math::BigInt->new(10)->bpow($digits));
$scaled->badd($half_ulp)->brsft($wp);
my $text   = $scaled->bstr;
my $result = substr($text, 0, 1) . '.' . substr($text, 1);
if (@ARGV > 1) {
    open my $fh, '>', $ARGV[1] or die "Cannot open $ARGV[1]: $!";
    print {$fh} "$result\n";
    close $fh;
}
else {
    print "$result\n";
}
