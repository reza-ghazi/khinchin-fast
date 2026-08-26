! Khinchin's constant (OEIS A002210) via the accelerated zeta series.
!
! Same mathematics as ../khinchin_fast.c:
!
!   ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
!                  + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
!
!   h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j.
!
! Standard Fortran has no arbitrary-precision arithmetic, so this port
! binds FLINT/Arb directly through ISO_C_BINDING - the same library the
! C program uses - working only through Arb's pointer-based API
! (_arb_vec_init), so no C struct layouts need to be declared.  Each
! zeta(2n) comes from arb_zeta_ui at full precision, the strategy of the
! C program's KHINCHIN_BACKEND=arb cross-check backend, with the n-range
! split across OpenMP threads in equal-count blocks.  GNU MPFR (also
! bound below) formats the final decimal output.
!
! Build and run:
!   gfortran -O2 -fopenmp -o khinchin-f90 khinchin.f90 -lflint -lmpfr -lgmp
!   ./khinchin-f90 DIGITS [OUTPUT_FILE]    (default 100)
!
! With OUTPUT_FILE the decimal value plus a newline is written there
! (the same contract as ../khinchin_fast.c); otherwise stdout.

module mpfr_binding
    use iso_c_binding
    implicit none

    integer(c_int), parameter :: RNDN = 0_c_int   ! round to nearest

    ! Layout of __mpfr_struct on LP64 (mpfr.h); mpfr_t is one of these.
    type, bind(c) :: mpfr_t
        integer(c_long) :: prec
        integer(c_int)  :: sign
        integer(c_long) :: exp
        type(c_ptr)     :: d
    end type

    interface
        subroutine mpfr_init2(x, prec) bind(c)
            import :: mpfr_t, c_long
            type(mpfr_t), intent(inout) :: x
            integer(c_long), value :: prec
        end subroutine
        subroutine mpfr_clear(x) bind(c)
            import :: mpfr_t
            type(mpfr_t), intent(inout) :: x
        end subroutine
        type(c_ptr) function mpfr_get_str(buf, expo, base, n, x, rnd) bind(c)
            import :: mpfr_t, c_int, c_long, c_ptr, c_size_t
            type(c_ptr), value :: buf
            integer(c_long), intent(out) :: expo
            integer(c_int), value :: base
            integer(c_size_t), value :: n
            type(mpfr_t), intent(in) :: x
            integer(c_int), value :: rnd
        end function
        subroutine mpfr_free_str(p) bind(c)
            import :: c_ptr
            type(c_ptr), value :: p
        end subroutine
    end interface
end module mpfr_binding

module arb_binding
    use iso_c_binding
    implicit none

    ! Arb values are handled purely as opaque pointers from
    ! _arb_vec_init, so no struct layout is required.
    interface
        type(c_ptr) function arb_vec_init(n) bind(c, name='_arb_vec_init')
            import :: c_ptr, c_long
            integer(c_long), value :: n
        end function
        subroutine arb_vec_clear(v, n) bind(c, name='_arb_vec_clear')
            import :: c_ptr, c_long
            type(c_ptr), value :: v
            integer(c_long), value :: n
        end subroutine
        subroutine arb_zero(x) bind(c)
            import :: c_ptr
            type(c_ptr), value :: x
        end subroutine
        subroutine arb_one(x) bind(c)
            import :: c_ptr
            type(c_ptr), value :: x
        end subroutine
        subroutine arb_set_ui(x, u) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: x
            integer(c_long), value :: u
        end subroutine
        subroutine arb_add(r, a, b, prec) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: r, a, b
            integer(c_long), value :: prec
        end subroutine
        subroutine arb_sub(r, a, b, prec) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: r, a, b
            integer(c_long), value :: prec
        end subroutine
        subroutine arb_mul(r, a, b, prec) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: r, a, b
            integer(c_long), value :: prec
        end subroutine
        subroutine arb_div(r, a, b, prec) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: r, a, b
            integer(c_long), value :: prec
        end subroutine
        subroutine arb_sub_ui(r, a, u, prec) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: r, a
            integer(c_long), value :: u, prec
        end subroutine
        subroutine arb_div_ui(r, a, u, prec) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: r, a
            integer(c_long), value :: u, prec
        end subroutine
        subroutine arb_mul_ui(r, a, u, prec) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: r, a
            integer(c_long), value :: u, prec
        end subroutine
        subroutine arb_inv(r, a, prec) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: r, a
            integer(c_long), value :: prec
        end subroutine
        subroutine arb_pow_ui(r, a, e, prec) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: r, a
            integer(c_long), value :: e, prec
        end subroutine
        subroutine arb_log(r, a, prec) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: r, a
            integer(c_long), value :: prec
        end subroutine
        subroutine arb_log_ui(r, u, prec) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: r
            integer(c_long), value :: u, prec
        end subroutine
        subroutine arb_exp(r, a, prec) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: r, a
            integer(c_long), value :: prec
        end subroutine
        subroutine arb_zeta_ui(r, u, prec) bind(c)
            import :: c_ptr, c_long
            type(c_ptr), value :: r
            integer(c_long), value :: u, prec
        end subroutine
        subroutine arb_get_interval_mpfr(a, b, x) bind(c)
            import :: c_ptr
            type(c_ptr), value :: a, b, x
        end subroutine
    end interface
end module arb_binding

program khinchin_prog
    use iso_c_binding
    use mpfr_binding
    use arb_binding
    use omp_lib
    implicit none

    integer :: digits, prec, wp, bign, m, k, i, tid, nthreads
    integer(c_long) :: wpl, expo
    type(c_ptr) :: s, t1, t2
    type(c_ptr), allocatable :: partials(:)
    type(mpfr_t), target :: out_lo, out_hi
    character(len=32) :: arg
    character(len=1024) :: outpath
    integer :: outunit
    type(c_ptr) :: cstr
    character(kind=c_char), pointer :: chars(:)
    character(len=:), allocatable :: out

    digits = 100
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) digits
    end if
    outpath = ''
    if (command_argument_count() >= 2) call get_command_argument(2, outpath)

    prec = ceiling((digits + 2) * log(10.0d0) / log(2.0d0)) + 32
    wp = prec + ceiling(sqrt(real(prec, 8))) + 64
    wpl = int(wp, c_long)
    bign = max(3, int(real(wp, 8)**0.35d0))
    m = ceiling(wp * log(2.0d0) / (2 * log(real(bign, 8)))) + 1

    s = arb_vec_init(1_c_long)
    t1 = arb_vec_init(1_c_long)
    t2 = arb_vec_init(1_c_long)
    call arb_zero(s)

    ! Finite logarithmic correction.
    do k = 2, bign - 1
        call arb_set_ui(t1, int(k - 1, c_long))
        call arb_div_ui(t1, t1, int(k, c_long), wpl)
        call arb_log(t1, t1, wpl)
        call arb_set_ui(t2, int(k + 1, c_long))
        call arb_div_ui(t2, t2, int(k, c_long), wpl)
        call arb_log(t2, t2, wpl)
        call arb_mul(t1, t1, t2, wpl)
        call arb_sub(s, s, t1, wpl)
    end do

    ! Zeta range in equal-count blocks, one per thread, ascending n with
    ! incremental power table and harmonic weight, as in the C program's
    ! KHINCHIN_BACKEND=arb backend.
    nthreads = min(omp_get_max_threads(), m)
    allocate(partials(nthreads))
    do i = 1, nthreads
        partials(i) = arb_vec_init(1_c_long)
    end do
    !$omp parallel num_threads(nthreads) private(tid)
    tid = omp_get_thread_num() + 1
    call zeta_block(partials(tid), &
        1 + (m * (tid - 1)) / nthreads, (m * tid) / nthreads)
    !$omp end parallel
    do i = 1, nthreads
        call arb_add(s, s, partials(i), wpl)
        call arb_vec_clear(partials(i), 1_c_long)
    end do

    call arb_log_ui(t1, 2_c_long, wpl)
    call arb_div(s, s, t1, wpl)
    call arb_exp(s, s, wpl)

    ! Midpoint out through MPFR, then decimal formatting as before.
    call mpfr_init2(out_lo, wpl)
    call mpfr_init2(out_hi, wpl)
    call arb_get_interval_mpfr(c_loc(out_lo), c_loc(out_hi), s)
    cstr = mpfr_get_str(c_null_ptr, expo, 10_c_int, &
        int(digits + 1, c_size_t), out_lo, RNDN)
    call c_f_pointer(cstr, chars, [digits + 2])
    allocate(character(len=digits + 2) :: out)
    out(1:1) = chars(1)
    out(2:2) = '.'
    do i = 2, digits + 1
        out(i + 1:i + 1) = chars(i)
    end do
    if (len_trim(outpath) > 0) then
        open(newunit=outunit, file=trim(outpath), status='replace', &
            action='write')
        write(outunit, '(A)') out
        close(outunit)
    else
        print '(A)', out
    end if
    call mpfr_free_str(cstr)
    call mpfr_clear(out_hi)
    call mpfr_clear(out_lo)
    call arb_vec_clear(t2, 1_c_long)
    call arb_vec_clear(t1, 1_c_long)
    call arb_vec_clear(s, 1_c_long)

contains

    ! Accelerated terms for n in [first, last]: zeta(2n) via arb_zeta_ui,
    ! shared incremental power table, ascending harmonic weight.
    subroutine zeta_block(sb, first, last)
        type(c_ptr), intent(in) :: sb
        integer, intent(in) :: first, last
        type(c_ptr) :: zt, h, term, recip, powers
        type(c_ptr), allocatable :: pw(:)
        integer :: n, k, j

        zt = arb_vec_init(1_c_long)
        h = arb_vec_init(1_c_long)
        term = arb_vec_init(1_c_long)
        recip = arb_vec_init(1_c_long)
        call arb_zero(sb)

        ! h(first) = 1 - sum_{n<first} 1/(2n(2n+1)).
        call arb_one(h)
        do n = 1, first - 1
            call arb_set_ui(recip, int(2 * n, c_long))
            call arb_mul_ui(recip, recip, int(2 * n + 1, c_long), wpl)
            call arb_inv(recip, recip, wpl)
            call arb_sub(h, h, recip, wpl)
        end do

        ! Power table entry k starts at k^(-2(first-1)); the loop
        ! advances it once per n.
        allocate(pw(2:bign - 1))
        do k = 2, bign - 1
            pw(k) = arb_vec_init(1_c_long)
            call arb_set_ui(pw(k), int(k, c_long))
            call arb_pow_ui(pw(k), pw(k), int(2 * (first - 1), c_long), wpl)
            call arb_inv(pw(k), pw(k), wpl)
        end do

        do n = first, last
            call arb_zeta_ui(zt, int(2 * n, c_long), wpl)
            call arb_sub_ui(zt, zt, 1_c_long, wpl)
            do k = 2, bign - 1
                call arb_div_ui(pw(k), pw(k), int(k * k, c_long), wpl)
                call arb_sub(zt, zt, pw(k), wpl)
            end do
            call arb_div_ui(term, zt, int(n, c_long), wpl)
            call arb_mul(term, term, h, wpl)
            call arb_add(sb, sb, term, wpl)
            call arb_set_ui(recip, int(2 * n, c_long))
            call arb_mul_ui(recip, recip, int(2 * n + 1, c_long), wpl)
            call arb_inv(recip, recip, wpl)
            call arb_sub(h, h, recip, wpl)
        end do

        do k = 2, bign - 1
            call arb_vec_clear(pw(k), 1_c_long)
        end do
        call arb_vec_clear(recip, 1_c_long)
        call arb_vec_clear(term, 1_c_long)
        call arb_vec_clear(h, 1_c_long)
        call arb_vec_clear(zt, 1_c_long)
    end subroutine zeta_block

end program khinchin_prog
