! Khinchin's constant (OEIS A002210) via the accelerated zeta series.
!
! Same mathematics as ../khinchin_fast.c:
!
!   ln(2) ln(K0) = -sum_{k=2}^{N-1} ln((k-1)/k) ln((k+1)/k)
!                  + sum_{n>=1} (zeta(2n) - 1 - sum_{k=2}^{N-1} k^(-2n)) / n * h(n),
!
!   h(n) = sum_{j=1}^{2n-1} (-1)^(j+1)/j,
!
! with the even zeta values from the classical positive-term recurrence
!
!   (n + 1/2) zeta(2n) = sum_{j=1}^{n-1} zeta(2j) zeta(2n-2j),  zeta(2) = pi^2/6.
!
! Standard Fortran has no arbitrary-precision arithmetic, so this port
! binds GNU MPFR directly through ISO_C_BINDING - the same engine that
! backs Julia's BigFloat and Rust's rug in the sibling ports, so all of
! them produce identical digits.  OpenMP parallelises the recurrence
! convolutions (halved via their j <-> n-j symmetry) and the main loop's
! equal-count blocks, each seeded with h(first) summed directly and a
! power table started at k^(-2 first).  Guard bits absorb the recurrence
! drift; unlike the C program, nothing here is interval-certified.
!
! Build and run:
!   gfortran -O2 -fopenmp -o khinchin-f90 khinchin.f90 -lmpfr -lgmp
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
        integer(c_int) function mpfr_set_ui(x, u, rnd) bind(c)
            import :: mpfr_t, c_int, c_long
            type(mpfr_t), intent(inout) :: x
            integer(c_long), value :: u
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_add(r, a, b, rnd) bind(c)
            import :: mpfr_t, c_int
            type(mpfr_t), intent(inout) :: r
            type(mpfr_t), intent(in) :: a, b
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_sub(r, a, b, rnd) bind(c)
            import :: mpfr_t, c_int
            type(mpfr_t), intent(inout) :: r
            type(mpfr_t), intent(in) :: a, b
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_mul(r, a, b, rnd) bind(c)
            import :: mpfr_t, c_int
            type(mpfr_t), intent(inout) :: r
            type(mpfr_t), intent(in) :: a, b
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_div(r, a, b, rnd) bind(c)
            import :: mpfr_t, c_int
            type(mpfr_t), intent(inout) :: r
            type(mpfr_t), intent(in) :: a, b
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_sqr(r, a, rnd) bind(c)
            import :: mpfr_t, c_int
            type(mpfr_t), intent(inout) :: r
            type(mpfr_t), intent(in) :: a
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_mul_ui(r, a, u, rnd) bind(c)
            import :: mpfr_t, c_int, c_long
            type(mpfr_t), intent(inout) :: r
            type(mpfr_t), intent(in) :: a
            integer(c_long), value :: u
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_div_ui(r, a, u, rnd) bind(c)
            import :: mpfr_t, c_int, c_long
            type(mpfr_t), intent(inout) :: r
            type(mpfr_t), intent(in) :: a
            integer(c_long), value :: u
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_sub_ui(r, a, u, rnd) bind(c)
            import :: mpfr_t, c_int, c_long
            type(mpfr_t), intent(inout) :: r
            type(mpfr_t), intent(in) :: a
            integer(c_long), value :: u
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_ui_div(r, u, a, rnd) bind(c)
            import :: mpfr_t, c_int, c_long
            type(mpfr_t), intent(inout) :: r
            integer(c_long), value :: u
            type(mpfr_t), intent(in) :: a
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_ui_pow_ui(r, b, e, rnd) bind(c)
            import :: mpfr_t, c_int, c_long
            type(mpfr_t), intent(inout) :: r
            integer(c_long), value :: b, e
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_log(r, a, rnd) bind(c)
            import :: mpfr_t, c_int
            type(mpfr_t), intent(inout) :: r
            type(mpfr_t), intent(in) :: a
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_exp(r, a, rnd) bind(c)
            import :: mpfr_t, c_int
            type(mpfr_t), intent(inout) :: r
            type(mpfr_t), intent(in) :: a
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_const_pi(r, rnd) bind(c)
            import :: mpfr_t, c_int
            type(mpfr_t), intent(inout) :: r
            integer(c_int), value :: rnd
        end function
        integer(c_int) function mpfr_const_log2(r, rnd) bind(c)
            import :: mpfr_t, c_int
            type(mpfr_t), intent(inout) :: r
            integer(c_int), value :: rnd
        end function
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

program khinchin_prog
    use iso_c_binding
    use mpfr_binding
    use omp_lib
    implicit none

    integer :: digits, prec, wp, bign, m, n, k, j, half, w, b, i, rc
    integer(c_long) :: wpl, expo
    type(mpfr_t) :: s, t1, t2, acc
    type(mpfr_t), allocatable :: z(:), partials(:)
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

    call mpfr_init2(s, wpl);  rc = mpfr_set_ui(s, 0_c_long, RNDN)
    call mpfr_init2(t1, wpl)
    call mpfr_init2(t2, wpl)
    call mpfr_init2(acc, wpl)

    ! Finite logarithmic correction.
    do k = 2, bign - 1
        rc = mpfr_set_ui(t1, int(k - 1, c_long), RNDN)
        rc = mpfr_div_ui(t1, t1, int(k, c_long), RNDN)
        rc = mpfr_log(t1, t1, RNDN)
        rc = mpfr_set_ui(t2, int(k + 1, c_long), RNDN)
        rc = mpfr_div_ui(t2, t2, int(k, c_long), RNDN)
        rc = mpfr_log(t2, t2, RNDN)
        rc = mpfr_mul(t1, t1, t2, RNDN)
        rc = mpfr_sub(s, s, t1, RNDN)
    end do

    ! zeta(2), zeta(4), ..., zeta(2M) by the positive recurrence, using
    ! the symmetry of the convolution: only j <= (n-1)/2 is summed,
    ! doubled, plus the middle square when n is even.
    allocate(z(m))
    do n = 1, m
        call mpfr_init2(z(n), wpl)
    end do
    rc = mpfr_const_pi(z(1), RNDN)
    rc = mpfr_sqr(z(1), z(1), RNDN)
    rc = mpfr_div_ui(z(1), z(1), 6_c_long, RNDN)
    do n = 2, m
        half = (n - 1) / 2
        rc = mpfr_set_ui(acc, 0_c_long, RNDN)
        if (half >= 64) then
            !$omp parallel private(j, rc) shared(z, acc, half, n, wpl)
            block
                type(mpfr_t) :: a, prod
                integer :: lo, hi, tid, nth
                call mpfr_init2(a, wpl);  rc = mpfr_set_ui(a, 0_c_long, RNDN)
                call mpfr_init2(prod, wpl)
                tid = omp_get_thread_num()
                nth = omp_get_num_threads()
                lo = 1 + (half * tid) / nth
                hi = (half * (tid + 1)) / nth
                do j = lo, hi
                    rc = mpfr_mul(prod, z(j), z(n - j), RNDN)
                    rc = mpfr_add(a, a, prod, RNDN)
                end do
                !$omp critical
                rc = mpfr_add(acc, acc, a, RNDN)
                !$omp end critical
                call mpfr_clear(prod)
                call mpfr_clear(a)
            end block
            !$omp end parallel
        else
            do j = 1, half
                rc = mpfr_mul(t1, z(j), z(n - j), RNDN)
                rc = mpfr_add(acc, acc, t1, RNDN)
            end do
        end if
        rc = mpfr_mul_ui(acc, acc, 2_c_long, RNDN)
        if (mod(n, 2) == 0) then
            rc = mpfr_sqr(t1, z(n / 2), RNDN)
            rc = mpfr_add(acc, acc, t1, RNDN)
        end if
        rc = mpfr_mul_ui(acc, acc, 2_c_long, RNDN)
        rc = mpfr_div_ui(z(n), acc, int(2 * n + 1, c_long), RNDN)
    end do

    ! Main loop in equal-count blocks, one thread per block.
    w = max(1, min(omp_get_max_threads(), m / 8))
    allocate(partials(w))
    do b = 1, w
        call mpfr_init2(partials(b), wpl)
    end do
    !$omp parallel do schedule(static, 1) private(b)
    do b = 1, w
        call block_sum(partials(b), 1 + (m * (b - 1)) / w, (m * b) / w)
    end do
    !$omp end parallel do
    do b = 1, w
        rc = mpfr_add(s, s, partials(b), RNDN)
        call mpfr_clear(partials(b))
    end do

    rc = mpfr_const_log2(t1, RNDN)
    rc = mpfr_div(s, s, t1, RNDN)
    rc = mpfr_exp(s, s, RNDN)

    ! digits + 1 significant digits, correctly rounded: "2" + digits.
    cstr = mpfr_get_str(c_null_ptr, expo, 10_c_int, &
        int(digits + 1, c_size_t), s, RNDN)
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

    do n = 1, m
        call mpfr_clear(z(n))
    end do
    call mpfr_clear(acc)
    call mpfr_clear(t2)
    call mpfr_clear(t1)
    call mpfr_clear(s)

contains

    ! Accelerated terms for n in [first, last], accumulated into sb.
    subroutine block_sum(sb, first, last)
        type(mpfr_t), intent(inout) :: sb
        integer, intent(in) :: first, last
        type(mpfr_t) :: h, tail, tmp
        type(mpfr_t), allocatable :: pw(:)
        integer :: n, k, j, rc

        call mpfr_init2(h, wpl);    rc = mpfr_set_ui(h, 0_c_long, RNDN)
        call mpfr_init2(tail, wpl)
        call mpfr_init2(tmp, wpl)
        rc = mpfr_set_ui(sb, 0_c_long, RNDN)

        ! h(first) = sum_{j=1}^{2 first - 1} (-1)^(j+1)/j.
        do j = 1, 2 * first - 1
            rc = mpfr_set_ui(tmp, 1_c_long, RNDN)
            rc = mpfr_div_ui(tmp, tmp, int(j, c_long), RNDN)
            if (mod(j, 2) == 1) then
                rc = mpfr_add(h, h, tmp, RNDN)
            else
                rc = mpfr_sub(h, h, tmp, RNDN)
            end if
        end do

        allocate(pw(2:bign - 1))
        do k = 2, bign - 1
            call mpfr_init2(pw(k), wpl)
            rc = mpfr_ui_pow_ui(pw(k), int(k, c_long), &
                int(2 * first, c_long), RNDN)
            rc = mpfr_ui_div(pw(k), 1_c_long, pw(k), RNDN)
        end do

        do n = first, last
            rc = mpfr_sub_ui(tail, z(n), 1_c_long, RNDN)
            do k = 2, bign - 1
                rc = mpfr_sub(tail, tail, pw(k), RNDN)
            end do
            rc = mpfr_mul(tail, tail, h, RNDN)
            rc = mpfr_div_ui(tail, tail, int(n, c_long), RNDN)
            rc = mpfr_add(sb, sb, tail, RNDN)
            rc = mpfr_set_ui(tmp, 1_c_long, RNDN)
            rc = mpfr_div_ui(tmp, tmp, int(2 * n, c_long) &
                * int(2 * n + 1, c_long), RNDN)
            rc = mpfr_sub(h, h, tmp, RNDN)
            do k = 2, bign - 1
                rc = mpfr_div_ui(pw(k), pw(k), int(k * k, c_long), RNDN)
            end do
        end do

        do k = 2, bign - 1
            call mpfr_clear(pw(k))
        end do
        call mpfr_clear(tmp)
        call mpfr_clear(tail)
        call mpfr_clear(h)
    end subroutine block_sum

end program khinchin_prog
