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
! runs in IEEE quadruple precision (REAL128, ~33 significant digits) and
! prints 30 digits after the decimal point, all correct against this
! repository's reference files.  It demonstrates the algorithm, not the
! record chase - for thousands of digits use ../khinchin_fast.c.
!
! Build and run:  gfortran -O2 -o khinchin-f90 khinchin.f90 && ./khinchin-f90

program khinchin
    use iso_fortran_env, only: real128
    implicit none
    integer, parameter :: qp = real128
    integer, parameter :: bits = 113          ! REAL128 mantissa
    integer :: bign, m, n, j, k
    real(qp) :: s, h, tail
    real(qp), allocatable :: z(:), powers(:)
    character(len=64) :: text

    bign = max(3, int(real(bits, qp)**0.35_qp))
    m = ceiling(bits * log(2.0_qp) / (2 * log(real(bign, qp)))) + 1

    ! Finite logarithmic correction.
    s = 0
    do k = 2, bign - 1
        s = s - log((k - 1.0_qp) / k) * log((k + 1.0_qp) / k)
    end do

    ! zeta(2), zeta(4), ..., zeta(2M) by the positive recurrence.
    allocate(z(m))
    z(1) = (4 * atan(1.0_qp))**2 / 6
    do n = 2, m
        tail = 0
        do j = 1, n - 1
            tail = tail + z(j) * z(n - j)
        end do
        z(n) = 2 * tail / (2 * n + 1)
    end do

    allocate(powers(2:bign - 1))
    do k = 2, bign - 1
        powers(k) = 1 / real(k, qp)**2
    end do

    h = 1
    do n = 1, m
        tail = z(n) - 1
        do k = 2, bign - 1
            tail = tail - powers(k)
        end do
        s = s + tail * h / n
        h = h - 1 / (2 * n * (2.0_qp * n + 1))
        do k = 2, bign - 1
            powers(k) = powers(k) / real(k * k, qp)
        end do
    end do

    s = exp(s / log(2.0_qp))
    write(text, '(F0.30)') s
    write(*, '(A)') trim(text)
end program khinchin
