module T_tensor_damp_erf
    implicit none
    private

    public Tn_damp_erf
    public T0_damp_erf
    public T1_damp_erf
    public T2_damp_erf
    public T3_damp_erf
    public T4_damp_erf
    public T5_damp_erf
    public T6_damp_erf
contains
    pure subroutine Tn_damp_erf(n, x, y, z, a, T)
        integer, intent(in) :: n
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(in) :: a(:)
        real(8), intent(inout) :: T(size(x), (n + 1) * (n + 2) * (n + 3) / 6)
        select case (n)
        case (0)
            call T0_damp_erf(x, y, z, a, T(:, 1:1))
        case (1)
            call T0_damp_erf(x, y, z, a, T(:, 1:1))
            call T1_damp_erf(x, y, z, a, T(:, 2:4))
        case (2)
            call T0_damp_erf(x, y, z, a, T(:, 1:1))
            call T1_damp_erf(x, y, z, a, T(:, 2:4))
            call T2_damp_erf(x, y, z, a, T(:, 5:10))
        case (3)
            call T0_damp_erf(x, y, z, a, T(:, 1:1))
            call T1_damp_erf(x, y, z, a, T(:, 2:4))
            call T2_damp_erf(x, y, z, a, T(:, 5:10))
            call T3_damp_erf(x, y, z, a, T(:, 11:20))
        case (4)
            call T0_damp_erf(x, y, z, a, T(:, 1:1))
            call T1_damp_erf(x, y, z, a, T(:, 2:4))
            call T2_damp_erf(x, y, z, a, T(:, 5:10))
            call T3_damp_erf(x, y, z, a, T(:, 11:20))
            call T4_damp_erf(x, y, z, a, T(:, 21:35))
        case (5)
            call T0_damp_erf(x, y, z, a, T(:, 1:1))
            call T1_damp_erf(x, y, z, a, T(:, 2:4))
            call T2_damp_erf(x, y, z, a, T(:, 5:10))
            call T3_damp_erf(x, y, z, a, T(:, 11:20))
            call T4_damp_erf(x, y, z, a, T(:, 21:35))
            call T5_damp_erf(x, y, z, a, T(:, 36:56))
        case (6)
            call T0_damp_erf(x, y, z, a, T(:, 1:1))
            call T1_damp_erf(x, y, z, a, T(:, 2:4))
            call T2_damp_erf(x, y, z, a, T(:, 5:10))
            call T3_damp_erf(x, y, z, a, T(:, 11:20))
            call T4_damp_erf(x, y, z, a, T(:, 21:35))
            call T5_damp_erf(x, y, z, a, T(:, 36:56))
            call T6_damp_erf(x, y, z, a, T(:, 57:84))
        case default
            call T0_damp_erf(x, y, z, a, T(:, 1:1))
            call T1_damp_erf(x, y, z, a, T(:, 2:4))
            call T2_damp_erf(x, y, z, a, T(:, 5:10))
            call T3_damp_erf(x, y, z, a, T(:, 11:20))
            call T4_damp_erf(x, y, z, a, T(:, 21:35))
            call T5_damp_erf(x, y, z, a, T(:, 36:56))
            call T6_damp_erf(x, y, z, a, T(:, 57:84))
        end select
    end subroutine Tn_damp_erf
    pure subroutine T0_damp_erf(x, y, z, a, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(in) :: a(:)
        real(8), intent(inout) :: T(size(x), 1)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: ai
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            T(i, 1) = (xi**2 + yi**2 + zi**2)**(-0.5d0) * erf(sqrt(ai) * sqrt(xi**2 + yi**2 + zi**2))
        end do
    end subroutine T0_damp_erf
    pure subroutine T1_damp_erf(x, y, z, a, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(in) :: a(:)
        real(8), intent(inout) :: T(size(x), 4)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: ai
        real(8) :: x2
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x2 = 1.12837916709551d0 * sqrt(ai) * exp(-ai * (xi**2 + yi**2 + zi**2)) / (xi**2 + yi**2 + zi**2) - (xi**2 + yi**2 + zi**2)**(-1.5d0) * erf(sqrt(ai) * sqrt(xi**2 + yi**2 + zi**2))
            T(i, 1) = x2 * xi
            T(i, 2) = x2 * yi
            T(i, 3) = x2 * zi
        end do
    end subroutine T1_damp_erf
    pure subroutine T2_damp_erf(x, y, z, a, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(in) :: a(:)
        real(8), intent(inout) :: T(size(x), 10)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: ai
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x6
        real(8) :: x8
        real(8) :: x9
        real(8) :: x12
        real(8) :: x13
        real(8) :: x14
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 0.564189583547756d0 * exp(-ai * x3)
            x6 = sqrt(ai)
            x8 = 4.0d0 * x4 * x6 / x3**2
            x9 = 1d0 / x3
            x12 = x3**(-1.5d0) * erf(sqrt(x3) * x6)
            x13 = 2.0d0 * ai
            x14 = 2.0d0 * x4 * x6 * x9
            T(i, 1) = -x0 * x8 + x12 * (3.0d0 * x0 * x9 - 1.0d0) - x14 * (x0 * x13 + x0 * x9 - 1.0d0)
            T(i, 2) = -xi * yi * (4.0d0 * ai**1.5d0 * x4 * x9 - 3.0d0 * x3**(-2.5d0) * erf(sqrt(x3) * x6) + 6.0d0 * x4 * x6 / x3**2)
            T(i, 3) = -xi * zi * (4.0d0 * ai**1.5d0 * x4 * x9 - 3.0d0 * x3**(-2.5d0) * erf(sqrt(x3) * x6) + 6.0d0 * x4 * x6 / x3**2)
            T(i, 4) = -x1 * x8 + x12 * (3.0d0 * x1 * x9 - 1.0d0) - x14 * (x1 * x13 + x1 * x9 - 1.0d0)
            T(i, 5) = -yi * zi * (4.0d0 * ai**1.5d0 * x4 * x9 - 3.0d0 * x3**(-2.5d0) * erf(sqrt(x3) * x6) + 6.0d0 * x4 * x6 / x3**2)
            T(i, 6) = x12 * (3.0d0 * x2 * x9 - 1.0d0) - x14 * (x13 * x2 + x2 * x9 - 1.0d0) - x2 * x8
        end do
    end subroutine T2_damp_erf
    pure subroutine T3_damp_erf(x, y, z, a, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(in) :: a(:)
        real(8), intent(inout) :: T(size(x), 20)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: ai
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x6
        real(8) :: x8
        real(8) :: x9
        real(8) :: x12
        real(8) :: x14
        real(8) :: x15
        real(8) :: x16
        real(8) :: x18
        real(8) :: x19
        real(8) :: x20
        real(8) :: x21
        real(8) :: x22
        real(8) :: x23
        real(8) :: x33
        real(8) :: x35
        real(8) :: x36
        real(8) :: x39
        real(8) :: x40
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0 / x3
            x6 = sqrt(ai)
            x8 = 3.0d0 * x3**(-2.5d0) * erf(sqrt(x3) * x6)
            x9 = 3.0d0 * x4
            x12 = 0.564189583547756d0 * exp(-ai * x3)
            x14 = x12 * x6 / x3**2
            x15 = 6.0d0 * x14
            x16 = 2.0d0 * ai
            x18 = 3.0d0 / x3**2
            x19 = 4.0d0 * ai**2
            x20 = 4.0d0 * ai
            x21 = -6.0d0 * ai - x9
            x22 = x12 * x4
            x23 = 2.0d0 * x22 * x6
            x33 = 4.0d0 * ai**1.5d0 * x22
            x35 = 15.0d0 * x3**(-3.5d0) * erf(sqrt(x3) * x6)
            x36 = x1 * x4
            x39 = 20.0d0 * ai**1.5d0 * x12 / x3**2
            x40 = 30.0d0 * x12 * x6 / x3**3
            T(i, 1) = xi * (x15 * (x0 * x9 - 1.0d0) + x15 * (x0 * x16 + x0 * x4 - 1.0d0) + x23 * (x0 * x18 + x0 * x19 + x0 * x20 * x4 + x21) - x8 * (5.0d0 * x0 * x4 - 3.0d0))
            T(i, 2) = yi * (8.0d0 * ai**1.5d0 * x0 * x12 / x3**2 + 20.0d0 * x0 * x12 * x6 / x3**3 - 6.0d0 * x0 * x3**(-3.5d0) * erf(sqrt(x3) * x6) + 2.0d0 * x14 * (x0 * x9 - 1.0d0) + 4.0d0 * x14 * (x0 * x16 + x0 * x4 - 1.0d0) + x33 * (x0 * x16 + x0 * x4 - 1.0d0) - x8 * (x0 * x9 - 1.0d0))
            T(i, 3) = zi * (8.0d0 * ai**1.5d0 * x0 * x12 / x3**2 + 20.0d0 * x0 * x12 * x6 / x3**3 - 6.0d0 * x0 * x3**(-3.5d0) * erf(sqrt(x3) * x6) + 2.0d0 * x14 * (x0 * x9 - 1.0d0) + 4.0d0 * x14 * (x0 * x16 + x0 * x4 - 1.0d0) + x33 * (x0 * x16 + x0 * x4 - 1.0d0) - x8 * (x0 * x9 - 1.0d0))
            T(i, 4) = xi * (8.0d0 * ai**2.5d0 * x12 * x36 - x1 * x35 + x1 * x39 + x1 * x40 - x15 - x33 + x8)
            T(i, 5) = xi * yi * zi * (8.0d0 * ai**2.5d0 * x22 - x35 + x39 + x40)
            T(i, 6) = xi * (8.0d0 * ai**2.5d0 * x12 * x2 * x4 - x15 - x2 * x35 + x2 * x39 + x2 * x40 - x33 + x8)
            T(i, 7) = yi * (x15 * (x1 * x9 - 1.0d0) + x15 * (x1 * x16 + x36 - 1.0d0) + x23 * (x1 * x18 + x1 * x19 + x20 * x36 + x21) - x8 * (5.0d0 * x36 - 3.0d0))
            T(i, 8) = zi * (8.0d0 * ai**1.5d0 * x1 * x12 / x3**2 + 20.0d0 * x1 * x12 * x6 / x3**3 - 6.0d0 * x1 * x3**(-3.5d0) * erf(sqrt(x3) * x6) + 2.0d0 * x14 * (x1 * x9 - 1.0d0) + 4.0d0 * x14 * (x1 * x16 + x36 - 1.0d0) + x33 * (x1 * x16 + x36 - 1.0d0) - x8 * (x1 * x9 - 1.0d0))
            T(i, 9) = yi * (8.0d0 * ai**2.5d0 * x12 * x2 * x4 - x15 - x2 * x35 + x2 * x39 + x2 * x40 - x33 + x8)
            T(i, 10) = zi * (x15 * (x2 * x9 - 1.0d0) + x15 * (x16 * x2 + x2 * x4 - 1.0d0) + x23 * (x18 * x2 + x19 * x2 + x2 * x20 * x4 + x21) - x8 * (5.0d0 * x2 * x4 - 3.0d0))
        end do
    end subroutine T3_damp_erf
    pure subroutine T4_damp_erf(x, y, z, a, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(in) :: a(:)
        real(8), intent(inout) :: T(size(x), 35)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: ai
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x7
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x14
        real(8) :: x15
        real(8) :: x16
        real(8) :: x17
        real(8) :: x18
        real(8) :: x19
        real(8) :: x21
        real(8) :: x22
        real(8) :: x23
        real(8) :: x24
        real(8) :: x25
        real(8) :: x29
        real(8) :: x31
        real(8) :: x32
        real(8) :: x33
        real(8) :: x34
        real(8) :: x35
        real(8) :: x36
        real(8) :: x38
        real(8) :: x39
        real(8) :: x41
        real(8) :: x43
        real(8) :: x44
        real(8) :: x46
        real(8) :: x47
        real(8) :: x51
        real(8) :: x53
        real(8) :: x54
        real(8) :: x55
        real(8) :: x58
        real(8) :: x59
        real(8) :: x60
        real(8) :: x61
        real(8) :: x63
        real(8) :: x65
        real(8) :: x67
        real(8) :: x68
        real(8) :: x71
        real(8) :: x80
        real(8) :: x81
        real(8) :: x83
        real(8) :: x98
        real(8) :: x100
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0 / x3
            x5 = x0 * x4
            x7 = x3**(-3)
            x9 = 0.564189583547756d0 * exp(-ai * x3)
            x10 = sqrt(ai) * x9
            x11 = x10 * x7
            x12 = 24.0d0 * x11
            x14 = x3**(-2)
            x15 = 35.0d0 * x14
            x16 = erf(sqrt(ai) * sqrt(x3))
            x17 = 3.0d0 * x16 * x3**(-2.5d0)
            x18 = 3.0d0 * x4
            x19 = x0 * x18 - 1.0d0
            x21 = 2.0d0 * ai * x0 + x5 - 1.0d0
            x22 = x10 * x14
            x23 = 12.0d0 * x22
            x24 = x0 * x14
            x25 = ai**2
            x29 = -6.0d0 * ai - x18
            x31 = 8.0d0 * x10
            x32 = 8.0d0 * ai**3
            x33 = 15.0d0 * x7
            x34 = 12.0d0 * x25 * x4
            x35 = 18.0d0 * ai * x14
            x36 = 6.0d0 * ai + x18
            x38 = 2.0d0 * x10 * x4
            x39 = x16 * x3**(-4.5d0)
            x41 = x10 / x3**4
            x43 = x16 * x3**(-3.5d0)
            x44 = 15.0d0 * x43
            x46 = ai**1.5d0 * x9
            x47 = x14 * x46
            x51 = x4 * x46
            x53 = x19 * x44
            x54 = 60.0d0 * x0 * x39
            x55 = 152.0d0 * x0 * x41
            x58 = 80.0d0 * x0 * x46 * x7
            x59 = 16.0d0 * x21
            x60 = x11 * x59
            x61 = x1 * x14
            x63 = ai**2.5d0 * x9
            x65 = 16.0d0 * x24 * x63
            x67 = 14.0d0 * x11 * x19
            x68 = x1 * x4
            x71 = 4.0d0 * x19
            x80 = yi * zi
            x81 = x14 * x2
            x83 = x2 * x4
            x98 = x1 * x18 - 1.0d0
            x100 = 2.0d0 * ai * x1 + x68 - 1.0d0
            T(i, 1) = -x0 * x12 * (5.0d0 * x5 - 3.0d0) + x17 * (x15 * xi**4 - 30.0d0 * x5 + 3.0d0) - x19 * x21 * x23 - x24 * x31 * (4.0d0 * ai * x5 + 4.0d0 * x0 * x25 + 3.0d0 * x24 + x29) - x38 * (-24.0d0 * ai * x5 - 24.0d0 * x0 * x25 - 18.0d0 * x24 + x32 * xi**4 + x33 * xi**4 + x34 * xi**4 + x35 * xi**4 + x36)
            T(i, 2) = -xi * yi * (-30.0d0 * x0 * x39 + 48.0d0 * x0 * x41 + 6.0d0 * x11 * (5.0d0 * x5 - 3.0d0) + 4.0d0 * x11 * (4.0d0 * ai * x0 + 6.0d0 * x5 - 3.0d0) + x12 * x19 + x12 * x21 + 12.0d0 * x19 * x47 + 12.0d0 * x21 * x47 + 4.0d0 * x22 * (4.0d0 * ai * x5 + 4.0d0 * x0 * x25 + 3.0d0 * x24 + x29) - x44 * (5.0d0 * x5 - 3.0d0) + 4.0d0 * x51 * (4.0d0 * ai * x5 + 4.0d0 * x0 * x25 + 3.0d0 * x24 + x29))
            T(i, 3) = -xi * zi * (-30.0d0 * x0 * x39 + 48.0d0 * x0 * x41 + 6.0d0 * x11 * (5.0d0 * x5 - 3.0d0) + 4.0d0 * x11 * (4.0d0 * ai * x0 + 6.0d0 * x5 - 3.0d0) + x12 * x19 + x12 * x21 + 12.0d0 * x19 * x47 + 12.0d0 * x21 * x47 + 4.0d0 * x22 * (4.0d0 * ai * x5 + 4.0d0 * x0 * x25 + 3.0d0 * x24 + x29) - x44 * (5.0d0 * x5 - 3.0d0) + 4.0d0 * x51 * (4.0d0 * ai * x5 + 4.0d0 * x0 * x25 + 3.0d0 * x24 + x29))
            T(i, 4) = 20.0d0 * x0 * x11 - 6.0d0 * x0 * x43 + x1 * x53 + x1 * x54 - x1 * x55 - x1 * x58 - x1 * x60 - x1 * x65 - x1 * x67 + 2.0d0 * x10 * x14 * x19 - x17 * x19 + 4.0d0 * x21 * x22 + 4.0d0 * x21 * x51 - 8.0d0 * x21 * x63 * x68 + 8.0d0 * x24 * x46 - x46 * x59 * x61 - x46 * x61 * x71
            T(i, 5) = -x80 * (8.0d0 * x21 * x4 * x63 + x47 * x59 + x47 * x71 - x53 - x54 + x55 + x58 + x60 + x65 + x67)
            T(i, 6) = 20.0d0 * x0 * x11 - 6.0d0 * x0 * x43 + 2.0d0 * x10 * x14 * x19 - x17 * x19 + x2 * x53 + x2 * x54 - x2 * x55 - x2 * x58 - x2 * x60 - x2 * x65 - x2 * x67 + 4.0d0 * x21 * x22 + 4.0d0 * x21 * x51 - 8.0d0 * x21 * x63 * x83 + 8.0d0 * x24 * x46 - x46 * x59 * x81 - x46 * x71 * x81
            T(i, 7) = -xi * yi * (16.0d0 * ai**3.5d0 * x68 * x9 - 105.0d0 * x1 * x39 + 210.0d0 * x1 * x41 + 140.0d0 * x1 * x46 * x7 - 90.0d0 * x11 - 24.0d0 * x4 * x63 + 45.0d0 * x43 - 60.0d0 * x47 + 56.0d0 * x61 * x63)
            T(i, 8) = -xi * zi * (16.0d0 * ai**3.5d0 * x68 * x9 - 105.0d0 * x1 * x39 + 210.0d0 * x1 * x41 + 140.0d0 * x1 * x46 * x7 - 30.0d0 * x11 - 8.0d0 * x4 * x63 + x44 - 20.0d0 * x47 + 56.0d0 * x61 * x63)
            T(i, 9) = -xi * yi * (16.0d0 * ai**3.5d0 * x83 * x9 - 30.0d0 * x11 - 105.0d0 * x2 * x39 + 210.0d0 * x2 * x41 + 140.0d0 * x2 * x46 * x7 - 8.0d0 * x4 * x63 + x44 - 20.0d0 * x47 + 56.0d0 * x63 * x81)
            T(i, 10) = -xi * zi * (16.0d0 * ai**3.5d0 * x83 * x9 - 90.0d0 * x11 - 105.0d0 * x2 * x39 + 210.0d0 * x2 * x41 + 140.0d0 * x2 * x46 * x7 - 24.0d0 * x4 * x63 + 45.0d0 * x43 - 60.0d0 * x47 + 56.0d0 * x63 * x81)
            T(i, 11) = -x1 * x12 * (5.0d0 * x68 - 3.0d0) - x100 * x23 * x98 + x17 * (x15 * yi**4 - 30.0d0 * x68 + 3.0d0) - x31 * x61 * (4.0d0 * ai * x68 + 4.0d0 * x1 * x25 + x29 + 3.0d0 * x61) - x38 * (-24.0d0 * ai * x68 - 24.0d0 * x1 * x25 + x32 * yi**4 + x33 * yi**4 + x34 * yi**4 + x35 * yi**4 + x36 - 18.0d0 * x61)
            T(i, 12) = -x80 * (-30.0d0 * x1 * x39 + 48.0d0 * x1 * x41 + x100 * x12 + 12.0d0 * x100 * x47 + 6.0d0 * x11 * (5.0d0 * x68 - 3.0d0) + 4.0d0 * x11 * (4.0d0 * ai * x1 + 6.0d0 * x68 - 3.0d0) + x12 * x98 + 4.0d0 * x22 * (4.0d0 * ai * x68 + 4.0d0 * x1 * x25 + x29 + 3.0d0 * x61) - x44 * (5.0d0 * x68 - 3.0d0) + 12.0d0 * x47 * x98 + 4.0d0 * x51 * (4.0d0 * ai * x68 + 4.0d0 * x1 * x25 + x29 + 3.0d0 * x61))
            T(i, 13) = 20.0d0 * x1 * x11 + 60.0d0 * x1 * x2 * x39 - 152.0d0 * x1 * x2 * x41 - 80.0d0 * x1 * x2 * x46 * x7 - 6.0d0 * x1 * x43 + 2.0d0 * x10 * x14 * x98 - 16.0d0 * x100 * x11 * x2 + 4.0d0 * x100 * x22 - 16.0d0 * x100 * x46 * x81 + 4.0d0 * x100 * x51 - 8.0d0 * x100 * x63 * x83 - 14.0d0 * x11 * x2 * x98 - x17 * x98 + x2 * x44 * x98 - 16.0d0 * x2 * x61 * x63 + 8.0d0 * x46 * x61 - 4.0d0 * x46 * x81 * x98
            T(i, 14) = -x80 * (16.0d0 * ai**3.5d0 * x83 * x9 - 90.0d0 * x11 - 105.0d0 * x2 * x39 + 210.0d0 * x2 * x41 + 140.0d0 * x2 * x46 * x7 - 24.0d0 * x4 * x63 + 45.0d0 * x43 - 60.0d0 * x47 + 56.0d0 * x63 * x81)
            T(i, 15) = -x12 * x2 * (5.0d0 * x83 - 3.0d0) + x17 * (x15 * zi**4 - 30.0d0 * x83 + 3.0d0) - x23 * (x18 * x2 - 1.0d0) * (2.0d0 * ai * x2 + x83 - 1.0d0) - x31 * x81 * (4.0d0 * ai * x83 + 4.0d0 * x2 * x25 + x29 + 3.0d0 * x81) - x38 * (-24.0d0 * ai * x83 - 24.0d0 * x2 * x25 + x32 * zi**4 + x33 * zi**4 + x34 * zi**4 + x35 * zi**4 + x36 - 18.0d0 * x81)
        end do
    end subroutine T4_damp_erf
    pure subroutine T5_damp_erf(x, y, z, a, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(in) :: a(:)
        real(8), intent(inout) :: T(size(x), 56)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: ai
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x6
        real(8) :: x7
        real(8) :: x8
        real(8) :: x10
        real(8) :: x12
        real(8) :: x14
        real(8) :: x15
        real(8) :: x16
        real(8) :: x17
        real(8) :: x18
        real(8) :: x19
        real(8) :: x20
        real(8) :: x21
        real(8) :: x22
        real(8) :: x23
        real(8) :: x24
        real(8) :: x25
        real(8) :: x26
        real(8) :: x30
        real(8) :: x31
        real(8) :: x32
        real(8) :: x33
        real(8) :: x34
        real(8) :: x35
        real(8) :: x36
        real(8) :: x38
        real(8) :: x39
        real(8) :: x40
        real(8) :: x43
        real(8) :: x44
        real(8) :: x45
        real(8) :: x46
        real(8) :: x47
        real(8) :: x48
        real(8) :: x49
        real(8) :: x50
        real(8) :: x51
        real(8) :: x52
        real(8) :: x54
        real(8) :: x55
        real(8) :: x56
        real(8) :: x58
        real(8) :: x59
        real(8) :: x61
        real(8) :: x62
        real(8) :: x63
        real(8) :: x64
        real(8) :: x65
        real(8) :: x72
        real(8) :: x73
        real(8) :: x75
        real(8) :: x77
        real(8) :: x79
        real(8) :: x80
        real(8) :: x81
        real(8) :: x84
        real(8) :: x86
        real(8) :: x88
        real(8) :: x89
        real(8) :: x90
        real(8) :: x91
        real(8) :: x94
        real(8) :: x97
        real(8) :: x98
        real(8) :: x101
        real(8) :: x102
        real(8) :: x103
        real(8) :: x104
        real(8) :: x105
        real(8) :: x106
        real(8) :: x107
        real(8) :: x109
        real(8) :: x110
        real(8) :: x111
        real(8) :: x112
        real(8) :: x113
        real(8) :: x114
        real(8) :: x115
        real(8) :: x117
        real(8) :: x119
        real(8) :: x120
        real(8) :: x121
        real(8) :: x122
        real(8) :: x126
        real(8) :: x131
        real(8) :: x132
        real(8) :: x134
        real(8) :: x136
        real(8) :: x137
        real(8) :: x138
        real(8) :: x139
        real(8) :: x143
        real(8) :: x148
        real(8) :: x153
        real(8) :: x155
        real(8) :: x166
        real(8) :: x168
        real(8) :: x169
        real(8) :: x171
        real(8) :: x172
        real(8) :: x176
        real(8) :: x177
        real(8) :: x178
        real(8) :: x179
        real(8) :: x180
        real(8) :: x181
        real(8) :: x182
        real(8) :: x189
        real(8) :: x190
        real(8) :: x191
        real(8) :: x196
        real(8) :: x197
        real(8) :: x198
        real(8) :: x200
        real(8) :: x204
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0 / x3
            x5 = x0 * x4
            x6 = x3**(-2)
            x7 = xi**4
            x8 = x6 * x7
            x10 = erf(sqrt(ai) * sqrt(x3))
            x12 = 15.0d0 * x10 * x3**(-3.5d0)
            x14 = x3**(-3)
            x15 = 0.564189583547756d0 * exp(-ai * x3)
            x16 = sqrt(ai) * x15
            x17 = x14 * x16
            x18 = 30.0d0 * x17
            x19 = 5.0d0 * x5 - 3.0d0
            x20 = ai * x0
            x21 = 2.0d0 * x20 + x5 - 1.0d0
            x22 = 60.0d0 * x17
            x23 = 3.0d0 * x4
            x24 = x0 * x23 - 1.0d0
            x25 = x0 * x6
            x26 = ai**2
            x30 = -6.0d0 * ai - x23
            x31 = 4.0d0 * ai * x5 + 4.0d0 * x0 * x26 + 3.0d0 * x25 + x30
            x32 = x16 * x6
            x33 = 20.0d0 * x32
            x34 = ai**3
            x35 = 8.0d0 * x34
            x36 = 15.0d0 * x14
            x38 = 12.0d0 * x26
            x39 = 18.0d0 * ai
            x40 = 6.0d0 * ai + x23
            x43 = 150.0d0 * x14
            x44 = 80.0d0 * x34
            x45 = 16.0d0 * ai**4
            x46 = x3**(-4)
            x47 = 105.0d0 * x46
            x48 = 180.0d0 * ai
            x49 = 120.0d0 * x26
            x50 = 32.0d0 * x34
            x51 = 72.0d0 * x26
            x52 = 120.0d0 * ai * x14
            x54 = 60.0d0 * ai * x4 + 60.0d0 * x26 + 45.0d0 * x6
            x55 = 2.0d0 * x16 * x4
            x56 = x16 / x3**5
            x58 = x10 * x3**(-4.5d0)
            x59 = x0 * x58
            x61 = x16 * x46
            x62 = x0 * x61
            x63 = 24.0d0 * x24
            x64 = ai**1.5d0 * x15
            x65 = x14 * x64
            x72 = 6.0d0 * x17
            x73 = 48.0d0 * x21
            x75 = x6 * x64
            x77 = 24.0d0 * x21
            x79 = 16.0d0 * x31
            x80 = x64 * x79
            x81 = 4.0d0 * x26
            x84 = 4.0d0 * x4 * x64
            x86 = x10 * x3**(-5.5d0)
            x88 = 420.0d0 * x0 * x86
            x89 = 105.0d0 * x58
            x90 = x19 * x89
            x91 = x1 * x4
            x94 = ai * x1
            x97 = 4.0d0 * x17
            x98 = ai**2.5d0 * x15
            x101 = 12.0d0 * x19 * x65
            x102 = 16.0d0 * x61
            x103 = x102 * (4.0d0 * x20 + 6.0d0 * x5 - 3.0d0)
            x104 = x17 * x79
            x105 = 16.0d0 * x65 * (4.0d0 * x20 + 6.0d0 * x5 - 3.0d0)
            x106 = x1 * x6
            x107 = x106 * x98
            x109 = 66.0d0 * x19 * x61
            x110 = 96.0d0 * x65
            x111 = x110 * x24
            x112 = 144.0d0 * x61
            x113 = x112 * x24
            x114 = x112 * x21
            x115 = x46 * x64
            x117 = 192.0d0 * x0 * x115
            x119 = 696.0d0 * x0 * x56
            x120 = 12.0d0 * x75
            x121 = -x120 * x24
            x122 = x110 * x21
            x126 = 4.0d0 * x32
            x131 = 8.0d0 * x4 * x98
            x132 = xi * yi * zi
            x134 = x2 * x4
            x136 = ai * x2
            x137 = x2 * x6
            x138 = x137 * x98
            x139 = 45.0d0 * x10 * x3**(-3.5d0)
            x143 = 24.0d0 * x4 * x98
            x148 = ai**3.5d0 * x15
            x153 = x1 * x65
            x155 = x1 * x61
            x166 = x2 * x89
            x168 = x2 * x65
            x169 = x2 * x61
            x171 = yi**4
            x172 = 945.0d0 * x86
            x176 = x171 * x4
            x177 = 32.0d0 * ai**4.5d0 * x15
            x178 = x171 * x6
            x179 = 144.0d0 * x148
            x180 = 504.0d0 * x14 * x98
            x181 = 1260.0d0 * x115
            x182 = 1890.0d0 * x56
            x189 = x1 * x2
            x190 = x106 * x2
            x191 = zi**4
            x196 = 5.0d0 * x91 - 3.0d0
            x197 = x91 + 2.0d0 * x94 - 1.0d0
            x198 = x1 * x23 - 1.0d0
            x200 = 4.0d0 * ai * x91 + x1 * x81 + 3.0d0 * x106 + x30
            x204 = x1 * x58
            T(i, 1) = xi * (-x12 * (-70.0d0 * x5 + 63.0d0 * x8 + 15.0d0) + x18 * (-30.0d0 * x5 + 35.0d0 * x8 + 3.0d0) + x19 * x21 * x22 + x24 * x31 * x33 + 10.0d0 * x32 * (-24.0d0 * ai * x5 - 24.0d0 * x0 * x26 - 18.0d0 * x25 + x35 * x7 + x36 * x7 + x38 * x4 * x7 + x39 * x8 + x40) + x55 * (-x0 * x43 - x0 * x44 - x25 * x48 + x4 * x50 * x7 + x45 * x7 + x47 * x7 - x49 * x5 + x51 * x8 + x52 * x7 + x54))
            T(i, 2) = yi * (32.0d0 * x0 * x17 * x31 + 48.0d0 * x0 * x19 * x65 - x12 * (-30.0d0 * x5 + 35.0d0 * x8 + 3.0d0) + x17 * x24 * x73 + 12.0d0 * x17 * (12.0d0 * ai * x4 * x7 - 8.0d0 * x20 - 12.0d0 * x5 + x7 * x81 + 15.0d0 * x8 + 1.0d0) + 144.0d0 * x19 * x62 + 72.0d0 * x21 * x62 + x24 * x75 * x77 + x25 * x80 + 4.0d0 * x32 * (-24.0d0 * ai * x5 - 24.0d0 * x0 * x26 - 18.0d0 * x25 + x35 * x7 + x36 * x7 + x38 * x4 * x7 + x39 * x8 + x40) + 240.0d0 * x56 * x7 - 60.0d0 * x59 * (7.0d0 * x5 - 3.0d0) + x62 * x63 + 16.0d0 * x62 * (4.0d0 * x20 + 6.0d0 * x5 - 3.0d0) + x72 * (-30.0d0 * x5 + 35.0d0 * x8 + 3.0d0) + x84 * (-24.0d0 * ai * x5 - 24.0d0 * x0 * x26 - 18.0d0 * x25 + x35 * x7 + x36 * x7 + x38 * x4 * x7 + x39 * x8 + x40))
            T(i, 3) = zi * (32.0d0 * x0 * x17 * x31 + 48.0d0 * x0 * x19 * x65 - x12 * (-30.0d0 * x5 + 35.0d0 * x8 + 3.0d0) + x17 * x24 * x73 + 12.0d0 * x17 * (12.0d0 * ai * x4 * x7 - 8.0d0 * x20 - 12.0d0 * x5 + x7 * x81 + 15.0d0 * x8 + 1.0d0) + 144.0d0 * x19 * x62 + 72.0d0 * x21 * x62 + x24 * x75 * x77 + x25 * x80 + 4.0d0 * x32 * (-24.0d0 * ai * x5 - 24.0d0 * x0 * x26 - 18.0d0 * x25 + x35 * x7 + x36 * x7 + x38 * x4 * x7 + x39 * x8 + x40) + 240.0d0 * x56 * x7 - 60.0d0 * x59 * (7.0d0 * x5 - 3.0d0) + x62 * x63 + 16.0d0 * x62 * (4.0d0 * x20 + 6.0d0 * x5 - 3.0d0) + x72 * (-30.0d0 * x5 + 35.0d0 * x8 + 3.0d0) + x84 * (-24.0d0 * ai * x5 - 24.0d0 * x0 * x26 - 18.0d0 * x25 + x35 * x7 + x36 * x7 + x38 * x4 * x7 + x39 * x8 + x40))
            T(i, 4) = xi * (x1 * x101 + x1 * x103 + x1 * x104 + x1 * x105 + x1 * x109 + x1 * x111 + x1 * x113 + x1 * x114 + x1 * x117 + x1 * x119 + x1 * x122 - x1 * x88 - x1 * x90 + x106 * x80 + x107 * x63 + x107 * x77 + x12 * x19 - x120 * x21 + x121 - x126 * x31 - 24.0d0 * x17 * x21 - 24.0d0 * x17 * x24 - x19 * x72 - x31 * x84 + 8.0d0 * x31 * x91 * x98 + 30.0d0 * x59 - 48.0d0 * x62 + x97 * (36.0d0 * x1 * x25 - 4.0d0 * x20 + 16.0d0 * x5 * x94 - 6.0d0 * x5 - 12.0d0 * x91 + 3.0d0))
            T(i, 5) = x132 * (x101 + x102 * (4.0d0 * x20 + 9.0d0 * x5 - 3.0d0) + x103 + x104 + x105 + x109 + x111 + x113 + x114 + x117 + x119 + x122 + x131 * x31 + 24.0d0 * x21 * x6 * x98 + 24.0d0 * x24 * x6 * x98 + x75 * x79 - x88 - x90)
            T(i, 6) = xi * (x101 * x2 + x103 * x2 + x104 * x2 + x105 * x2 + x109 * x2 + x111 * x2 + x113 * x2 + x114 * x2 + x117 * x2 + x119 * x2 + x12 * x19 - x120 * x21 + x121 + x122 * x2 - x126 * x31 + 8.0d0 * x134 * x31 * x98 + x137 * x80 + x138 * x63 + x138 * x77 - 24.0d0 * x17 * x21 - 24.0d0 * x17 * x24 - x19 * x72 - x2 * x88 - x2 * x90 - x31 * x84 + 30.0d0 * x59 - 48.0d0 * x62 + x97 * (-12.0d0 * x134 + 16.0d0 * x136 * x5 + 36.0d0 * x2 * x25 - 4.0d0 * x20 - 6.0d0 * x5 + 3.0d0))
            T(i, 7) = yi * (840.0d0 * x0 * x1 * x115 + 240.0d0 * x0 * x1 * x14 * x98 + 1452.0d0 * x0 * x1 * x56 - 630.0d0 * x0 * x1 * x86 - 240.0d0 * x0 * x65 + x1 * x122 + 32.0d0 * x1 * x148 * x25 - x1 * x24 * x89 + 8.0d0 * x107 * x24 + x107 * x73 + x121 + x139 * x24 - x143 * x21 + 16.0d0 * x148 * x21 * x91 + 44.0d0 * x153 * x24 + 96.0d0 * x155 * x21 + 114.0d0 * x155 * x24 - 42.0d0 * x17 * x24 - x17 * x73 - 48.0d0 * x25 * x98 + 180.0d0 * x59 - 456.0d0 * x62 - x73 * x75)
            T(i, 8) = zi * (840.0d0 * x0 * x1 * x115 + 240.0d0 * x0 * x1 * x14 * x98 + 1452.0d0 * x0 * x1 * x56 - 630.0d0 * x0 * x1 * x86 - 80.0d0 * x0 * x65 + x1 * x122 + 32.0d0 * x1 * x148 * x25 - x1 * x24 * x89 + 8.0d0 * x107 * x24 + x107 * x73 + x12 * x24 - x131 * x21 + 16.0d0 * x148 * x21 * x91 + 44.0d0 * x153 * x24 + 96.0d0 * x155 * x21 + 114.0d0 * x155 * x24 - 16.0d0 * x17 * x21 - 14.0d0 * x17 * x24 - 16.0d0 * x21 * x75 - 4.0d0 * x24 * x75 - 16.0d0 * x25 * x98 + 60.0d0 * x59 - 152.0d0 * x62)
            T(i, 9) = yi * (840.0d0 * x0 * x115 * x2 + 240.0d0 * x0 * x14 * x2 * x98 + 1452.0d0 * x0 * x2 * x56 - 630.0d0 * x0 * x2 * x86 - 80.0d0 * x0 * x65 + x12 * x24 + x122 * x2 - x131 * x21 + 16.0d0 * x134 * x148 * x21 + 8.0d0 * x138 * x24 + x138 * x73 + 32.0d0 * x148 * x2 * x25 - x166 * x24 + 44.0d0 * x168 * x24 + 96.0d0 * x169 * x21 + 114.0d0 * x169 * x24 - 16.0d0 * x17 * x21 - 14.0d0 * x17 * x24 - 16.0d0 * x21 * x75 - 4.0d0 * x24 * x75 - 16.0d0 * x25 * x98 + 60.0d0 * x59 - 152.0d0 * x62)
            T(i, 10) = zi * (840.0d0 * x0 * x115 * x2 + 240.0d0 * x0 * x14 * x2 * x98 + 1452.0d0 * x0 * x2 * x56 - 630.0d0 * x0 * x2 * x86 - 240.0d0 * x0 * x65 + x121 + x122 * x2 + 16.0d0 * x134 * x148 * x21 + 8.0d0 * x138 * x24 + x138 * x73 + x139 * x24 - x143 * x21 + 32.0d0 * x148 * x2 * x25 - x166 * x24 + 44.0d0 * x168 * x24 + 96.0d0 * x169 * x21 + 114.0d0 * x169 * x24 - 42.0d0 * x17 * x24 - x17 * x73 - 48.0d0 * x25 * x98 + 180.0d0 * x59 - 456.0d0 * x62 - x73 * x75)
            T(i, 11) = xi * (630.0d0 * x1 * x58 - 1260.0d0 * x1 * x61 - 336.0d0 * x107 - x139 + x143 - 96.0d0 * x148 * x91 - 840.0d0 * x153 + 90.0d0 * x17 - x171 * x172 + x171 * x180 + x171 * x181 + x171 * x182 + x176 * x177 + x178 * x179 + 60.0d0 * x75)
            T(i, 12) = x132 * (1260.0d0 * x1 * x115 - x1 * x172 + x1 * x180 + x1 * x182 + x106 * x179 - 48.0d0 * x148 * x4 + x177 * x91 + 315.0d0 * x58 - 168.0d0 * x6 * x98 - 630.0d0 * x61 - 420.0d0 * x65)
            T(i, 13) = xi * (-x1 * x172 * x2 - 210.0d0 * x1 * x61 + x1 * x89 - 56.0d0 * x107 - x12 + x131 - 16.0d0 * x134 * x148 - 56.0d0 * x138 - 16.0d0 * x148 * x91 - 140.0d0 * x153 + x166 - 140.0d0 * x168 + x177 * x2 * x91 + x179 * x190 + x18 + x180 * x189 + x181 * x189 + x182 * x189 - 210.0d0 * x2 * x61 + 20.0d0 * x75)
            T(i, 14) = x132 * (x134 * x177 + x137 * x179 - 48.0d0 * x148 * x4 - x172 * x2 + x180 * x2 + x181 * x2 + x182 * x2 + 315.0d0 * x58 - 168.0d0 * x6 * x98 - 630.0d0 * x61 - 420.0d0 * x65)
            T(i, 15) = xi * (-96.0d0 * x134 * x148 - 336.0d0 * x138 - x139 + x143 - 840.0d0 * x168 - 1260.0d0 * x169 + 90.0d0 * x17 - x172 * x191 + x177 * x191 * x4 + x179 * x191 * x6 + x180 * x191 + x181 * x191 + x182 * x191 + 630.0d0 * x2 * x58 + 60.0d0 * x75)
            T(i, 16) = yi * (-x12 * (63.0d0 * x178 - 70.0d0 * x91 + 15.0d0) + x18 * (35.0d0 * x178 - 30.0d0 * x91 + 3.0d0) + x196 * x197 * x22 + x198 * x200 * x33 + 10.0d0 * x32 * (-24.0d0 * ai * x91 - 24.0d0 * x1 * x26 - 18.0d0 * x106 + x171 * x35 + x171 * x36 + x176 * x38 + x178 * x39 + x40) + x55 * (-x1 * x43 - x1 * x44 - x106 * x48 + x171 * x45 + x171 * x47 + x171 * x52 + x176 * x50 + x178 * x51 - x49 * x91 + x54))
            T(i, 17) = zi * (x1 * x102 * (6.0d0 * x91 + 4.0d0 * x94 - 3.0d0) + x1 * x112 * x196 + 32.0d0 * x1 * x17 * x200 + 16.0d0 * x106 * x200 * x64 - x12 * (35.0d0 * x178 - 30.0d0 * x91 + 3.0d0) + x126 * (-24.0d0 * ai * x91 - 24.0d0 * x1 * x26 - 18.0d0 * x106 + x171 * x35 + x171 * x36 + x176 * x38 + x178 * x39 + x40) + 48.0d0 * x153 * x196 + 72.0d0 * x155 * x197 + 24.0d0 * x155 * x198 + 48.0d0 * x17 * x197 * x198 + 12.0d0 * x17 * (12.0d0 * ai * x171 * x4 + x171 * x81 + 15.0d0 * x178 - 12.0d0 * x91 - 8.0d0 * x94 + 1.0d0) + 240.0d0 * x171 * x56 + 24.0d0 * x197 * x198 * x75 - 60.0d0 * x204 * (7.0d0 * x91 - 3.0d0) + x72 * (35.0d0 * x178 - 30.0d0 * x91 + 3.0d0) + x84 * (-24.0d0 * ai * x91 - 24.0d0 * x1 * x26 - 18.0d0 * x106 + x171 * x35 + x171 * x36 + x176 * x38 + x178 * x39 + x40))
            T(i, 18) = yi * (x102 * x2 * (6.0d0 * x91 + 4.0d0 * x94 - 3.0d0) + x110 * x197 * x2 + x110 * x198 * x2 + x112 * x197 * x2 + x112 * x198 * x2 + 192.0d0 * x115 * x189 + x12 * x196 - x120 * x197 - x120 * x198 - x126 * x200 + 8.0d0 * x134 * x200 * x98 + 16.0d0 * x137 * x200 * x64 + 24.0d0 * x138 * x197 + 24.0d0 * x138 * x198 - 48.0d0 * x155 - x166 * x196 + 12.0d0 * x168 * x196 + 16.0d0 * x168 * (6.0d0 * x91 + 4.0d0 * x94 - 3.0d0) - 24.0d0 * x17 * x197 - 24.0d0 * x17 * x198 + 16.0d0 * x17 * x2 * x200 + 696.0d0 * x189 * x56 - 420.0d0 * x189 * x86 + 66.0d0 * x196 * x2 * x61 - x196 * x72 - x200 * x84 + 30.0d0 * x204 + x97 * (-12.0d0 * x134 + 16.0d0 * x136 * x91 + 36.0d0 * x190 - 6.0d0 * x91 - 4.0d0 * x94 + 3.0d0))
            T(i, 19) = zi * (-48.0d0 * x107 + x110 * x197 * x2 + 840.0d0 * x115 * x189 - x120 * x198 + 16.0d0 * x134 * x148 * x197 + 48.0d0 * x138 * x197 + 8.0d0 * x138 * x198 + x139 * x198 + 240.0d0 * x14 * x189 * x98 - x143 * x197 + 32.0d0 * x148 * x190 - 240.0d0 * x153 - 456.0d0 * x155 - x166 * x198 + 44.0d0 * x168 * x198 + 96.0d0 * x169 * x197 + 114.0d0 * x169 * x198 - 48.0d0 * x17 * x197 - 42.0d0 * x17 * x198 + 1452.0d0 * x189 * x56 - 630.0d0 * x189 * x86 - 48.0d0 * x197 * x75 + 180.0d0 * x204)
            T(i, 20) = yi * (-96.0d0 * x134 * x148 - 336.0d0 * x138 - x139 + x143 - 840.0d0 * x168 - 1260.0d0 * x169 + 90.0d0 * x17 - x172 * x191 + x177 * x191 * x4 + x179 * x191 * x6 + x180 * x191 + x181 * x191 + x182 * x191 + 630.0d0 * x2 * x58 + 60.0d0 * x75)
            T(i, 21) = zi * (-x12 * (-70.0d0 * x134 + 63.0d0 * x191 * x6 + 15.0d0) + x18 * (-30.0d0 * x134 + 35.0d0 * x191 * x6 + 3.0d0) + x22 * (5.0d0 * x134 - 3.0d0) * (x134 + 2.0d0 * x136 - 1.0d0) + 10.0d0 * x32 * (-24.0d0 * ai * x134 - 18.0d0 * x137 + x191 * x35 + x191 * x36 + x191 * x38 * x4 + x191 * x39 * x6 - 24.0d0 * x2 * x26 + x40) + x33 * (x2 * x23 - 1.0d0) * (4.0d0 * ai * x134 + 3.0d0 * x137 + x2 * x81 + x30) + x55 * (-x134 * x49 - x137 * x48 + x191 * x4 * x50 + x191 * x45 + x191 * x47 + x191 * x51 * x6 + x191 * x52 - x2 * x43 - x2 * x44 + x54))
        end do
    end subroutine T5_damp_erf
    pure subroutine T6_damp_erf(x, y, z, a, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(in) :: a(:)
        real(8), intent(inout) :: T(size(x), 84)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: ai
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x6
        real(8) :: x7
        real(8) :: x8
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x13
        real(8) :: x14
        real(8) :: x15
        real(8) :: x17
        real(8) :: x18
        real(8) :: x19
        real(8) :: x21
        real(8) :: x22
        real(8) :: x23
        real(8) :: x24
        real(8) :: x26
        real(8) :: x27
        real(8) :: x28
        real(8) :: x29
        real(8) :: x30
        real(8) :: x32
        real(8) :: x34
        real(8) :: x35
        real(8) :: x36
        real(8) :: x37
        real(8) :: x39
        real(8) :: x40
        real(8) :: x41
        real(8) :: x42
        real(8) :: x44
        real(8) :: x46
        real(8) :: x47
        real(8) :: x48
        real(8) :: x49
        real(8) :: x50
        real(8) :: x52
        real(8) :: x61
        real(8) :: x62
        real(8) :: x64
        real(8) :: x66
        real(8) :: x67
        real(8) :: x68
        real(8) :: x69
        real(8) :: x70
        real(8) :: x71
        real(8) :: x72
        real(8) :: x73
        real(8) :: x74
        real(8) :: x75
        real(8) :: x76
        real(8) :: x77
        real(8) :: x78
        real(8) :: x79
        real(8) :: x81
        real(8) :: x82
        real(8) :: x84
        real(8) :: x88
        real(8) :: x89
        real(8) :: x92
        real(8) :: x93
        real(8) :: x95
        real(8) :: x96
        real(8) :: x100
        real(8) :: x101
        real(8) :: x107
        real(8) :: x112
        real(8) :: x113
        real(8) :: x114
        real(8) :: x121
        real(8) :: x123
        real(8) :: x125
        real(8) :: x126
        real(8) :: x129
        real(8) :: x130
        real(8) :: x134
        real(8) :: x135
        real(8) :: x136
        real(8) :: x140
        real(8) :: x141
        real(8) :: x142
        real(8) :: x144
        real(8) :: x146
        real(8) :: x148
        real(8) :: x149
        real(8) :: x150
        real(8) :: x151
        real(8) :: x152
        real(8) :: x153
        real(8) :: x156
        real(8) :: x157
        real(8) :: x159
        real(8) :: x161
        real(8) :: x164
        real(8) :: x165
        real(8) :: x166
        real(8) :: x168
        real(8) :: x169
        real(8) :: x170
        real(8) :: x171
        real(8) :: x173
        real(8) :: x174
        real(8) :: x175
        real(8) :: x176
        real(8) :: x177
        real(8) :: x178
        real(8) :: x179
        real(8) :: x181
        real(8) :: x182
        real(8) :: x183
        real(8) :: x184
        real(8) :: x185
        real(8) :: x186
        real(8) :: x187
        real(8) :: x188
        real(8) :: x189
        real(8) :: x190
        real(8) :: x191
        real(8) :: x192
        real(8) :: x193
        real(8) :: x194
        real(8) :: x195
        real(8) :: x197
        real(8) :: x200
        real(8) :: x203
        real(8) :: x205
        real(8) :: x206
        real(8) :: x207
        real(8) :: x209
        real(8) :: x210
        real(8) :: x211
        real(8) :: x213
        real(8) :: x214
        real(8) :: x216
        real(8) :: x217
        real(8) :: x218
        real(8) :: x219
        real(8) :: x220
        real(8) :: x221
        real(8) :: x222
        real(8) :: x223
        real(8) :: x224
        real(8) :: x225
        real(8) :: x226
        real(8) :: x230
        real(8) :: x231
        real(8) :: x232
        real(8) :: x234
        real(8) :: x237
        real(8) :: x239
        real(8) :: x241
        real(8) :: x242
        real(8) :: x243
        real(8) :: x244
        real(8) :: x245
        real(8) :: x247
        real(8) :: x249
        real(8) :: x252
        real(8) :: x253
        real(8) :: x254
        real(8) :: x256
        real(8) :: x258
        real(8) :: x262
        real(8) :: x268
        real(8) :: x270
        real(8) :: x278
        real(8) :: x280
        real(8) :: x281
        real(8) :: x283
        real(8) :: x284
        real(8) :: x285
        real(8) :: x286
        real(8) :: x287
        real(8) :: x295
        real(8) :: x298
        real(8) :: x300
        real(8) :: x301
        real(8) :: x302
        real(8) :: x303
        real(8) :: x306
        real(8) :: x309
        real(8) :: x310
        real(8) :: x311
        real(8) :: x314
        real(8) :: x315
        real(8) :: x316
        real(8) :: x317
        real(8) :: x319
        real(8) :: x321
        real(8) :: x326
        real(8) :: x330
        real(8) :: x332
        real(8) :: x342
        real(8) :: x351
        real(8) :: x352
        real(8) :: x362
        real(8) :: x366
        real(8) :: x367
        real(8) :: x371
        real(8) :: x374
        real(8) :: x376
        real(8) :: x377
        real(8) :: x380
        real(8) :: x381
        real(8) :: x384
        real(8) :: x386
        real(8) :: x389
        real(8) :: x390
        real(8) :: x391
        real(8) :: x392
        real(8) :: x393
        real(8) :: x396
        real(8) :: x398
        real(8) :: x403
        real(8) :: x405
        real(8) :: x406
        real(8) :: x407
        real(8) :: x409
        real(8) :: x411
        real(8) :: x413
        real(8) :: x417
        real(8) :: x432
        real(8) :: x438
        real(8) :: x444
        real(8) :: x450
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0 / x3
            x5 = x0 * x4
            x6 = x3**(-2)
            x7 = xi**4
            x8 = x6 * x7
            x10 = x3**(-4)
            x11 = sqrt(ai)
            x12 = 0.564189583547756d0 * exp(-ai * x3)
            x13 = x11 * x12
            x14 = x10 * x13
            x15 = x0 * x14
            x17 = x3**(-3)
            x18 = 231.0d0 * x17
            x19 = erf(x11 * sqrt(x3))
            x21 = 45.0d0 * x19 * x3**(-3.5d0)
            x22 = -30.0d0 * x5 + 35.0d0 * x8 + 3.0d0
            x23 = ai * x0
            x24 = 2.0d0 * x23 + x5 - 1.0d0
            x26 = x11 * x12 * x17
            x27 = 90.0d0 * x26
            x28 = 5.0d0 * x5 - 3.0d0
            x29 = x0 * x6
            x30 = ai**2
            x32 = ai * x5
            x34 = 3.0d0 * x4
            x35 = -6.0d0 * ai - x34
            x36 = 4.0d0 * x0 * x30 + 3.0d0 * x29 + 4.0d0 * x32 + x35
            x37 = x0 * x17
            x39 = 3.0d0 * x5 - 1.0d0
            x40 = ai**3
            x41 = x40 * x7
            x42 = x17 * x7
            x44 = x30 * x4
            x46 = 6.0d0 * ai + x34
            x47 = 18.0d0 * ai * x8 - 24.0d0 * x0 * x30 - 18.0d0 * x29 - 24.0d0 * x32 + 8.0d0 * x41 + 15.0d0 * x42 + 12.0d0 * x44 * x7 + x46
            x48 = x12 * x6
            x49 = x11 * x48
            x50 = 30.0d0 * x49
            x52 = ai**4
            x61 = 60.0d0 * x30
            x62 = ai * x4
            x64 = 45.0d0 * x6 + x61 + 60.0d0 * x62
            x66 = 12.0d0 * x13
            x67 = 32.0d0 * ai**5
            x68 = x3**(-5)
            x69 = 945.0d0 * x68
            x70 = 80.0d0 * x4 * x52
            x71 = 240.0d0 * x40 * x6
            x72 = 600.0d0 * x17 * x30
            x73 = 1050.0d0 * ai * x10
            x74 = -45.0d0 * x6 - x61 - 60.0d0 * x62
            x75 = x12 * x4
            x76 = 2.0d0 * x11 * x75
            x77 = 9.0d0 * x5
            x78 = x19 * x3**(-5.5d0)
            x79 = x0 * x78
            x81 = x13 * x68
            x82 = x0 * x81
            x84 = 7.0d0 * x5
            x88 = x19 * x3**(-4.5d0)
            x89 = 105.0d0 * x88
            x92 = ai**1.5d0
            x93 = x12 * x17 * x92
            x95 = x14 * x28
            x96 = x28 * x93
            x100 = 4.0d0 * x23 + 6.0d0 * x5 - 3.0d0
            x101 = x14 * x39
            x107 = 4.0d0 * x30
            x112 = x26 * x39
            x113 = x48 * x92
            x114 = x113 * x39
            x121 = x75 * x92
            x123 = x1 * x4
            x125 = -12.0d0 * x123
            x126 = x1 * x29
            x129 = 60.0d0 * x0 * x88
            x130 = x22 * x89
            x134 = ai * x1
            x135 = x134 * x5
            x136 = x123 * x30
            x140 = 12.0d0 * x26
            x141 = 600.0d0 * x79 * (x84 - 3.0d0)
            x142 = x13 / x3**6
            x144 = 4128.0d0 * x142 * x7
            x146 = x12 * x92
            x148 = 960.0d0 * x146 * x68 * x7
            x149 = 66.0d0 * x14 * x22
            x150 = 48.0d0 * x14 * (x107 * x7 - 8.0d0 * x23 - 12.0d0 * x5 + 12.0d0 * x62 * x7 + 15.0d0 * x8 + 1.0d0)
            x151 = 48.0d0 * x107 * x7 - 384.0d0 * x23 - 576.0d0 * x5 + 576.0d0 * x62 * x7 + 720.0d0 * x8 + 48.0d0
            x152 = x1 * x17
            x153 = x146 * x152
            x156 = x125 + 36.0d0 * x126 + 16.0d0 * x135 - 4.0d0 * x23 - 6.0d0 * x5 + 3.0d0
            x157 = 16.0d0 * x14
            x159 = 16.0d0 * x47
            x161 = x1 * x6
            x164 = ai**2.5d0
            x165 = x12 * x164
            x166 = 8.0d0 * x165
            x168 = 1152.0d0 * x28 * x82
            x169 = 864.0d0 * x24 * x82
            x170 = x10 * x146
            x171 = x0 * x170
            x173 = 576.0d0 * x171 * x28
            x174 = 288.0d0 * x39
            x175 = x174 * x82
            x176 = 288.0d0 * x24
            x177 = x101 * x176
            x178 = x171 * x176
            x179 = 240.0d0 * x13
            x181 = x0 * x179 * x68 * (x84 - 3.0d0)
            x182 = 192.0d0 * x15 * x36
            x183 = 192.0d0 * x24
            x184 = x153 * x39
            x185 = 128.0d0 * x100 * x82
            x186 = x146 * x37
            x187 = 128.0d0 * x186 * x36
            x188 = 96.0d0 * x171 * x39
            x189 = x165 * x37
            x190 = 96.0d0 * x189 * x28
            x191 = 64.0d0 * x100 * x171
            x192 = x161 * x165
            x193 = 48.0d0 * x24
            x194 = x193 * x39
            x195 = 32.0d0 * x36
            x197 = 15.0d0 * x19 * x3**(-3.5d0)
            x200 = 16.0d0 * x36
            x203 = 256.0d0 * x23 + 64.0d0 * x77 - 192.0d0
            x205 = x164 * x48
            x206 = x165 * x29
            x207 = 48.0d0 * x14
            x209 = 8.0d0 * x164 * x75
            x210 = yi * zi
            x211 = x2 * x4
            x213 = x2 * x29
            x214 = 4.0d0 * x211
            x216 = ai * x2
            x217 = x216 * x5
            x218 = x211 * x30
            x219 = x17 * x2
            x220 = x146 * x219
            x221 = 3.0d0 - 12.0d0 * x211
            x222 = 36.0d0 * x213 + 16.0d0 * x217 + x221 - 4.0d0 * x23 - 6.0d0 * x5
            x223 = x13 * x219
            x224 = x2 * x6
            x225 = x220 * x39
            x226 = x165 * x224
            x230 = 8.0d0 * x134
            x231 = -4.0d0 * x23 - x77 + 3.0d0
            x232 = 288.0d0 * x93
            x234 = 315.0d0 * x88
            x237 = 48.0d0 * x36
            x239 = 24.0d0 * x164 * x75
            x241 = x19 * x3**(-6.5d0)
            x242 = x0 * x241
            x243 = 5670.0d0 * x1
            x244 = x1 * x78
            x245 = 945.0d0 * x244
            x247 = ai**3.5d0 * x12
            x249 = x152 * x165
            x252 = x161 * x247
            x253 = 48.0d0 * x39
            x254 = x1 * x81
            x256 = x1 * x14
            x258 = x1 * x170
            x262 = 576.0d0 * x10 * x165
            x268 = 9612.0d0 * x142
            x270 = xi * yi
            x278 = xi * zi
            x280 = x2 * x81
            x281 = x170 * x2
            x283 = x2 * x78
            x284 = 945.0d0 * x283
            x285 = x165 * x219
            x286 = x224 * x247
            x287 = x14 * x2
            x295 = x1 * x88
            x298 = yi**4
            x300 = 7560.0d0 * x242
            x301 = 16656.0d0 * x142
            x302 = x0 * x301
            x303 = x298 * x68
            x306 = x10 * x298
            x309 = 1122.0d0 * x81
            x310 = 768.0d0 * x24
            x311 = x310 * x81
            x314 = 640.0d0 * x247 * x37
            x315 = 492.0d0 * x39
            x316 = x17 * x298
            x317 = 384.0d0 * x24
            x319 = x298 * x6
            x321 = 128.0d0 * x24
            x326 = ai**4.5d0 * x12
            x330 = 16.0d0 * x39
            x332 = 96.0d0 * x247
            x342 = x123 * x326
            x351 = 3264.0d0 * x10 * x165
            x352 = 10080.0d0 * x0 * x146 * x68
            x362 = x1 * x2
            x366 = x170 * x362
            x367 = x161 * x2
            x371 = x211 * x326
            x374 = zi**4
            x376 = x374 * x68
            x377 = x10 * x374
            x380 = x17 * x374
            x381 = x374 * x6
            x384 = 64.0d0 * x374
            x386 = 10395.0d0 * x241
            x389 = 352.0d0 * x326
            x390 = 1584.0d0 * x247
            x391 = 5544.0d0 * x165
            x392 = 13860.0d0 * x146
            x393 = 20790.0d0 * x142
            x396 = 48.0d0 * ai**3.5d0 * x75 + 630.0d0 * x14 + 168.0d0 * x205 - x234 + 420.0d0 * x93
            x398 = x362 * x68
            x403 = 180.0d0 * x14
            x405 = -30.0d0 * x123 + 35.0d0 * x319 + 3.0d0
            x406 = x123 + 2.0d0 * x134 - 1.0d0
            x407 = 5.0d0 * x123 - 3.0d0
            x409 = 4.0d0 * ai * x123 + x1 * x107 + 3.0d0 * x161 + x35
            x411 = x1 * x34 - 1.0d0
            x413 = x298 * x40
            x417 = -24.0d0 * ai * x123 + 18.0d0 * ai * x319 - 24.0d0 * x1 * x30 - 18.0d0 * x161 + 12.0d0 * x298 * x44 + 15.0d0 * x316 + 8.0d0 * x413 + x46
            x432 = 6.0d0 * x123 + 4.0d0 * x134 - 3.0d0
            x438 = x113 * x411
            x444 = 48.0d0 * x406
            x450 = x123 * x216
            T(i, 1) = -120.0d0 * x13 * x28 * x36 * x37 - 180.0d0 * x15 * (-70.0d0 * x5 + 63.0d0 * x8 + 15.0d0) + x21 * (x18 * xi**6 + 105.0d0 * x5 - 315.0d0 * x8 - 5.0d0) - x22 * x24 * x27 - x29 * x66 * (120.0d0 * ai * x42 - 80.0d0 * x0 * x40 + 105.0d0 * x10 * x7 - 180.0d0 * x23 * x6 - 120.0d0 * x30 * x5 + 72.0d0 * x30 * x8 - 150.0d0 * x37 + 32.0d0 * x4 * x41 + 16.0d0 * x52 * x7 + x64) - x39 * x47 * x50 - x76 * (-1800.0d0 * ai * x42 + 360.0d0 * x0 * x40 - 1575.0d0 * x10 * x7 + 810.0d0 * x23 * x6 + 540.0d0 * x30 * x5 - 1080.0d0 * x30 * x8 + 675.0d0 * x37 - 480.0d0 * x4 * x41 - 240.0d0 * x52 * x7 + x67 * xi**6 + x69 * xi**6 + x70 * xi**6 + x71 * xi**6 + x72 * xi**6 + x73 * xi**6 + x74)
            T(i, 2) = -xi * yi * (40.0d0 * x100 * x101 + 80.0d0 * x112 * x36 + 20.0d0 * x113 * x47 + 40.0d0 * x114 * x36 + 4.0d0 * x121 * (120.0d0 * ai * x42 - 80.0d0 * x0 * x40 + 105.0d0 * x10 * x7 - 180.0d0 * x23 * x6 - 120.0d0 * x30 * x5 + 72.0d0 * x30 * x8 - 150.0d0 * x37 + 32.0d0 * x4 * x41 + 16.0d0 * x52 * x7 + x64) + 180.0d0 * x14 * x22 + 30.0d0 * x14 * (-70.0d0 * x5 + 63.0d0 * x8 + 15.0d0) + 60.0d0 * x14 * (x107 * x7 - 8.0d0 * x23 - 12.0d0 * x5 + 12.0d0 * x62 * x7 + 15.0d0 * x8 + 1.0d0) + 120.0d0 * x15 * x36 + 60.0d0 * x22 * x93 + 600.0d0 * x24 * x82 + 360.0d0 * x24 * x95 + 120.0d0 * x24 * x96 + 40.0d0 * x26 * x47 + 8.0d0 * x26 * (180.0d0 * ai * x8 + 30.0d0 * ai - x0 * x61 - 225.0d0 * x29 - 180.0d0 * x32 + 45.0d0 * x4 + 16.0d0 * x41 + 210.0d0 * x42 + 72.0d0 * x44 * x7) + 120.0d0 * x28 * x82 + 4.0d0 * x49 * (120.0d0 * ai * x42 - 80.0d0 * x0 * x40 + 105.0d0 * x10 * x7 - 180.0d0 * x23 * x6 - 120.0d0 * x30 * x5 + 72.0d0 * x30 * x8 - 150.0d0 * x37 + 32.0d0 * x4 * x41 + 16.0d0 * x52 * x7 + x64) - 420.0d0 * x79 * (x77 - 5.0d0) + 600.0d0 * x82 * (x84 - 3.0d0) - x89 * (-70.0d0 * x5 + 63.0d0 * x8 + 15.0d0))
            T(i, 3) = -xi * zi * (40.0d0 * x100 * x101 + 80.0d0 * x112 * x36 + 20.0d0 * x113 * x47 + 40.0d0 * x114 * x36 + 4.0d0 * x121 * (120.0d0 * ai * x42 - 80.0d0 * x0 * x40 + 105.0d0 * x10 * x7 - 180.0d0 * x23 * x6 - 120.0d0 * x30 * x5 + 72.0d0 * x30 * x8 - 150.0d0 * x37 + 32.0d0 * x4 * x41 + 16.0d0 * x52 * x7 + x64) + 180.0d0 * x14 * x22 + 30.0d0 * x14 * (-70.0d0 * x5 + 63.0d0 * x8 + 15.0d0) + 60.0d0 * x14 * (x107 * x7 - 8.0d0 * x23 - 12.0d0 * x5 + 12.0d0 * x62 * x7 + 15.0d0 * x8 + 1.0d0) + 120.0d0 * x15 * x36 + 60.0d0 * x22 * x93 + 600.0d0 * x24 * x82 + 360.0d0 * x24 * x95 + 120.0d0 * x24 * x96 + 40.0d0 * x26 * x47 + 8.0d0 * x26 * (180.0d0 * ai * x8 + 30.0d0 * ai - x0 * x61 - 225.0d0 * x29 - 180.0d0 * x32 + 45.0d0 * x4 + 16.0d0 * x41 + 210.0d0 * x42 + 72.0d0 * x44 * x7) + 120.0d0 * x28 * x82 + 4.0d0 * x49 * (120.0d0 * ai * x42 - 80.0d0 * x0 * x40 + 105.0d0 * x10 * x7 - 180.0d0 * x23 * x6 - 120.0d0 * x30 * x5 + 72.0d0 * x30 * x8 - 150.0d0 * x37 + 32.0d0 * x4 * x41 + 16.0d0 * x52 * x7 + x64) - 420.0d0 * x79 * (x77 - 5.0d0) + 600.0d0 * x82 * (x84 - 3.0d0) - x89 * (-70.0d0 * x5 + 63.0d0 * x8 + 15.0d0))
            T(i, 4) = -x0 * x156 * x157 + 144.0d0 * x0 * x95 + x1 * x130 + x1 * x141 - x1 * x144 - x1 * x148 - x1 * x149 - x1 * x150 - x1 * x168 - x1 * x169 - x1 * x173 - x1 * x175 - x1 * x177 - x1 * x178 - x1 * x181 - x1 * x182 - x1 * x185 - x1 * x187 - x1 * x188 - x1 * x190 - x1 * x191 + 24.0d0 * x114 * x24 + 4.0d0 * x121 * x47 - x123 * x166 * x47 - x126 * x165 * x195 + x129 * (x125 + 42.0d0 * x126 - x84 + 3.0d0) - x13 * x152 * x159 + x13 * x195 * x37 - x140 * (120.0d0 * x1 * x42 - x107 * x7 + 4.0d0 * x123 - 72.0d0 * x126 + 72.0d0 * x134 * x8 - 32.0d0 * x135 + 16.0d0 * x136 * x7 + 8.0d0 * x23 + 12.0d0 * x5 - 12.0d0 * x62 * x7 - 15.0d0 * x8 - 1.0d0) - x146 * x159 * x161 + x146 * x200 * x29 + 72.0d0 * x15 * x24 + 24.0d0 * x15 * x39 - x151 * x153 - 12.0d0 * x153 * x22 + x179 * x68 * x7 - x183 * x184 + 48.0d0 * x186 * x28 - x192 * x194 + x193 * x26 * x39 - x197 * x22 + 6.0d0 * x22 * x26 + 4.0d0 * x47 * x49
            T(i, 5) = -x210 * (16.0d0 * x113 * x47 - x130 - x141 + x144 + x148 + x149 + x150 + x151 * x93 + x168 + x169 + x173 + x175 + x177 + x178 + x181 + x182 + x183 * x39 * x93 + x185 + x187 + x188 + x190 + x191 + x194 * x205 + x195 * x206 + x203 * x82 + x207 * (x107 * x7 - 8.0d0 * x23 - 18.0d0 * x5 + 18.0d0 * x62 * x7 + 30.0d0 * x8 + 1.0d0) + x209 * x47 + 12.0d0 * x22 * x93 + 16.0d0 * x26 * x47 - 360.0d0 * x79 * (x84 - 2.0d0))
            T(i, 6) = -x0 * x157 * x222 + 144.0d0 * x0 * x95 + 24.0d0 * x114 * x24 + 4.0d0 * x121 * x47 + x129 * (-12.0d0 * x211 + 42.0d0 * x213 - x84 + 3.0d0) + x13 * x195 * x37 + x130 * x2 - x140 * (-x107 * x7 + 120.0d0 * x2 * x42 - 72.0d0 * x213 + x214 + 72.0d0 * x216 * x8 - 32.0d0 * x217 + 16.0d0 * x218 * x7 + 8.0d0 * x23 + 12.0d0 * x5 - 12.0d0 * x62 * x7 - 15.0d0 * x8 - 1.0d0) + x141 * x2 - x144 * x2 - x146 * x159 * x224 + x146 * x200 * x29 - x148 * x2 - x149 * x2 + 72.0d0 * x15 * x24 + 24.0d0 * x15 * x39 - x150 * x2 - x151 * x220 - x159 * x223 - x165 * x195 * x213 - x166 * x211 * x47 - x168 * x2 - x169 * x2 - x173 * x2 - x175 * x2 - x177 * x2 - x178 * x2 + x179 * x68 * x7 - x181 * x2 - x182 * x2 - x183 * x225 - x185 * x2 + 48.0d0 * x186 * x28 - x187 * x2 - x188 * x2 - x190 * x2 - x191 * x2 + x193 * x26 * x39 - x194 * x226 - x197 * x22 - 12.0d0 * x22 * x220 + 6.0d0 * x22 * x26 + 4.0d0 * x47 * x49
            T(i, 7) = -x270 * (3816.0d0 * x0 * x1 * x146 * x68 + x0 * x1 * x262 + x0 * x1 * x268 - 24.0d0 * x100 * x14 + 48.0d0 * x100 * x249 + 96.0d0 * x100 * x254 + 96.0d0 * x100 * x258 - 24.0d0 * x100 * x93 - x113 * x237 + x123 * x200 * x247 + 24.0d0 * x14 * x156 - 432.0d0 * x14 * x24 - 432.0d0 * x14 * x39 + 96.0d0 * x153 * x36 + 24.0d0 * x156 * x93 - 576.0d0 * x171 + x174 * x249 + x176 * x249 + x192 * x237 + x193 * x252 - 72.0d0 * x205 * x24 - 72.0d0 * x205 * x39 + x207 * (-6.0d0 * x123 + 24.0d0 * x126 + x230 * x5 + x231) - x232 * x24 - x232 * x39 + x234 * x28 - x237 * x26 - x239 * x36 + 1152.0d0 * x24 * x254 + 864.0d0 * x24 * x258 - x242 * x243 - x245 * x28 + 24.0d0 * x249 * x28 + x252 * x253 + 738.0d0 * x254 * x28 + 1152.0d0 * x254 * x39 + 96.0d0 * x256 * x36 + 204.0d0 * x258 * x28 + 864.0d0 * x258 * x39 + 1260.0d0 * x79 - 2088.0d0 * x82 - 198.0d0 * x95 - 36.0d0 * x96)
            T(i, 8) = -x278 * (3816.0d0 * x0 * x1 * x146 * x68 + x0 * x1 * x262 + x0 * x1 * x268 - 8.0d0 * x100 * x14 + 48.0d0 * x100 * x249 + 96.0d0 * x100 * x254 + 96.0d0 * x100 * x258 - 8.0d0 * x100 * x93 - x113 * x200 + x123 * x200 * x247 + 8.0d0 * x14 * x156 - 144.0d0 * x14 * x24 - 144.0d0 * x14 * x39 + 96.0d0 * x153 * x36 + 8.0d0 * x156 * x93 + x157 * (-18.0d0 * x123 + 72.0d0 * x126 + 24.0d0 * x135 + x231) - 192.0d0 * x171 + x174 * x249 + x176 * x249 + x192 * x237 + x193 * x252 - x200 * x26 + x203 * x254 + x203 * x258 - 24.0d0 * x205 * x24 - 24.0d0 * x205 * x39 - x209 * x36 + 1152.0d0 * x24 * x254 + 864.0d0 * x24 * x258 - 96.0d0 * x24 * x93 - x242 * x243 - x245 * x28 + 24.0d0 * x249 * x28 + x252 * x253 + 738.0d0 * x254 * x28 + 1152.0d0 * x254 * x39 + 96.0d0 * x256 * x36 + 204.0d0 * x258 * x28 + 864.0d0 * x258 * x39 + x28 * x89 - 96.0d0 * x39 * x93 + 420.0d0 * x79 - 696.0d0 * x82 - 66.0d0 * x95 - 12.0d0 * x96)
            T(i, 9) = -x270 * (3816.0d0 * x0 * x146 * x2 * x68 + x0 * x2 * x262 + x0 * x2 * x268 - 8.0d0 * x100 * x14 + 96.0d0 * x100 * x280 + 96.0d0 * x100 * x281 + 48.0d0 * x100 * x285 - 8.0d0 * x100 * x93 - x113 * x200 + 8.0d0 * x14 * x222 - 144.0d0 * x14 * x24 - 144.0d0 * x14 * x39 + x157 * (-18.0d0 * x211 + 72.0d0 * x213 + 24.0d0 * x217 + x231) - 192.0d0 * x171 + x174 * x285 + x176 * x285 + x193 * x286 - 5670.0d0 * x2 * x242 + x200 * x211 * x247 - x200 * x26 + x203 * x280 + x203 * x281 - 24.0d0 * x205 * x24 - 24.0d0 * x205 * x39 - x209 * x36 + 96.0d0 * x220 * x36 + 8.0d0 * x222 * x93 + x226 * x237 + 1152.0d0 * x24 * x280 + 864.0d0 * x24 * x281 - 96.0d0 * x24 * x93 + x253 * x286 + 738.0d0 * x28 * x280 + 204.0d0 * x28 * x281 - x28 * x284 + 24.0d0 * x28 * x285 + x28 * x89 + 1152.0d0 * x280 * x39 + 864.0d0 * x281 * x39 + 96.0d0 * x287 * x36 - 96.0d0 * x39 * x93 + 420.0d0 * x79 - 696.0d0 * x82 - 66.0d0 * x95 - 12.0d0 * x96)
            T(i, 10) = -x278 * (3816.0d0 * x0 * x146 * x2 * x68 + x0 * x2 * x262 + x0 * x2 * x268 - 24.0d0 * x100 * x14 + 96.0d0 * x100 * x280 + 96.0d0 * x100 * x281 + 48.0d0 * x100 * x285 - 24.0d0 * x100 * x93 - x113 * x237 + 24.0d0 * x14 * x222 - 432.0d0 * x14 * x24 - 432.0d0 * x14 * x39 - 576.0d0 * x171 + x174 * x285 + x176 * x285 + x193 * x286 - 5670.0d0 * x2 * x242 + x200 * x211 * x247 - 72.0d0 * x205 * x24 - 72.0d0 * x205 * x39 + x207 * (-6.0d0 * x211 + 24.0d0 * x213 + 8.0d0 * x217 + x231) + 96.0d0 * x220 * x36 + 24.0d0 * x222 * x93 + x226 * x237 - x232 * x24 - x232 * x39 + x234 * x28 - x237 * x26 - x239 * x36 + 1152.0d0 * x24 * x280 + 864.0d0 * x24 * x281 + x253 * x286 + 738.0d0 * x28 * x280 + 204.0d0 * x28 * x281 - x28 * x284 + 24.0d0 * x28 * x285 + 1152.0d0 * x280 * x39 + 864.0d0 * x281 * x39 + 96.0d0 * x287 * x36 + 1260.0d0 * x79 - 2088.0d0 * x82 - 198.0d0 * x95 - 36.0d0 * x96)
            T(i, 11) = -32.0d0 * ai**4.5d0 * x24 * x298 * x75 - 10080.0d0 * x0 * x146 * x303 - 3264.0d0 * x0 * x165 * x306 + 180.0d0 * x0 * x88 + 684.0d0 * x1 * x101 + 5040.0d0 * x1 * x171 + 1440.0d0 * x1 * x189 - 3780.0d0 * x1 * x79 + 8712.0d0 * x1 * x82 - 42.0d0 * x112 - x113 * x193 - 12.0d0 * x114 + x123 * x24 * x332 + 192.0d0 * x126 * x247 - x146 * x306 * x310 - x146 * x306 * x315 - 456.0d0 * x15 + 576.0d0 * x153 * x24 - x165 * x316 * x317 - 120.0d0 * x165 * x316 * x39 + x176 * x192 + 264.0d0 * x184 - 240.0d0 * x186 + 48.0d0 * x192 * x39 - x193 * x26 - 48.0d0 * x206 + x21 * x39 - x239 * x24 + 576.0d0 * x24 * x256 - x247 * x319 * x321 - x247 * x319 * x330 - 64.0d0 * x29 * x298 * x326 - 630.0d0 * x295 * x39 + x298 * x300 - x298 * x302 - x298 * x309 * x39 - x298 * x311 - x298 * x314 + 945.0d0 * x298 * x39 * x78
            T(i, 12) = -x210 * (-48.0d0 * ai**3.5d0 * x24 * x75 + x0 * x1 * x301 + x0 * x1 * x351 - x1 * x300 + x1 * x314 + x1 * x352 - 342.0d0 * x101 + 64.0d0 * x126 * x326 - x14 * x176 - 2520.0d0 * x171 - 720.0d0 * x189 - 144.0d0 * x205 * x24 - 24.0d0 * x205 * x39 - x232 * x24 + x234 * x39 + 32.0d0 * x24 * x342 - x245 * x39 + x249 * x317 + 120.0d0 * x249 * x39 + x252 * x321 + x252 * x330 + x254 * x310 + 1122.0d0 * x254 * x39 + x258 * x310 + x258 * x315 - x29 * x332 - 132.0d0 * x39 * x93 + 1890.0d0 * x79 - 4356.0d0 * x82)
            T(i, 13) = -x0 * x351 * x362 + 114.0d0 * x1 * x101 + 840.0d0 * x1 * x171 + 240.0d0 * x1 * x189 + 7560.0d0 * x1 * x2 * x242 - x1 * x39 * x89 - 630.0d0 * x1 * x79 + 1452.0d0 * x1 * x82 + 114.0d0 * x101 * x2 - 14.0d0 * x112 - 16.0d0 * x113 * x24 - 4.0d0 * x114 + 16.0d0 * x123 * x24 * x247 - 64.0d0 * x126 * x2 * x326 + 32.0d0 * x126 * x247 + x129 - 152.0d0 * x15 + 96.0d0 * x153 * x24 + x166 * x224 * x39 + 840.0d0 * x171 * x2 + 44.0d0 * x184 - 80.0d0 * x186 + 240.0d0 * x189 * x2 + x192 * x193 + 8.0d0 * x192 * x39 + x193 * x226 + x197 * x39 - 32.0d0 * x2 * x24 * x342 + x2 * x245 * x39 - x2 * x249 * x317 - 120.0d0 * x2 * x249 * x39 - x2 * x39 * x89 - 630.0d0 * x2 * x79 + 1452.0d0 * x2 * x82 - 16.0d0 * x206 - x209 * x24 + 16.0d0 * x211 * x24 * x247 + 32.0d0 * x213 * x247 + 96.0d0 * x220 * x24 + 44.0d0 * x225 + 96.0d0 * x24 * x256 - 16.0d0 * x24 * x26 + 96.0d0 * x24 * x287 - x247 * x321 * x367 - x247 * x330 * x367 - x302 * x362 - x309 * x362 * x39 - x310 * x366 - x311 * x362 - x314 * x362 - x315 * x366 - x352 * x362
            T(i, 14) = -x210 * (-48.0d0 * ai**3.5d0 * x24 * x75 + x0 * x2 * x301 + x0 * x2 * x351 - 342.0d0 * x101 - x14 * x176 - 2520.0d0 * x171 - 720.0d0 * x189 - 7560.0d0 * x2 * x242 + x2 * x314 + x2 * x352 - 144.0d0 * x205 * x24 - 24.0d0 * x205 * x39 + 64.0d0 * x213 * x326 - x232 * x24 + x234 * x39 + 32.0d0 * x24 * x371 + x280 * x310 + 1122.0d0 * x280 * x39 + x281 * x310 + x281 * x315 - x284 * x39 + x285 * x317 + 120.0d0 * x285 * x39 + x286 * x321 + x286 * x330 - x29 * x332 - 132.0d0 * x39 * x93 + 1890.0d0 * x79 - 4356.0d0 * x82)
            T(i, 15) = -32.0d0 * ai**4.5d0 * x24 * x374 * x75 - 10080.0d0 * x0 * x146 * x376 - 3264.0d0 * x0 * x165 * x377 + 180.0d0 * x0 * x88 + 684.0d0 * x101 * x2 - 42.0d0 * x112 - x113 * x193 - 12.0d0 * x114 - x146 * x310 * x377 - x146 * x315 * x377 - 456.0d0 * x15 - x165 * x317 * x380 - 120.0d0 * x165 * x380 * x39 + 5040.0d0 * x171 * x2 + x176 * x226 - 240.0d0 * x186 + 1440.0d0 * x189 * x2 - x193 * x26 - 630.0d0 * x2 * x39 * x88 - 3780.0d0 * x2 * x79 + 8712.0d0 * x2 * x82 - 48.0d0 * x206 + x21 * x39 + x211 * x24 * x332 + 192.0d0 * x213 * x247 + 576.0d0 * x220 * x24 + 264.0d0 * x225 + x226 * x253 - x239 * x24 + 576.0d0 * x24 * x287 - x247 * x321 * x381 - x247 * x330 * x381 - x29 * x326 * x384 + x300 * x374 - x302 * x374 - x309 * x374 * x39 - x311 * x374 - x314 * x374 + 945.0d0 * x374 * x39 * x78
            T(i, 16) = -x270 * (240.0d0 * ai**3.5d0 * x75 + 64.0d0 * ai**5.5d0 * x298 * x75 + 3150.0d0 * x14 + 840.0d0 * x205 + 9450.0d0 * x244 - 5040.0d0 * x249 - 1440.0d0 * x252 - 18900.0d0 * x254 - 12600.0d0 * x258 - x298 * x386 + x298 * x393 + x303 * x392 + x306 * x391 + x316 * x390 + x319 * x389 - 320.0d0 * x342 - 1575.0d0 * x88 + 2100.0d0 * x93)
            T(i, 17) = -x278 * (64.0d0 * ai**5.5d0 * x298 * x75 - 7560.0d0 * x1 * x170 + 5670.0d0 * x244 - 3024.0d0 * x249 - 864.0d0 * x252 - 11340.0d0 * x254 - x298 * x386 + x298 * x393 + x303 * x392 + x306 * x391 + x316 * x390 + x319 * x389 - 192.0d0 * x342 + x396)
            T(i, 18) = -x270 * (64.0d0 * ai**5.5d0 * x12 * x123 * x2 + x10 * x362 * x391 + x152 * x2 * x390 - 3780.0d0 * x170 * x2 + x245 - 504.0d0 * x249 - 144.0d0 * x252 - 1890.0d0 * x254 - 1260.0d0 * x258 - 5670.0d0 * x280 + 2835.0d0 * x283 - 1512.0d0 * x285 - 432.0d0 * x286 - 32.0d0 * x342 - x362 * x386 + x362 * x393 + x367 * x389 - 96.0d0 * x371 + x392 * x398 + x396)
            T(i, 19) = -x278 * (64.0d0 * ai**5.5d0 * x12 * x123 * x2 + x10 * x362 * x391 + x152 * x2 * x390 - x243 * x81 + 2835.0d0 * x244 - 1512.0d0 * x249 - 432.0d0 * x252 - 3780.0d0 * x258 - 1890.0d0 * x280 - 1260.0d0 * x281 + x284 - 504.0d0 * x285 - 144.0d0 * x286 - 96.0d0 * x342 - x362 * x386 + x362 * x393 + x367 * x389 - 32.0d0 * x371 + x392 * x398 + x396)
            T(i, 20) = -x270 * (ai**5.5d0 * x384 * x75 - 11340.0d0 * x280 - 7560.0d0 * x281 + 5670.0d0 * x283 - 3024.0d0 * x285 - 864.0d0 * x286 - 192.0d0 * x371 - x374 * x386 + x374 * x393 + x376 * x392 + x377 * x391 + x380 * x390 + x381 * x389 + x396)
            T(i, 21) = -x278 * (240.0d0 * ai**3.5d0 * x75 + ai**5.5d0 * x384 * x75 + 3150.0d0 * x14 + 840.0d0 * x205 - 18900.0d0 * x280 - 12600.0d0 * x281 + 9450.0d0 * x283 - 5040.0d0 * x285 - 1440.0d0 * x286 - 320.0d0 * x371 - x374 * x386 + x374 * x393 + x376 * x392 + x377 * x391 + x380 * x390 + x381 * x389 - 1575.0d0 * x88 + 2100.0d0 * x93)
            T(i, 22) = -x1 * x403 * (-70.0d0 * x123 + 63.0d0 * x319 + 15.0d0) - 120.0d0 * x13 * x152 * x407 * x409 - x161 * x66 * (120.0d0 * ai * x316 - 80.0d0 * x1 * x40 - 180.0d0 * x134 * x6 - 120.0d0 * x136 - 150.0d0 * x152 + 16.0d0 * x298 * x52 + 72.0d0 * x30 * x319 + 105.0d0 * x306 + 32.0d0 * x4 * x413 + x64) + x21 * (105.0d0 * x123 + x18 * yi**6 - 315.0d0 * x319 - 5.0d0) - x27 * x405 * x406 - x411 * x417 * x50 - x76 * (-1800.0d0 * ai * x316 + 360.0d0 * x1 * x40 + 810.0d0 * x134 * x6 + 540.0d0 * x136 + 675.0d0 * x152 - 240.0d0 * x298 * x52 - 1080.0d0 * x30 * x319 - 1575.0d0 * x306 - 480.0d0 * x4 * x413 + x67 * yi**6 + x69 * yi**6 + x70 * yi**6 + x71 * yi**6 + x72 * yi**6 + x73 * yi**6 + x74)
            T(i, 23) = -x210 * (20.0d0 * x113 * x417 + 4.0d0 * x121 * (120.0d0 * ai * x316 - 80.0d0 * x1 * x40 - 180.0d0 * x134 * x6 - 120.0d0 * x136 - 150.0d0 * x152 + 16.0d0 * x298 * x52 + 72.0d0 * x30 * x319 + 105.0d0 * x306 + 32.0d0 * x4 * x413 + x64) + 360.0d0 * x14 * x406 * x407 + 40.0d0 * x14 * x411 * x432 + 30.0d0 * x14 * (-70.0d0 * x123 + 63.0d0 * x319 + 15.0d0) + 60.0d0 * x14 * (x107 * x298 + x125 - x230 + 12.0d0 * x298 * x62 + 15.0d0 * x319 + 1.0d0) - 420.0d0 * x244 * (9.0d0 * x123 - 5.0d0) + 600.0d0 * x254 * x406 + 120.0d0 * x254 * x407 + 600.0d0 * x254 * (7.0d0 * x123 - 3.0d0) + 120.0d0 * x256 * x409 + 80.0d0 * x26 * x409 * x411 + 40.0d0 * x26 * x417 + 8.0d0 * x26 * (-180.0d0 * ai * x123 + 180.0d0 * ai * x319 + 30.0d0 * ai - x1 * x61 - 225.0d0 * x161 + 72.0d0 * x298 * x44 + 210.0d0 * x316 + 45.0d0 * x4 + 16.0d0 * x413) + x403 * x405 + 60.0d0 * x405 * x93 + 120.0d0 * x406 * x407 * x93 + 40.0d0 * x409 * x438 + 4.0d0 * x49 * (120.0d0 * ai * x316 - 80.0d0 * x1 * x40 - 180.0d0 * x134 * x6 - 120.0d0 * x136 - 150.0d0 * x152 + 16.0d0 * x298 * x52 + 72.0d0 * x30 * x319 + 105.0d0 * x306 + 32.0d0 * x4 * x413 + x64) - x89 * (-70.0d0 * x123 + 63.0d0 * x319 + 15.0d0))
            T(i, 24) = 144.0d0 * x1 * x14 * x407 - x1 * x157 * (-6.0d0 * x123 - 4.0d0 * x134 + x221 + 36.0d0 * x367 + 16.0d0 * x450) + 4.0d0 * x121 * x417 + 32.0d0 * x13 * x152 * x409 - 66.0d0 * x14 * x2 * x405 - 192.0d0 * x14 * x362 * x409 - x140 * (-x107 * x298 + 12.0d0 * x123 + 120.0d0 * x2 * x316 + x214 + 72.0d0 * x216 * x319 + 16.0d0 * x218 * x298 + x230 - 12.0d0 * x298 * x62 - 15.0d0 * x319 - 72.0d0 * x367 - 32.0d0 * x450 - 1.0d0) - 4128.0d0 * x142 * x2 * x298 + 16.0d0 * x146 * x161 * x409 - 960.0d0 * x146 * x2 * x303 - 16.0d0 * x146 * x224 * x417 - 128.0d0 * x153 * x2 * x409 + 48.0d0 * x153 * x407 - 32.0d0 * x165 * x367 * x409 - x166 * x211 * x417 + x179 * x303 - x179 * x398 * (7.0d0 * x123 - 3.0d0) - x197 * x405 - x2 * x207 * (x107 * x298 + x125 - x230 + 12.0d0 * x298 * x62 + 15.0d0 * x319 + 1.0d0) + 600.0d0 * x2 * x244 * (7.0d0 * x123 - 3.0d0) - 96.0d0 * x2 * x249 * x407 + x2 * x405 * x89 - 12.0d0 * x220 * x405 - 192.0d0 * x220 * x406 * x411 - 48.0d0 * x220 * (x107 * x298 + x125 - x230 + 12.0d0 * x298 * x62 + 15.0d0 * x319 + 1.0d0) - 16.0d0 * x223 * x417 - x226 * x411 * x444 + 72.0d0 * x256 * x406 + 24.0d0 * x256 * x411 + 6.0d0 * x26 * x405 + x26 * x411 * x444 - 288.0d0 * x287 * x406 * x411 + 60.0d0 * x295 * (-7.0d0 * x123 + x221 + 42.0d0 * x367) - 864.0d0 * x362 * x406 * x81 - 1152.0d0 * x362 * x407 * x81 - 288.0d0 * x362 * x411 * x81 - 128.0d0 * x362 * x432 * x81 - 288.0d0 * x366 * x406 - 576.0d0 * x366 * x407 - 96.0d0 * x366 * x411 - 64.0d0 * x366 * x432 + 24.0d0 * x406 * x438 + 4.0d0 * x417 * x49
            T(i, 25) = -x210 * (-48.0d0 * x113 * x409 - 432.0d0 * x14 * x406 - 198.0d0 * x14 * x407 - 432.0d0 * x14 * x411 - 24.0d0 * x14 * x432 + 24.0d0 * x14 * (-6.0d0 * x123 - 4.0d0 * x134 + x221 + 36.0d0 * x367 + 16.0d0 * x450) + 3816.0d0 * x146 * x398 - x2 * x241 * x243 - 72.0d0 * x205 * x406 - 72.0d0 * x205 * x411 + x207 * (-9.0d0 * x123 - 4.0d0 * x134 - 6.0d0 * x211 + 24.0d0 * x367 + 8.0d0 * x450 + 3.0d0) + 16.0d0 * x211 * x247 * x409 + 96.0d0 * x220 * x409 + 48.0d0 * x226 * x409 - x232 * x406 - x232 * x411 + x234 * x407 - x239 * x409 + 1260.0d0 * x244 - 2088.0d0 * x254 - 576.0d0 * x258 - 48.0d0 * x26 * x409 + x262 * x362 + x268 * x362 + 1152.0d0 * x280 * x406 + 738.0d0 * x280 * x407 + 1152.0d0 * x280 * x411 + 96.0d0 * x280 * x432 + 864.0d0 * x281 * x406 + 204.0d0 * x281 * x407 + 864.0d0 * x281 * x411 + 96.0d0 * x281 * x432 - x284 * x407 + 288.0d0 * x285 * x406 + 24.0d0 * x285 * x407 + 288.0d0 * x285 * x411 + 48.0d0 * x285 * x432 + 48.0d0 * x286 * x411 + x286 * x444 + 96.0d0 * x287 * x409 - 36.0d0 * x407 * x93 - 24.0d0 * x432 * x93 + 24.0d0 * x93 * (-6.0d0 * x123 - 4.0d0 * x134 + x221 + 36.0d0 * x367 + 16.0d0 * x450))
            T(i, 26) = -32.0d0 * ai**4.5d0 * x374 * x406 * x75 - 10080.0d0 * x1 * x146 * x376 - 3264.0d0 * x1 * x165 * x377 + 7560.0d0 * x1 * x241 * x374 - x1 * x301 * x374 - x113 * x444 - 768.0d0 * x146 * x377 * x406 - 492.0d0 * x146 * x377 * x411 - 640.0d0 * x152 * x247 * x374 - 240.0d0 * x153 - x161 * x326 * x384 - 384.0d0 * x165 * x380 * x406 - 120.0d0 * x165 * x380 * x411 - 48.0d0 * x192 - 3780.0d0 * x2 * x244 + 1440.0d0 * x2 * x249 - 630.0d0 * x2 * x411 * x88 + x21 * x411 + x211 * x332 * x406 + 576.0d0 * x220 * x406 + 264.0d0 * x220 * x411 + 288.0d0 * x226 * x406 + 48.0d0 * x226 * x411 - x239 * x406 + 192.0d0 * x247 * x367 - 128.0d0 * x247 * x381 * x406 - 16.0d0 * x247 * x381 * x411 - 456.0d0 * x256 - 42.0d0 * x26 * x411 - x26 * x444 + 576.0d0 * x287 * x406 + 684.0d0 * x287 * x411 + 180.0d0 * x295 - x309 * x374 * x411 + 8712.0d0 * x362 * x81 + 5040.0d0 * x366 - 768.0d0 * x374 * x406 * x81 + 945.0d0 * x374 * x411 * x78 - 12.0d0 * x438
            T(i, 27) = -x210 * (240.0d0 * ai**3.5d0 * x75 + ai**5.5d0 * x384 * x75 + 3150.0d0 * x14 + 840.0d0 * x205 - 18900.0d0 * x280 - 12600.0d0 * x281 + 9450.0d0 * x283 - 5040.0d0 * x285 - 1440.0d0 * x286 - 320.0d0 * x371 - x374 * x386 + x374 * x393 + x376 * x392 + x377 * x391 + x380 * x390 + x381 * x389 - 1575.0d0 * x88 + 2100.0d0 * x93)
            T(i, 28) = -x2 * x403 * (-70.0d0 * x211 + 63.0d0 * x381 + 15.0d0) + x21 * (x18 * zi**6 + 105.0d0 * x211 - 315.0d0 * x381 - 5.0d0) - 120.0d0 * x223 * (5.0d0 * x211 - 3.0d0) * (ai * x214 + x107 * x2 + 3.0d0 * x224 + x35) - x224 * x66 * (120.0d0 * ai * x380 - 80.0d0 * x2 * x40 - 180.0d0 * x216 * x6 - 120.0d0 * x218 - 150.0d0 * x219 + 72.0d0 * x30 * x381 + 32.0d0 * x374 * x4 * x40 + 16.0d0 * x374 * x52 + 105.0d0 * x377 + x64) - x27 * (-30.0d0 * x211 + 35.0d0 * x381 + 3.0d0) * (x211 + 2.0d0 * x216 - 1.0d0) - x50 * (x2 * x34 - 1.0d0) * (-24.0d0 * ai * x211 + 18.0d0 * ai * x381 - 24.0d0 * x2 * x30 - 18.0d0 * x224 + 8.0d0 * x374 * x40 + 12.0d0 * x374 * x44 + 15.0d0 * x380 + x46) - x76 * (-1800.0d0 * ai * x380 + 360.0d0 * x2 * x40 + 810.0d0 * x216 * x6 + 540.0d0 * x218 + 675.0d0 * x219 - 1080.0d0 * x30 * x381 - 480.0d0 * x374 * x4 * x40 - 240.0d0 * x374 * x52 - 1575.0d0 * x377 + x67 * zi**6 + x69 * zi**6 + x70 * zi**6 + x71 * zi**6 + x72 * zi**6 + x73 * zi**6 + x74)
        end do
    end subroutine T6_damp_erf
end module T_tensor_damp_erf
