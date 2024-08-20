module T_tensor_damp_erf
    use precision
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
        integer(ip), intent(in) :: n
        real(rp), intent(in) :: x(:), y(size(x)), z(size(x))
        real(rp), intent(in) :: a(:)
        real(rp), intent(inout) :: T(size(x), (n + 1) * (n + 2) * (n + 3) / 6)
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
        real(rp), intent(in) :: x(:), y(size(x)), z(size(x))
        real(rp), intent(in) :: a(:)
        real(rp), intent(inout) :: T(size(x), 1)
        integer(ip) :: i
        real(rp) :: xi, yi, zi
        real(rp) :: ai
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            T(i, 1) = (xi**2 + yi**2 + zi**2)**(-0.5_rp) * erf(sqrt(ai) * sqrt(xi**2 + yi**2 + zi**2))
        end do
    end subroutine T0_damp_erf
    pure subroutine T1_damp_erf(x, y, z, a, T)
        real(rp), intent(in) :: x(:), y(size(x)), z(size(x))
        real(rp), intent(in) :: a(:)
        real(rp), intent(inout) :: T(size(x), 4)
        integer(ip) :: i
        real(rp) :: xi, yi, zi
        real(rp) :: ai
        real(rp) :: x2
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x2 = 1.12837916709551_rp * sqrt(ai) * exp(-ai * (xi**2 + yi**2 + zi**2)) / (xi**2 + yi**2 + zi**2) - (xi**2 + yi**2&
& + zi**2)**(-1.5_rp) * erf(sqrt(ai) * sqrt(xi**2 + yi**2 + zi**2))
            T(i, 1) = x2 * xi
            T(i, 2) = x2 * yi
            T(i, 3) = x2 * zi
        end do
    end subroutine T1_damp_erf
    pure subroutine T2_damp_erf(x, y, z, a, T)
        real(rp), intent(in) :: x(:), y(size(x)), z(size(x))
        real(rp), intent(in) :: a(:)
        real(rp), intent(inout) :: T(size(x), 10)
        integer(ip) :: i
        real(rp) :: xi, yi, zi
        real(rp) :: ai
        real(rp) :: x0
        real(rp) :: x1
        real(rp) :: x2
        real(rp) :: x3
        real(rp) :: x4
        real(rp) :: x6
        real(rp) :: x8
        real(rp) :: x9
        real(rp) :: x12
        real(rp) :: x13
        real(rp) :: x14
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 0.564189583547756_rp * exp(-ai * x3)
            x6 = sqrt(ai)
            x8 = 4.0_rp * x4 * x6 / x3**2
            x9 = 1_rp / x3
            x12 = x3**(-1.5_rp) * erf(sqrt(x3) * x6)
            x13 = 2.0_rp * ai
            x14 = 2.0_rp * x4 * x6 * x9
            T(i, 1) = -x0 * x8 + x12 * (3.0_rp * x0 * x9 - 1.0_rp) - x14 * (x0 * x13 + x0 * x9 - 1.0_rp)
            T(i, 2) = -xi * yi * (4.0_rp * ai**1.5_rp * x4 * x9 - 3.0_rp * x3**(-2.5_rp) * erf(sqrt(x3) * x6) + 6.0_rp * x4 * x6 / &
&x3**2)
            T(i, 3) = -xi * zi * (4.0_rp * ai**1.5_rp * x4 * x9 - 3.0_rp * x3**(-2.5_rp) * erf(sqrt(x3) * x6) + 6.0_rp * x4 * x6 / &
&x3**2)
            T(i, 4) = -x1 * x8 + x12 * (3.0_rp * x1 * x9 - 1.0_rp) - x14 * (x1 * x13 + x1 * x9 - 1.0_rp)
            T(i, 5) = -yi * zi * (4.0_rp * ai**1.5_rp * x4 * x9 - 3.0_rp * x3**(-2.5_rp) * erf(sqrt(x3) * x6) + 6.0_rp * x4 * x6 / &
&x3**2)
            T(i, 6) = x12 * (3.0_rp * x2 * x9 - 1.0_rp) - x14 * (x13 * x2 + x2 * x9 - 1.0_rp) - x2 * x8
        end do
    end subroutine T2_damp_erf
    pure subroutine T3_damp_erf(x, y, z, a, T)
        real(rp), intent(in) :: x(:), y(size(x)), z(size(x))
        real(rp), intent(in) :: a(:)
        real(rp), intent(inout) :: T(size(x), 20)
        integer(ip) :: i
        real(rp) :: xi, yi, zi
        real(rp) :: ai
        real(rp) :: x0
        real(rp) :: x1
        real(rp) :: x2
        real(rp) :: x3
        real(rp) :: x4
        real(rp) :: x6
        real(rp) :: x8
        real(rp) :: x9
        real(rp) :: x12
        real(rp) :: x14
        real(rp) :: x15
        real(rp) :: x16
        real(rp) :: x18
        real(rp) :: x19
        real(rp) :: x20
        real(rp) :: x21
        real(rp) :: x22
        real(rp) :: x23
        real(rp) :: x33
        real(rp) :: x35
        real(rp) :: x36
        real(rp) :: x39
        real(rp) :: x40
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1_rp / x3
            x6 = sqrt(ai)
            x8 = 3.0_rp * x3**(-2.5_rp) * erf(sqrt(x3) * x6)
            x9 = 3.0_rp * x4
            x12 = 0.564189583547756_rp * exp(-ai * x3)
            x14 = x12 * x6 / x3**2
            x15 = 6.0_rp * x14
            x16 = 2.0_rp * ai
            x18 = 3.0_rp / x3**2
            x19 = 4.0_rp * ai**2
            x20 = 4.0_rp * ai
            x21 = -6.0_rp * ai - x9
            x22 = x12 * x4
            x23 = 2.0_rp * x22 * x6
            x33 = 4.0_rp * ai**1.5_rp * x22
            x35 = 15.0_rp * x3**(-3.5_rp) * erf(sqrt(x3) * x6)
            x36 = x1 * x4
            x39 = 20.0_rp * ai**1.5_rp * x12 / x3**2
            x40 = 30.0_rp * x12 * x6 / x3**3
            T(i, 1) = xi * (x15 * (x0 * x9 - 1.0_rp) + x15 * (x0 * x16 + x0 * x4 - 1.0_rp) + x23 * (x0 * x18 + x0 * x19 + x0 * x&
&20 * x4 + x21) - x8 * (5.0_rp * x0 * x4 - 3.0_rp))
            T(i, 2) = yi * (8.0_rp * ai**1.5_rp * x0 * x12 / x3**2 + 20.0_rp * x0 * x12 * x6 / x3**3 - 6.0_rp * x0 * x3**(-3.5_rp) &
&* erf(sqrt(x3) * x6) + 2.0_rp * x14 * (x0 * x9 - 1.0_rp) + 4.0_rp * x14 * (x0 * x16 + x0 * x4 - 1.0_rp) + x33 * (x0 * x16 + x0 * x&
&4 - 1.0_rp) - x8 * (x0 * x9 - 1.0_rp))
            T(i, 3) = zi * (8.0_rp * ai**1.5_rp * x0 * x12 / x3**2 + 20.0_rp * x0 * x12 * x6 / x3**3 - 6.0_rp * x0 * x3**(-3.5_rp) &
&* erf(sqrt(x3) * x6) + 2.0_rp * x14 * (x0 * x9 - 1.0_rp) + 4.0_rp * x14 * (x0 * x16 + x0 * x4 - 1.0_rp) + x33 * (x0 * x16 + x0 * x&
&4 - 1.0_rp) - x8 * (x0 * x9 - 1.0_rp))
            T(i, 4) = xi * (8.0_rp * ai**2.5_rp * x12 * x36 - x1 * x35 + x1 * x39 + x1 * x40 - x15 - x33 + x8)
            T(i, 5) = xi * yi * zi * (8.0_rp * ai**2.5_rp * x22 - x35 + x39 + x40)
            T(i, 6) = xi * (8.0_rp * ai**2.5_rp * x12 * x2 * x4 - x15 - x2 * x35 + x2 * x39 + x2 * x40 - x33 + x8)
            T(i, 7) = yi * (x15 * (x1 * x9 - 1.0_rp) + x15 * (x1 * x16 + x36 - 1.0_rp) + x23 * (x1 * x18 + x1 * x19 + x20 * x36 &
&+ x21) - x8 * (5.0_rp * x36 - 3.0_rp))
            T(i, 8) = zi * (8.0_rp * ai**1.5_rp * x1 * x12 / x3**2 + 20.0_rp * x1 * x12 * x6 / x3**3 - 6.0_rp * x1 * x3**(-3.5_rp) &
&* erf(sqrt(x3) * x6) + 2.0_rp * x14 * (x1 * x9 - 1.0_rp) + 4.0_rp * x14 * (x1 * x16 + x36 - 1.0_rp) + x33 * (x1 * x16 + x36 - 1.0_rp) - x8 * (x1 * x9 - 1.0_rp))
            T(i, 9) = yi * (8.0_rp * ai**2.5_rp * x12 * x2 * x4 - x15 - x2 * x35 + x2 * x39 + x2 * x40 - x33 + x8)
            T(i, 10) = zi * (x15 * (x2 * x9 - 1.0_rp) + x15 * (x16 * x2 + x2 * x4 - 1.0_rp) + x23 * (x18 * x2 + x19 * x2 + x2 * &
&x20 * x4 + x21) - x8 * (5.0_rp * x2 * x4 - 3.0_rp))
        end do
    end subroutine T3_damp_erf
    pure subroutine T4_damp_erf(x, y, z, a, T)
        real(rp), intent(in) :: x(:), y(size(x)), z(size(x))
        real(rp), intent(in) :: a(:)
        real(rp), intent(inout) :: T(size(x), 35)
        integer(ip) :: i
        real(rp) :: xi, yi, zi
        real(rp) :: ai
        real(rp) :: x0
        real(rp) :: x1
        real(rp) :: x2
        real(rp) :: x3
        real(rp) :: x4
        real(rp) :: x5
        real(rp) :: x7
        real(rp) :: x9
        real(rp) :: x10
        real(rp) :: x11
        real(rp) :: x12
        real(rp) :: x14
        real(rp) :: x15
        real(rp) :: x16
        real(rp) :: x17
        real(rp) :: x18
        real(rp) :: x19
        real(rp) :: x21
        real(rp) :: x22
        real(rp) :: x23
        real(rp) :: x24
        real(rp) :: x25
        real(rp) :: x29
        real(rp) :: x31
        real(rp) :: x32
        real(rp) :: x33
        real(rp) :: x34
        real(rp) :: x35
        real(rp) :: x36
        real(rp) :: x38
        real(rp) :: x39
        real(rp) :: x41
        real(rp) :: x43
        real(rp) :: x44
        real(rp) :: x46
        real(rp) :: x47
        real(rp) :: x51
        real(rp) :: x53
        real(rp) :: x54
        real(rp) :: x55
        real(rp) :: x58
        real(rp) :: x59
        real(rp) :: x60
        real(rp) :: x61
        real(rp) :: x63
        real(rp) :: x65
        real(rp) :: x67
        real(rp) :: x68
        real(rp) :: x71
        real(rp) :: x80
        real(rp) :: x81
        real(rp) :: x83
        real(rp) :: x98
        real(rp) :: x100
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1_rp / x3
            x5 = x0 * x4
            x7 = x3**(-3)
            x9 = 0.564189583547756_rp * exp(-ai * x3)
            x10 = sqrt(ai) * x9
            x11 = x10 * x7
            x12 = 24.0_rp * x11
            x14 = x3**(-2)
            x15 = 35.0_rp * x14
            x16 = erf(sqrt(ai) * sqrt(x3))
            x17 = 3.0_rp * x16 * x3**(-2.5_rp)
            x18 = 3.0_rp * x4
            x19 = x0 * x18 - 1.0_rp
            x21 = 2.0_rp * ai * x0 + x5 - 1.0_rp
            x22 = x10 * x14
            x23 = 12.0_rp * x22
            x24 = x0 * x14
            x25 = ai**2
            x29 = -6.0_rp * ai - x18
            x31 = 8.0_rp * x10
            x32 = 8.0_rp * ai**3
            x33 = 15.0_rp * x7
            x34 = 12.0_rp * x25 * x4
            x35 = 18.0_rp * ai * x14
            x36 = 6.0_rp * ai + x18
            x38 = 2.0_rp * x10 * x4
            x39 = x16 * x3**(-4.5_rp)
            x41 = x10 / x3**4
            x43 = x16 * x3**(-3.5_rp)
            x44 = 15.0_rp * x43
            x46 = ai**1.5_rp * x9
            x47 = x14 * x46
            x51 = x4 * x46
            x53 = x19 * x44
            x54 = 60.0_rp * x0 * x39
            x55 = 152.0_rp * x0 * x41
            x58 = 80.0_rp * x0 * x46 * x7
            x59 = 16.0_rp * x21
            x60 = x11 * x59
            x61 = x1 * x14
            x63 = ai**2.5_rp * x9
            x65 = 16.0_rp * x24 * x63
            x67 = 14.0_rp * x11 * x19
            x68 = x1 * x4
            x71 = 4.0_rp * x19
            x80 = yi * zi
            x81 = x14 * x2
            x83 = x2 * x4
            x98 = x1 * x18 - 1.0_rp
            x100 = 2.0_rp * ai * x1 + x68 - 1.0_rp
            T(i, 1) = -x0 * x12 * (5.0_rp * x5 - 3.0_rp) + x17 * (x15 * xi**4 - 30.0_rp * x5 + 3.0_rp) - x19 * x21 * x23 - x24 * x&
&31 * (4.0_rp * ai * x5 + 4.0_rp * x0 * x25 + 3.0_rp * x24 + x29) - x38 * (-24.0_rp * ai * x5 - 24.0_rp * x0 * x25 - 18.0_rp * x24 + &
&x32 * xi**4 + x33 * xi**4 + x34 * xi**4 + x35 * xi**4 + x36)
            T(i, 2) = -xi * yi * (-30.0_rp * x0 * x39 + 48.0_rp * x0 * x41 + 6.0_rp * x11 * (5.0_rp * x5 - 3.0_rp) + 4.0_rp * x11 * &
&(4.0_rp * ai * x0 + 6.0_rp * x5 - 3.0_rp) + x12 * x19 + x12 * x21 + 12.0_rp * x19 * x47 + 12.0_rp * x21 * x47 + 4.0_rp * x22 * (4.0_rp * ai * x5 + 4.0_rp * x0 * x25 + 3.0_rp * x24 + x29) - x44 * (5.0_rp * x5 - 3.0_rp) + 4.0_rp * x51 * (4.0_rp * ai * x5 + 4.0_rp * x0&
& * x25 + 3.0_rp * x24 + x29))
            T(i, 3) = -xi * zi * (-30.0_rp * x0 * x39 + 48.0_rp * x0 * x41 + 6.0_rp * x11 * (5.0_rp * x5 - 3.0_rp) + 4.0_rp * x11 * &
&(4.0_rp * ai * x0 + 6.0_rp * x5 - 3.0_rp) + x12 * x19 + x12 * x21 + 12.0_rp * x19 * x47 + 12.0_rp * x21 * x47 + 4.0_rp * x22 * (4.0_rp * ai * x5 + 4.0_rp * x0 * x25 + 3.0_rp * x24 + x29) - x44 * (5.0_rp * x5 - 3.0_rp) + 4.0_rp * x51 * (4.0_rp * ai * x5 + 4.0_rp * x0&
& * x25 + 3.0_rp * x24 + x29))
            T(i, 4) = 20.0_rp * x0 * x11 - 6.0_rp * x0 * x43 + x1 * x53 + x1 * x54 - x1 * x55 - x1 * x58 - x1 * x60 - x1 * x65 -&
& x1 * x67 + 2.0_rp * x10 * x14 * x19 - x17 * x19 + 4.0_rp * x21 * x22 + 4.0_rp * x21 * x51 - 8.0_rp * x21 * x63 * x68 + 8.0_rp * x2&
&4 * x46 - x46 * x59 * x61 - x46 * x61 * x71
            T(i, 5) = -x80 * (8.0_rp * x21 * x4 * x63 + x47 * x59 + x47 * x71 - x53 - x54 + x55 + x58 + x60 + x65 + x67)
            T(i, 6) = 20.0_rp * x0 * x11 - 6.0_rp * x0 * x43 + 2.0_rp * x10 * x14 * x19 - x17 * x19 + x2 * x53 + x2 * x54 - x2 * &
&x55 - x2 * x58 - x2 * x60 - x2 * x65 - x2 * x67 + 4.0_rp * x21 * x22 + 4.0_rp * x21 * x51 - 8.0_rp * x21 * x63 * x83 + 8.0_rp * x2&
&4 * x46 - x46 * x59 * x81 - x46 * x71 * x81
            T(i, 7) = -xi * yi * (16.0_rp * ai**3.5_rp * x68 * x9 - 105.0_rp * x1 * x39 + 210.0_rp * x1 * x41 + 140.0_rp * x1 * x46&
                                 & * x7 - 90.0_rp * x11 - 24.0_rp * x4 * x63 + 45.0_rp * x43 - 60.0_rp * x47 + 56.0_rp * x61 * x63)
            T(i, 8) = -xi * zi * (16.0_rp * ai**3.5_rp * x68 * x9 - 105.0_rp * x1 * x39 + 210.0_rp * x1 * x41 + 140.0_rp * x1 * x46&
                                 & * x7 - 30.0_rp * x11 - 8.0_rp * x4 * x63 + x44 - 20.0_rp * x47 + 56.0_rp * x61 * x63)
            T(i, 9) = -xi * yi * (16.0_rp * ai**3.5_rp * x83 * x9 - 30.0_rp * x11 - 105.0_rp * x2 * x39 + 210.0_rp * x2 * x41 + 140&
&.0_rp * x2 * x46 * x7 - 8.0_rp * x4 * x63 + x44 - 20.0_rp * x47 + 56.0_rp * x63 * x81)
            T(i, 10) = -xi * zi * (16.0_rp * ai**3.5_rp * x83 * x9 - 90.0_rp * x11 - 105.0_rp * x2 * x39 + 210.0_rp * x2 * x41 + 14&
&0.0_rp * x2 * x46 * x7 - 24.0_rp * x4 * x63 + 45.0_rp * x43 - 60.0_rp * x47 + 56.0_rp * x63 * x81)
            T(i, 11) = -x1 * x12 * (5.0_rp * x68 - 3.0_rp) - x100 * x23 * x98 + x17 * (x15 * yi**4 - 30.0_rp * x68 + 3.0_rp) - x31&
                      & * x61 * (4.0_rp * ai * x68 + 4.0_rp * x1 * x25 + x29 + 3.0_rp * x61) - x38 * (-24.0_rp * ai * x68 - 24.0_rp * x1 * x25 + x32 * yi*&
                                                                                                  &*4 + x33 * yi**4 + x34 * yi**4 + x35 * yi**4 + x36 - 18.0_rp * x61)
            T(i, 12) = -x80 * (-30.0_rp * x1 * x39 + 48.0_rp * x1 * x41 + x100 * x12 + 12.0_rp * x100 * x47 + 6.0_rp * x11 * (5.0_rp * x68 - 3.0_rp) + 4.0_rp * x11 * (4.0_rp * ai * x1 + 6.0_rp * x68 - 3.0_rp) + x12 * x98 + 4.0_rp * x22 * (4.0_rp * ai * x68 + 4.0_rp&
& * x1 * x25 + x29 + 3.0_rp * x61) - x44 * (5.0_rp * x68 - 3.0_rp) + 12.0_rp * x47 * x98 + 4.0_rp * x51 * (4.0_rp * ai * x68 + 4.0_rp &
&* x1 * x25 + x29 + 3.0_rp * x61))
            T(i, 13) = 20.0_rp * x1 * x11 + 60.0_rp * x1 * x2 * x39 - 152.0_rp * x1 * x2 * x41 - 80.0_rp * x1 * x2 * x46 * x7 - 6.&
&0_rp * x1 * x43 + 2.0_rp * x10 * x14 * x98 - 16.0_rp * x100 * x11 * x2 + 4.0_rp * x100 * x22 - 16.0_rp * x100 * x46 * x81 + 4.0_rp *&
& x100 * x51 - 8.0_rp * x100 * x63 * x83 - 14.0_rp * x11 * x2 * x98 - x17 * x98 + x2 * x44 * x98 - 16.0_rp * x2 * x61 * x63 + 8.0_rp * x46 * x61 - 4.0_rp * x46 * x81 * x98
            T(i, 14) = -x80 * (16.0_rp * ai**3.5_rp * x83 * x9 - 90.0_rp * x11 - 105.0_rp * x2 * x39 + 210.0_rp * x2 * x41 + 140.0_rp * x2 * x46 * x7 - 24.0_rp * x4 * x63 + 45.0_rp * x43 - 60.0_rp * x47 + 56.0_rp * x63 * x81)
            T(i, 15) = -x12 * x2 * (5.0_rp * x83 - 3.0_rp) + x17 * (x15 * zi**4 - 30.0_rp * x83 + 3.0_rp) - x23 * (x18 * x2 - 1.0_rp) * (2.0_rp * ai * x2 + x83 - 1.0_rp) - x31 * x81 * (4.0_rp * ai * x83 + 4.0_rp * x2 * x25 + x29 + 3.0_rp * x81) - x38 * (-24.0_rp &
&* ai * x83 - 24.0_rp * x2 * x25 + x32 * zi**4 + x33 * zi**4 + x34 * zi**4 + x35 * zi**4 + x36 - 18.0_rp * x81)
        end do
    end subroutine T4_damp_erf
    pure subroutine T5_damp_erf(x, y, z, a, T)
        real(rp), intent(in) :: x(:), y(size(x)), z(size(x))
        real(rp), intent(in) :: a(:)
        real(rp), intent(inout) :: T(size(x), 56)
        integer(ip) :: i
        real(rp) :: xi, yi, zi
        real(rp) :: ai
        real(rp) :: x0
        real(rp) :: x1
        real(rp) :: x2
        real(rp) :: x3
        real(rp) :: x4
        real(rp) :: x5
        real(rp) :: x6
        real(rp) :: x7
        real(rp) :: x8
        real(rp) :: x10
        real(rp) :: x12
        real(rp) :: x14
        real(rp) :: x15
        real(rp) :: x16
        real(rp) :: x17
        real(rp) :: x18
        real(rp) :: x19
        real(rp) :: x20
        real(rp) :: x21
        real(rp) :: x22
        real(rp) :: x23
        real(rp) :: x24
        real(rp) :: x25
        real(rp) :: x26
        real(rp) :: x30
        real(rp) :: x31
        real(rp) :: x32
        real(rp) :: x33
        real(rp) :: x34
        real(rp) :: x35
        real(rp) :: x36
        real(rp) :: x38
        real(rp) :: x39
        real(rp) :: x40
        real(rp) :: x43
        real(rp) :: x44
        real(rp) :: x45
        real(rp) :: x46
        real(rp) :: x47
        real(rp) :: x48
        real(rp) :: x49
        real(rp) :: x50
        real(rp) :: x51
        real(rp) :: x52
        real(rp) :: x54
        real(rp) :: x55
        real(rp) :: x56
        real(rp) :: x58
        real(rp) :: x59
        real(rp) :: x61
        real(rp) :: x62
        real(rp) :: x63
        real(rp) :: x64
        real(rp) :: x65
        real(rp) :: x72
        real(rp) :: x73
        real(rp) :: x75
        real(rp) :: x77
        real(rp) :: x79
        real(rp) :: x80
        real(rp) :: x81
        real(rp) :: x84
        real(rp) :: x86
        real(rp) :: x88
        real(rp) :: x89
        real(rp) :: x90
        real(rp) :: x91
        real(rp) :: x94
        real(rp) :: x97
        real(rp) :: x98
        real(rp) :: x101
        real(rp) :: x102
        real(rp) :: x103
        real(rp) :: x104
        real(rp) :: x105
        real(rp) :: x106
        real(rp) :: x107
        real(rp) :: x109
        real(rp) :: x110
        real(rp) :: x111
        real(rp) :: x112
        real(rp) :: x113
        real(rp) :: x114
        real(rp) :: x115
        real(rp) :: x117
        real(rp) :: x119
        real(rp) :: x120
        real(rp) :: x121
        real(rp) :: x122
        real(rp) :: x126
        real(rp) :: x131
        real(rp) :: x132
        real(rp) :: x134
        real(rp) :: x136
        real(rp) :: x137
        real(rp) :: x138
        real(rp) :: x139
        real(rp) :: x143
        real(rp) :: x148
        real(rp) :: x153
        real(rp) :: x155
        real(rp) :: x166
        real(rp) :: x168
        real(rp) :: x169
        real(rp) :: x171
        real(rp) :: x172
        real(rp) :: x176
        real(rp) :: x177
        real(rp) :: x178
        real(rp) :: x179
        real(rp) :: x180
        real(rp) :: x181
        real(rp) :: x182
        real(rp) :: x189
        real(rp) :: x190
        real(rp) :: x191
        real(rp) :: x196
        real(rp) :: x197
        real(rp) :: x198
        real(rp) :: x200
        real(rp) :: x204
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1_rp / x3
            x5 = x0 * x4
            x6 = x3**(-2)
            x7 = xi**4
            x8 = x6 * x7
            x10 = erf(sqrt(ai) * sqrt(x3))
            x12 = 15.0_rp * x10 * x3**(-3.5_rp)
            x14 = x3**(-3)
            x15 = 0.564189583547756_rp * exp(-ai * x3)
            x16 = sqrt(ai) * x15
            x17 = x14 * x16
            x18 = 30.0_rp * x17
            x19 = 5.0_rp * x5 - 3.0_rp
            x20 = ai * x0
            x21 = 2.0_rp * x20 + x5 - 1.0_rp
            x22 = 60.0_rp * x17
            x23 = 3.0_rp * x4
            x24 = x0 * x23 - 1.0_rp
            x25 = x0 * x6
            x26 = ai**2
            x30 = -6.0_rp * ai - x23
            x31 = 4.0_rp * ai * x5 + 4.0_rp * x0 * x26 + 3.0_rp * x25 + x30
            x32 = x16 * x6
            x33 = 20.0_rp * x32
            x34 = ai**3
            x35 = 8.0_rp * x34
            x36 = 15.0_rp * x14
            x38 = 12.0_rp * x26
            x39 = 18.0_rp * ai
            x40 = 6.0_rp * ai + x23
            x43 = 150.0_rp * x14
            x44 = 80.0_rp * x34
            x45 = 16.0_rp * ai**4
            x46 = x3**(-4)
            x47 = 105.0_rp * x46
            x48 = 180.0_rp * ai
            x49 = 120.0_rp * x26
            x50 = 32.0_rp * x34
            x51 = 72.0_rp * x26
            x52 = 120.0_rp * ai * x14
            x54 = 60.0_rp * ai * x4 + 60.0_rp * x26 + 45.0_rp * x6
            x55 = 2.0_rp * x16 * x4
            x56 = x16 / x3**5
            x58 = x10 * x3**(-4.5_rp)
            x59 = x0 * x58
            x61 = x16 * x46
            x62 = x0 * x61
            x63 = 24.0_rp * x24
            x64 = ai**1.5_rp * x15
            x65 = x14 * x64
            x72 = 6.0_rp * x17
            x73 = 48.0_rp * x21
            x75 = x6 * x64
            x77 = 24.0_rp * x21
            x79 = 16.0_rp * x31
            x80 = x64 * x79
            x81 = 4.0_rp * x26
            x84 = 4.0_rp * x4 * x64
            x86 = x10 * x3**(-5.5_rp)
            x88 = 420.0_rp * x0 * x86
            x89 = 105.0_rp * x58
            x90 = x19 * x89
            x91 = x1 * x4
            x94 = ai * x1
            x97 = 4.0_rp * x17
            x98 = ai**2.5_rp * x15
            x101 = 12.0_rp * x19 * x65
            x102 = 16.0_rp * x61
            x103 = x102 * (4.0_rp * x20 + 6.0_rp * x5 - 3.0_rp)
            x104 = x17 * x79
            x105 = 16.0_rp * x65 * (4.0_rp * x20 + 6.0_rp * x5 - 3.0_rp)
            x106 = x1 * x6
            x107 = x106 * x98
            x109 = 66.0_rp * x19 * x61
            x110 = 96.0_rp * x65
            x111 = x110 * x24
            x112 = 144.0_rp * x61
            x113 = x112 * x24
            x114 = x112 * x21
            x115 = x46 * x64
            x117 = 192.0_rp * x0 * x115
            x119 = 696.0_rp * x0 * x56
            x120 = 12.0_rp * x75
            x121 = -x120 * x24
            x122 = x110 * x21
            x126 = 4.0_rp * x32
            x131 = 8.0_rp * x4 * x98
            x132 = xi * yi * zi
            x134 = x2 * x4
            x136 = ai * x2
            x137 = x2 * x6
            x138 = x137 * x98
            x139 = 45.0_rp * x10 * x3**(-3.5_rp)
            x143 = 24.0_rp * x4 * x98
            x148 = ai**3.5_rp * x15
            x153 = x1 * x65
            x155 = x1 * x61
            x166 = x2 * x89
            x168 = x2 * x65
            x169 = x2 * x61
            x171 = yi**4
            x172 = 945.0_rp * x86
            x176 = x171 * x4
            x177 = 32.0_rp * ai**4.5_rp * x15
            x178 = x171 * x6
            x179 = 144.0_rp * x148
            x180 = 504.0_rp * x14 * x98
            x181 = 1260.0_rp * x115
            x182 = 1890.0_rp * x56
            x189 = x1 * x2
            x190 = x106 * x2
            x191 = zi**4
            x196 = 5.0_rp * x91 - 3.0_rp
            x197 = x91 + 2.0_rp * x94 - 1.0_rp
            x198 = x1 * x23 - 1.0_rp
            x200 = 4.0_rp * ai * x91 + x1 * x81 + 3.0_rp * x106 + x30
            x204 = x1 * x58
            T(i, 1) = xi * (-x12 * (-70.0_rp * x5 + 63.0_rp * x8 + 15.0_rp) + x18 * (-30.0_rp * x5 + 35.0_rp * x8 + 3.0_rp) + x19 * &
&x21 * x22 + x24 * x31 * x33 + 10.0_rp * x32 * (-24.0_rp * ai * x5 - 24.0_rp * x0 * x26 - 18.0_rp * x25 + x35 * x7 + x36 * x7 + x38&
& * x4 * x7 + x39 * x8 + x40) + x55 * (-x0 * x43 - x0 * x44 - x25 * x48 + x4 * x50 * x7 + x45 * x7 + x47 * x7 - x49 * x5 + x51 &
&* x8 + x52 * x7 + x54))
            T(i, 2) = yi * (32.0_rp * x0 * x17 * x31 + 48.0_rp * x0 * x19 * x65 - x12 * (-30.0_rp * x5 + 35.0_rp * x8 + 3.0_rp) + x&
&17 * x24 * x73 + 12.0_rp * x17 * (12.0_rp * ai * x4 * x7 - 8.0_rp * x20 - 12.0_rp * x5 + x7 * x81 + 15.0_rp * x8 + 1.0_rp) + 144.0_rp&
& * x19 * x62 + 72.0_rp * x21 * x62 + x24 * x75 * x77 + x25 * x80 + 4.0_rp * x32 * (-24.0_rp * ai * x5 - 24.0_rp * x0 * x26 - 18.0_rp * x25 + x35 * x7 + x36 * x7 + x38 * x4 * x7 + x39 * x8 + x40) + 240.0_rp * x56 * x7 - 60.0_rp * x59 * (7.0_rp * x5 - 3.0_rp) + x&
&62 * x63 + 16.0_rp * x62 * (4.0_rp * x20 + 6.0_rp * x5 - 3.0_rp) + x72 * (-30.0_rp * x5 + 35.0_rp * x8 + 3.0_rp) + x84 * (-24.0_rp * a&
&i * x5 - 24.0_rp * x0 * x26 - 18.0_rp * x25 + x35 * x7 + x36 * x7 + x38 * x4 * x7 + x39 * x8 + x40))
            T(i, 3) = zi * (32.0_rp * x0 * x17 * x31 + 48.0_rp * x0 * x19 * x65 - x12 * (-30.0_rp * x5 + 35.0_rp * x8 + 3.0_rp) + x&
&17 * x24 * x73 + 12.0_rp * x17 * (12.0_rp * ai * x4 * x7 - 8.0_rp * x20 - 12.0_rp * x5 + x7 * x81 + 15.0_rp * x8 + 1.0_rp) + 144.0_rp&
& * x19 * x62 + 72.0_rp * x21 * x62 + x24 * x75 * x77 + x25 * x80 + 4.0_rp * x32 * (-24.0_rp * ai * x5 - 24.0_rp * x0 * x26 - 18.0_rp * x25 + x35 * x7 + x36 * x7 + x38 * x4 * x7 + x39 * x8 + x40) + 240.0_rp * x56 * x7 - 60.0_rp * x59 * (7.0_rp * x5 - 3.0_rp) + x&
&62 * x63 + 16.0_rp * x62 * (4.0_rp * x20 + 6.0_rp * x5 - 3.0_rp) + x72 * (-30.0_rp * x5 + 35.0_rp * x8 + 3.0_rp) + x84 * (-24.0_rp * a&
&i * x5 - 24.0_rp * x0 * x26 - 18.0_rp * x25 + x35 * x7 + x36 * x7 + x38 * x4 * x7 + x39 * x8 + x40))
            T(i, 4) = xi * (x1 * x101 + x1 * x103 + x1 * x104 + x1 * x105 + x1 * x109 + x1 * x111 + x1 * x113 + x1 * x114 + x1&
& * x117 + x1 * x119 + x1 * x122 - x1 * x88 - x1 * x90 + x106 * x80 + x107 * x63 + x107 * x77 + x12 * x19 - x120 * x21 + x121 -&
& x126 * x31 - 24.0_rp * x17 * x21 - 24.0_rp * x17 * x24 - x19 * x72 - x31 * x84 + 8.0_rp * x31 * x91 * x98 + 30.0_rp * x59 - 48.0_rp * x62 + x97 * (36.0_rp * x1 * x25 - 4.0_rp * x20 + 16.0_rp * x5 * x94 - 6.0_rp * x5 - 12.0_rp * x91 + 3.0_rp))
            T(i, 5) = x132 * (x101 + x102 * (4.0_rp * x20 + 9.0_rp * x5 - 3.0_rp) + x103 + x104 + x105 + x109 + x111 + x113 + x11&
&4 + x117 + x119 + x122 + x131 * x31 + 24.0_rp * x21 * x6 * x98 + 24.0_rp * x24 * x6 * x98 + x75 * x79 - x88 - x90)
            T(i, 6) = xi * (x101 * x2 + x103 * x2 + x104 * x2 + x105 * x2 + x109 * x2 + x111 * x2 + x113 * x2 + x114 * x2 + x1&
&17 * x2 + x119 * x2 + x12 * x19 - x120 * x21 + x121 + x122 * x2 - x126 * x31 + 8.0_rp * x134 * x31 * x98 + x137 * x80 + x138 * &
&x63 + x138 * x77 - 24.0_rp * x17 * x21 - 24.0_rp * x17 * x24 - x19 * x72 - x2 * x88 - x2 * x90 - x31 * x84 + 30.0_rp * x59 - 48.0&
&_rp * x62 + x97 * (-12.0_rp * x134 + 16.0_rp * x136 * x5 + 36.0_rp * x2 * x25 - 4.0_rp * x20 - 6.0_rp * x5 + 3.0_rp))
            T(i, 7) = yi * (840.0_rp * x0 * x1 * x115 + 240.0_rp * x0 * x1 * x14 * x98 + 1452.0_rp * x0 * x1 * x56 - 630.0_rp * x0&
& * x1 * x86 - 240.0_rp * x0 * x65 + x1 * x122 + 32.0_rp * x1 * x148 * x25 - x1 * x24 * x89 + 8.0_rp * x107 * x24 + x107 * x73 + x&
&121 + x139 * x24 - x143 * x21 + 16.0_rp * x148 * x21 * x91 + 44.0_rp * x153 * x24 + 96.0_rp * x155 * x21 + 114.0_rp * x155 * x24 -&
& 42.0_rp * x17 * x24 - x17 * x73 - 48.0_rp * x25 * x98 + 180.0_rp * x59 - 456.0_rp * x62 - x73 * x75)
            T(i, 8) = zi * (840.0_rp * x0 * x1 * x115 + 240.0_rp * x0 * x1 * x14 * x98 + 1452.0_rp * x0 * x1 * x56 - 630.0_rp * x0&
& * x1 * x86 - 80.0_rp * x0 * x65 + x1 * x122 + 32.0_rp * x1 * x148 * x25 - x1 * x24 * x89 + 8.0_rp * x107 * x24 + x107 * x73 + x1&
&2 * x24 - x131 * x21 + 16.0_rp * x148 * x21 * x91 + 44.0_rp * x153 * x24 + 96.0_rp * x155 * x21 + 114.0_rp * x155 * x24 - 16.0_rp *&
& x17 * x21 - 14.0_rp * x17 * x24 - 16.0_rp * x21 * x75 - 4.0_rp * x24 * x75 - 16.0_rp * x25 * x98 + 60.0_rp * x59 - 152.0_rp * x62)
            T(i, 9) = yi * (840.0_rp * x0 * x115 * x2 + 240.0_rp * x0 * x14 * x2 * x98 + 1452.0_rp * x0 * x2 * x56 - 630.0_rp * x0&
& * x2 * x86 - 80.0_rp * x0 * x65 + x12 * x24 + x122 * x2 - x131 * x21 + 16.0_rp * x134 * x148 * x21 + 8.0_rp * x138 * x24 + x138 &
&* x73 + 32.0_rp * x148 * x2 * x25 - x166 * x24 + 44.0_rp * x168 * x24 + 96.0_rp * x169 * x21 + 114.0_rp * x169 * x24 - 16.0_rp * x1&
&7 * x21 - 14.0_rp * x17 * x24 - 16.0_rp * x21 * x75 - 4.0_rp * x24 * x75 - 16.0_rp * x25 * x98 + 60.0_rp * x59 - 152.0_rp * x62)
            T(i, 10) = zi * (840.0_rp * x0 * x115 * x2 + 240.0_rp * x0 * x14 * x2 * x98 + 1452.0_rp * x0 * x2 * x56 - 630.0_rp * x&
&0 * x2 * x86 - 240.0_rp * x0 * x65 + x121 + x122 * x2 + 16.0_rp * x134 * x148 * x21 + 8.0_rp * x138 * x24 + x138 * x73 + x139 * x&
&24 - x143 * x21 + 32.0_rp * x148 * x2 * x25 - x166 * x24 + 44.0_rp * x168 * x24 + 96.0_rp * x169 * x21 + 114.0_rp * x169 * x24 - 4&
&2.0_rp * x17 * x24 - x17 * x73 - 48.0_rp * x25 * x98 + 180.0_rp * x59 - 456.0_rp * x62 - x73 * x75)
            T(i, 11) = xi * (630.0_rp * x1 * x58 - 1260.0_rp * x1 * x61 - 336.0_rp * x107 - x139 + x143 - 96.0_rp * x148 * x91 - 8&
&40.0_rp * x153 + 90.0_rp * x17 - x171 * x172 + x171 * x180 + x171 * x181 + x171 * x182 + x176 * x177 + x178 * x179 + 60.0_rp * x7&
&5)
            T(i, 12) = x132 * (1260.0_rp * x1 * x115 - x1 * x172 + x1 * x180 + x1 * x182 + x106 * x179 - 48.0_rp * x148 * x4 + x&
&177 * x91 + 315.0_rp * x58 - 168.0_rp * x6 * x98 - 630.0_rp * x61 - 420.0_rp * x65)
            T(i, 13) = xi * (-x1 * x172 * x2 - 210.0_rp * x1 * x61 + x1 * x89 - 56.0_rp * x107 - x12 + x131 - 16.0_rp * x134 * x1&
&48 - 56.0_rp * x138 - 16.0_rp * x148 * x91 - 140.0_rp * x153 + x166 - 140.0_rp * x168 + x177 * x2 * x91 + x179 * x190 + x18 + x180&
& * x189 + x181 * x189 + x182 * x189 - 210.0_rp * x2 * x61 + 20.0_rp * x75)
            T(i, 14) = x132 * (x134 * x177 + x137 * x179 - 48.0_rp * x148 * x4 - x172 * x2 + x180 * x2 + x181 * x2 + x182 * x2 &
&+ 315.0_rp * x58 - 168.0_rp * x6 * x98 - 630.0_rp * x61 - 420.0_rp * x65)
            T(i, 15) = xi * (-96.0_rp * x134 * x148 - 336.0_rp * x138 - x139 + x143 - 840.0_rp * x168 - 1260.0_rp * x169 + 90.0_rp &
&* x17 - x172 * x191 + x177 * x191 * x4 + x179 * x191 * x6 + x180 * x191 + x181 * x191 + x182 * x191 + 630.0_rp * x2 * x58 + 60.&
&0_rp * x75)
            T(i, 16) = yi * (-x12 * (63.0_rp * x178 - 70.0_rp * x91 + 15.0_rp) + x18 * (35.0_rp * x178 - 30.0_rp * x91 + 3.0_rp) + x&
&196 * x197 * x22 + x198 * x200 * x33 + 10.0_rp * x32 * (-24.0_rp * ai * x91 - 24.0_rp * x1 * x26 - 18.0_rp * x106 + x171 * x35 + x&
&171 * x36 + x176 * x38 + x178 * x39 + x40) + x55 * (-x1 * x43 - x1 * x44 - x106 * x48 + x171 * x45 + x171 * x47 + x171 * x52 +&
& x176 * x50 + x178 * x51 - x49 * x91 + x54))
            T(i, 17) = zi * (x1 * x102 * (6.0_rp * x91 + 4.0_rp * x94 - 3.0_rp) + x1 * x112 * x196 + 32.0_rp * x1 * x17 * x200 + 1&
&6.0_rp * x106 * x200 * x64 - x12 * (35.0_rp * x178 - 30.0_rp * x91 + 3.0_rp) + x126 * (-24.0_rp * ai * x91 - 24.0_rp * x1 * x26 - 18&
&.0_rp * x106 + x171 * x35 + x171 * x36 + x176 * x38 + x178 * x39 + x40) + 48.0_rp * x153 * x196 + 72.0_rp * x155 * x197 + 24.0_rp &
&* x155 * x198 + 48.0_rp * x17 * x197 * x198 + 12.0_rp * x17 * (12.0_rp * ai * x171 * x4 + x171 * x81 + 15.0_rp * x178 - 12.0_rp * x&
&91 - 8.0_rp * x94 + 1.0_rp) + 240.0_rp * x171 * x56 + 24.0_rp * x197 * x198 * x75 - 60.0_rp * x204 * (7.0_rp * x91 - 3.0_rp) + x72 * &
&(35.0_rp * x178 - 30.0_rp * x91 + 3.0_rp) + x84 * (-24.0_rp * ai * x91 - 24.0_rp * x1 * x26 - 18.0_rp * x106 + x171 * x35 + x171 * x&
&36 + x176 * x38 + x178 * x39 + x40))
            T(i, 18) = yi * (x102 * x2 * (6.0_rp * x91 + 4.0_rp * x94 - 3.0_rp) + x110 * x197 * x2 + x110 * x198 * x2 + x112 * x1&
&97 * x2 + x112 * x198 * x2 + 192.0_rp * x115 * x189 + x12 * x196 - x120 * x197 - x120 * x198 - x126 * x200 + 8.0_rp * x134 * x20&
&0 * x98 + 16.0_rp * x137 * x200 * x64 + 24.0_rp * x138 * x197 + 24.0_rp * x138 * x198 - 48.0_rp * x155 - x166 * x196 + 12.0_rp * x1&
&68 * x196 + 16.0_rp * x168 * (6.0_rp * x91 + 4.0_rp * x94 - 3.0_rp) - 24.0_rp * x17 * x197 - 24.0_rp * x17 * x198 + 16.0_rp * x17 * x&
&2 * x200 + 696.0_rp * x189 * x56 - 420.0_rp * x189 * x86 + 66.0_rp * x196 * x2 * x61 - x196 * x72 - x200 * x84 + 30.0_rp * x204 + &
&x97 * (-12.0_rp * x134 + 16.0_rp * x136 * x91 + 36.0_rp * x190 - 6.0_rp * x91 - 4.0_rp * x94 + 3.0_rp))
            T(i, 19) = zi * (-48.0_rp * x107 + x110 * x197 * x2 + 840.0_rp * x115 * x189 - x120 * x198 + 16.0_rp * x134 * x148 * &
&x197 + 48.0_rp * x138 * x197 + 8.0_rp * x138 * x198 + x139 * x198 + 240.0_rp * x14 * x189 * x98 - x143 * x197 + 32.0_rp * x148 * x&
&190 - 240.0_rp * x153 - 456.0_rp * x155 - x166 * x198 + 44.0_rp * x168 * x198 + 96.0_rp * x169 * x197 + 114.0_rp * x169 * x198 - 48&
&.0_rp * x17 * x197 - 42.0_rp * x17 * x198 + 1452.0_rp * x189 * x56 - 630.0_rp * x189 * x86 - 48.0_rp * x197 * x75 + 180.0_rp * x204)
            T(i, 20) = yi * (-96.0_rp * x134 * x148 - 336.0_rp * x138 - x139 + x143 - 840.0_rp * x168 - 1260.0_rp * x169 + 90.0_rp &
&* x17 - x172 * x191 + x177 * x191 * x4 + x179 * x191 * x6 + x180 * x191 + x181 * x191 + x182 * x191 + 630.0_rp * x2 * x58 + 60.&
&0_rp * x75)
            T(i, 21) = zi * (-x12 * (-70.0_rp * x134 + 63.0_rp * x191 * x6 + 15.0_rp) + x18 * (-30.0_rp * x134 + 35.0_rp * x191 * x&
&6 + 3.0_rp) + x22 * (5.0_rp * x134 - 3.0_rp) * (x134 + 2.0_rp * x136 - 1.0_rp) + 10.0_rp * x32 * (-24.0_rp * ai * x134 - 18.0_rp * x13&
&7 + x191 * x35 + x191 * x36 + x191 * x38 * x4 + x191 * x39 * x6 - 24.0_rp * x2 * x26 + x40) + x33 * (x2 * x23 - 1.0_rp) * (4.0_rp&
& * ai * x134 + 3.0_rp * x137 + x2 * x81 + x30) + x55 * (-x134 * x49 - x137 * x48 + x191 * x4 * x50 + x191 * x45 + x191 * x47 + &
&x191 * x51 * x6 + x191 * x52 - x2 * x43 - x2 * x44 + x54))
        end do
    end subroutine T5_damp_erf
    pure subroutine T6_damp_erf(x, y, z, a, T)
        real(rp), intent(in) :: x(:), y(size(x)), z(size(x))
        real(rp), intent(in) :: a(:)
        real(rp), intent(inout) :: T(size(x), 84)
        integer(ip) :: i
        real(rp) :: xi, yi, zi
        real(rp) :: ai
        real(rp) :: x0
        real(rp) :: x1
        real(rp) :: x2
        real(rp) :: x3
        real(rp) :: x4
        real(rp) :: x5
        real(rp) :: x6
        real(rp) :: x7
        real(rp) :: x8
        real(rp) :: x10
        real(rp) :: x11
        real(rp) :: x12
        real(rp) :: x13
        real(rp) :: x14
        real(rp) :: x15
        real(rp) :: x17
        real(rp) :: x18
        real(rp) :: x19
        real(rp) :: x21
        real(rp) :: x22
        real(rp) :: x23
        real(rp) :: x24
        real(rp) :: x26
        real(rp) :: x27
        real(rp) :: x28
        real(rp) :: x29
        real(rp) :: x30
        real(rp) :: x32
        real(rp) :: x34
        real(rp) :: x35
        real(rp) :: x36
        real(rp) :: x37
        real(rp) :: x39
        real(rp) :: x40
        real(rp) :: x41
        real(rp) :: x42
        real(rp) :: x44
        real(rp) :: x46
        real(rp) :: x47
        real(rp) :: x48
        real(rp) :: x49
        real(rp) :: x50
        real(rp) :: x52
        real(rp) :: x61
        real(rp) :: x62
        real(rp) :: x64
        real(rp) :: x66
        real(rp) :: x67
        real(rp) :: x68
        real(rp) :: x69
        real(rp) :: x70
        real(rp) :: x71
        real(rp) :: x72
        real(rp) :: x73
        real(rp) :: x74
        real(rp) :: x75
        real(rp) :: x76
        real(rp) :: x77
        real(rp) :: x78
        real(rp) :: x79
        real(rp) :: x81
        real(rp) :: x82
        real(rp) :: x84
        real(rp) :: x88
        real(rp) :: x89
        real(rp) :: x92
        real(rp) :: x93
        real(rp) :: x95
        real(rp) :: x96
        real(rp) :: x100
        real(rp) :: x101
        real(rp) :: x107
        real(rp) :: x112
        real(rp) :: x113
        real(rp) :: x114
        real(rp) :: x121
        real(rp) :: x123
        real(rp) :: x125
        real(rp) :: x126
        real(rp) :: x129
        real(rp) :: x130
        real(rp) :: x134
        real(rp) :: x135
        real(rp) :: x136
        real(rp) :: x140
        real(rp) :: x141
        real(rp) :: x142
        real(rp) :: x144
        real(rp) :: x146
        real(rp) :: x148
        real(rp) :: x149
        real(rp) :: x150
        real(rp) :: x151
        real(rp) :: x152
        real(rp) :: x153
        real(rp) :: x156
        real(rp) :: x157
        real(rp) :: x159
        real(rp) :: x161
        real(rp) :: x164
        real(rp) :: x165
        real(rp) :: x166
        real(rp) :: x168
        real(rp) :: x169
        real(rp) :: x170
        real(rp) :: x171
        real(rp) :: x173
        real(rp) :: x174
        real(rp) :: x175
        real(rp) :: x176
        real(rp) :: x177
        real(rp) :: x178
        real(rp) :: x179
        real(rp) :: x181
        real(rp) :: x182
        real(rp) :: x183
        real(rp) :: x184
        real(rp) :: x185
        real(rp) :: x186
        real(rp) :: x187
        real(rp) :: x188
        real(rp) :: x189
        real(rp) :: x190
        real(rp) :: x191
        real(rp) :: x192
        real(rp) :: x193
        real(rp) :: x194
        real(rp) :: x195
        real(rp) :: x197
        real(rp) :: x200
        real(rp) :: x203
        real(rp) :: x205
        real(rp) :: x206
        real(rp) :: x207
        real(rp) :: x209
        real(rp) :: x210
        real(rp) :: x211
        real(rp) :: x213
        real(rp) :: x214
        real(rp) :: x216
        real(rp) :: x217
        real(rp) :: x218
        real(rp) :: x219
        real(rp) :: x220
        real(rp) :: x221
        real(rp) :: x222
        real(rp) :: x223
        real(rp) :: x224
        real(rp) :: x225
        real(rp) :: x226
        real(rp) :: x230
        real(rp) :: x231
        real(rp) :: x232
        real(rp) :: x234
        real(rp) :: x237
        real(rp) :: x239
        real(rp) :: x241
        real(rp) :: x242
        real(rp) :: x243
        real(rp) :: x244
        real(rp) :: x245
        real(rp) :: x247
        real(rp) :: x249
        real(rp) :: x252
        real(rp) :: x253
        real(rp) :: x254
        real(rp) :: x256
        real(rp) :: x258
        real(rp) :: x262
        real(rp) :: x268
        real(rp) :: x270
        real(rp) :: x278
        real(rp) :: x280
        real(rp) :: x281
        real(rp) :: x283
        real(rp) :: x284
        real(rp) :: x285
        real(rp) :: x286
        real(rp) :: x287
        real(rp) :: x295
        real(rp) :: x298
        real(rp) :: x300
        real(rp) :: x301
        real(rp) :: x302
        real(rp) :: x303
        real(rp) :: x306
        real(rp) :: x309
        real(rp) :: x310
        real(rp) :: x311
        real(rp) :: x314
        real(rp) :: x315
        real(rp) :: x316
        real(rp) :: x317
        real(rp) :: x319
        real(rp) :: x321
        real(rp) :: x326
        real(rp) :: x330
        real(rp) :: x332
        real(rp) :: x342
        real(rp) :: x351
        real(rp) :: x352
        real(rp) :: x362
        real(rp) :: x366
        real(rp) :: x367
        real(rp) :: x371
        real(rp) :: x374
        real(rp) :: x376
        real(rp) :: x377
        real(rp) :: x380
        real(rp) :: x381
        real(rp) :: x384
        real(rp) :: x386
        real(rp) :: x389
        real(rp) :: x390
        real(rp) :: x391
        real(rp) :: x392
        real(rp) :: x393
        real(rp) :: x396
        real(rp) :: x398
        real(rp) :: x403
        real(rp) :: x405
        real(rp) :: x406
        real(rp) :: x407
        real(rp) :: x409
        real(rp) :: x411
        real(rp) :: x413
        real(rp) :: x417
        real(rp) :: x432
        real(rp) :: x438
        real(rp) :: x444
        real(rp) :: x450
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1_rp / x3
            x5 = x0 * x4
            x6 = x3**(-2)
            x7 = xi**4
            x8 = x6 * x7
            x10 = x3**(-4)
            x11 = sqrt(ai)
            x12 = 0.564189583547756_rp * exp(-ai * x3)
            x13 = x11 * x12
            x14 = x10 * x13
            x15 = x0 * x14
            x17 = x3**(-3)
            x18 = 231.0_rp * x17
            x19 = erf(x11 * sqrt(x3))
            x21 = 45.0_rp * x19 * x3**(-3.5_rp)
            x22 = -30.0_rp * x5 + 35.0_rp * x8 + 3.0_rp
            x23 = ai * x0
            x24 = 2.0_rp * x23 + x5 - 1.0_rp
            x26 = x11 * x12 * x17
            x27 = 90.0_rp * x26
            x28 = 5.0_rp * x5 - 3.0_rp
            x29 = x0 * x6
            x30 = ai**2
            x32 = ai * x5
            x34 = 3.0_rp * x4
            x35 = -6.0_rp * ai - x34
            x36 = 4.0_rp * x0 * x30 + 3.0_rp * x29 + 4.0_rp * x32 + x35
            x37 = x0 * x17
            x39 = 3.0_rp * x5 - 1.0_rp
            x40 = ai**3
            x41 = x40 * x7
            x42 = x17 * x7
            x44 = x30 * x4
            x46 = 6.0_rp * ai + x34
            x47 = 18.0_rp * ai * x8 - 24.0_rp * x0 * x30 - 18.0_rp * x29 - 24.0_rp * x32 + 8.0_rp * x41 + 15.0_rp * x42 + 12.0_rp * x&
&44 * x7 + x46
            x48 = x12 * x6
            x49 = x11 * x48
            x50 = 30.0_rp * x49
            x52 = ai**4
            x61 = 60.0_rp * x30
            x62 = ai * x4
            x64 = 45.0_rp * x6 + x61 + 60.0_rp * x62
            x66 = 12.0_rp * x13
            x67 = 32.0_rp * ai**5
            x68 = x3**(-5)
            x69 = 945.0_rp * x68
            x70 = 80.0_rp * x4 * x52
            x71 = 240.0_rp * x40 * x6
            x72 = 600.0_rp * x17 * x30
            x73 = 1050.0_rp * ai * x10
            x74 = -45.0_rp * x6 - x61 - 60.0_rp * x62
            x75 = x12 * x4
            x76 = 2.0_rp * x11 * x75
            x77 = 9.0_rp * x5
            x78 = x19 * x3**(-5.5_rp)
            x79 = x0 * x78
            x81 = x13 * x68
            x82 = x0 * x81
            x84 = 7.0_rp * x5
            x88 = x19 * x3**(-4.5_rp)
            x89 = 105.0_rp * x88
            x92 = ai**1.5_rp
            x93 = x12 * x17 * x92
            x95 = x14 * x28
            x96 = x28 * x93
            x100 = 4.0_rp * x23 + 6.0_rp * x5 - 3.0_rp
            x101 = x14 * x39
            x107 = 4.0_rp * x30
            x112 = x26 * x39
            x113 = x48 * x92
            x114 = x113 * x39
            x121 = x75 * x92
            x123 = x1 * x4
            x125 = -12.0_rp * x123
            x126 = x1 * x29
            x129 = 60.0_rp * x0 * x88
            x130 = x22 * x89
            x134 = ai * x1
            x135 = x134 * x5
            x136 = x123 * x30
            x140 = 12.0_rp * x26
            x141 = 600.0_rp * x79 * (x84 - 3.0_rp)
            x142 = x13 / x3**6
            x144 = 4128.0_rp * x142 * x7
            x146 = x12 * x92
            x148 = 960.0_rp * x146 * x68 * x7
            x149 = 66.0_rp * x14 * x22
            x150 = 48.0_rp * x14 * (x107 * x7 - 8.0_rp * x23 - 12.0_rp * x5 + 12.0_rp * x62 * x7 + 15.0_rp * x8 + 1.0_rp)
            x151 = 48.0_rp * x107 * x7 - 384.0_rp * x23 - 576.0_rp * x5 + 576.0_rp * x62 * x7 + 720.0_rp * x8 + 48.0_rp
            x152 = x1 * x17
            x153 = x146 * x152
            x156 = x125 + 36.0_rp * x126 + 16.0_rp * x135 - 4.0_rp * x23 - 6.0_rp * x5 + 3.0_rp
            x157 = 16.0_rp * x14
            x159 = 16.0_rp * x47
            x161 = x1 * x6
            x164 = ai**2.5_rp
            x165 = x12 * x164
            x166 = 8.0_rp * x165
            x168 = 1152.0_rp * x28 * x82
            x169 = 864.0_rp * x24 * x82
            x170 = x10 * x146
            x171 = x0 * x170
            x173 = 576.0_rp * x171 * x28
            x174 = 288.0_rp * x39
            x175 = x174 * x82
            x176 = 288.0_rp * x24
            x177 = x101 * x176
            x178 = x171 * x176
            x179 = 240.0_rp * x13
            x181 = x0 * x179 * x68 * (x84 - 3.0_rp)
            x182 = 192.0_rp * x15 * x36
            x183 = 192.0_rp * x24
            x184 = x153 * x39
            x185 = 128.0_rp * x100 * x82
            x186 = x146 * x37
            x187 = 128.0_rp * x186 * x36
            x188 = 96.0_rp * x171 * x39
            x189 = x165 * x37
            x190 = 96.0_rp * x189 * x28
            x191 = 64.0_rp * x100 * x171
            x192 = x161 * x165
            x193 = 48.0_rp * x24
            x194 = x193 * x39
            x195 = 32.0_rp * x36
            x197 = 15.0_rp * x19 * x3**(-3.5_rp)
            x200 = 16.0_rp * x36
            x203 = 256.0_rp * x23 + 64.0_rp * x77 - 192.0_rp
            x205 = x164 * x48
            x206 = x165 * x29
            x207 = 48.0_rp * x14
            x209 = 8.0_rp * x164 * x75
            x210 = yi * zi
            x211 = x2 * x4
            x213 = x2 * x29
            x214 = 4.0_rp * x211
            x216 = ai * x2
            x217 = x216 * x5
            x218 = x211 * x30
            x219 = x17 * x2
            x220 = x146 * x219
            x221 = 3.0_rp - 12.0_rp * x211
            x222 = 36.0_rp * x213 + 16.0_rp * x217 + x221 - 4.0_rp * x23 - 6.0_rp * x5
            x223 = x13 * x219
            x224 = x2 * x6
            x225 = x220 * x39
            x226 = x165 * x224
            x230 = 8.0_rp * x134
            x231 = -4.0_rp * x23 - x77 + 3.0_rp
            x232 = 288.0_rp * x93
            x234 = 315.0_rp * x88
            x237 = 48.0_rp * x36
            x239 = 24.0_rp * x164 * x75
            x241 = x19 * x3**(-6.5_rp)
            x242 = x0 * x241
            x243 = 5670.0_rp * x1
            x244 = x1 * x78
            x245 = 945.0_rp * x244
            x247 = ai**3.5_rp * x12
            x249 = x152 * x165
            x252 = x161 * x247
            x253 = 48.0_rp * x39
            x254 = x1 * x81
            x256 = x1 * x14
            x258 = x1 * x170
            x262 = 576.0_rp * x10 * x165
            x268 = 9612.0_rp * x142
            x270 = xi * yi
            x278 = xi * zi
            x280 = x2 * x81
            x281 = x170 * x2
            x283 = x2 * x78
            x284 = 945.0_rp * x283
            x285 = x165 * x219
            x286 = x224 * x247
            x287 = x14 * x2
            x295 = x1 * x88
            x298 = yi**4
            x300 = 7560.0_rp * x242
            x301 = 16656.0_rp * x142
            x302 = x0 * x301
            x303 = x298 * x68
            x306 = x10 * x298
            x309 = 1122.0_rp * x81
            x310 = 768.0_rp * x24
            x311 = x310 * x81
            x314 = 640.0_rp * x247 * x37
            x315 = 492.0_rp * x39
            x316 = x17 * x298
            x317 = 384.0_rp * x24
            x319 = x298 * x6
            x321 = 128.0_rp * x24
            x326 = ai**4.5_rp * x12
            x330 = 16.0_rp * x39
            x332 = 96.0_rp * x247
            x342 = x123 * x326
            x351 = 3264.0_rp * x10 * x165
            x352 = 10080.0_rp * x0 * x146 * x68
            x362 = x1 * x2
            x366 = x170 * x362
            x367 = x161 * x2
            x371 = x211 * x326
            x374 = zi**4
            x376 = x374 * x68
            x377 = x10 * x374
            x380 = x17 * x374
            x381 = x374 * x6
            x384 = 64.0_rp * x374
            x386 = 10395.0_rp * x241
            x389 = 352.0_rp * x326
            x390 = 1584.0_rp * x247
            x391 = 5544.0_rp * x165
            x392 = 13860.0_rp * x146
            x393 = 20790.0_rp * x142
            x396 = 48.0_rp * ai**3.5_rp * x75 + 630.0_rp * x14 + 168.0_rp * x205 - x234 + 420.0_rp * x93
            x398 = x362 * x68
            x403 = 180.0_rp * x14
            x405 = -30.0_rp * x123 + 35.0_rp * x319 + 3.0_rp
            x406 = x123 + 2.0_rp * x134 - 1.0_rp
            x407 = 5.0_rp * x123 - 3.0_rp
            x409 = 4.0_rp * ai * x123 + x1 * x107 + 3.0_rp * x161 + x35
            x411 = x1 * x34 - 1.0_rp
            x413 = x298 * x40
            x417 = -24.0_rp * ai * x123 + 18.0_rp * ai * x319 - 24.0_rp * x1 * x30 - 18.0_rp * x161 + 12.0_rp * x298 * x44 + 15.0_rp&
                  & * x316 + 8.0_rp * x413 + x46
            x432 = 6.0_rp * x123 + 4.0_rp * x134 - 3.0_rp
            x438 = x113 * x411
            x444 = 48.0_rp * x406
            x450 = x123 * x216
            T(i, 1) = -120.0_rp * x13 * x28 * x36 * x37 - 180.0_rp * x15 * (-70.0_rp * x5 + 63.0_rp * x8 + 15.0_rp) + x21 * (x18 * &
&xi**6 + 105.0_rp * x5 - 315.0_rp * x8 - 5.0_rp) - x22 * x24 * x27 - x29 * x66 * (120.0_rp * ai * x42 - 80.0_rp * x0 * x40 + 105.0_rp&
& * x10 * x7 - 180.0_rp * x23 * x6 - 120.0_rp * x30 * x5 + 72.0_rp * x30 * x8 - 150.0_rp * x37 + 32.0_rp * x4 * x41 + 16.0_rp * x52 *&
& x7 + x64) - x39 * x47 * x50 - x76 * (-1800.0_rp * ai * x42 + 360.0_rp * x0 * x40 - 1575.0_rp * x10 * x7 + 810.0_rp * x23 * x6 + 5&
&40.0_rp * x30 * x5 - 1080.0_rp * x30 * x8 + 675.0_rp * x37 - 480.0_rp * x4 * x41 - 240.0_rp * x52 * x7 + x67 * xi**6 + x69 * xi**6 &
&+ x70 * xi**6 + x71 * xi**6 + x72 * xi**6 + x73 * xi**6 + x74)
            T(i, 2) = -xi * yi * (40.0_rp * x100 * x101 + 80.0_rp * x112 * x36 + 20.0_rp * x113 * x47 + 40.0_rp * x114 * x36 + 4.0&
&_rp * x121 * (120.0_rp * ai * x42 - 80.0_rp * x0 * x40 + 105.0_rp * x10 * x7 - 180.0_rp * x23 * x6 - 120.0_rp * x30 * x5 + 72.0_rp * &
&x30 * x8 - 150.0_rp * x37 + 32.0_rp * x4 * x41 + 16.0_rp * x52 * x7 + x64) + 180.0_rp * x14 * x22 + 30.0_rp * x14 * (-70.0_rp * x5 +&
& 63.0_rp * x8 + 15.0_rp) + 60.0_rp * x14 * (x107 * x7 - 8.0_rp * x23 - 12.0_rp * x5 + 12.0_rp * x62 * x7 + 15.0_rp * x8 + 1.0_rp) + 12&
&0.0_rp * x15 * x36 + 60.0_rp * x22 * x93 + 600.0_rp * x24 * x82 + 360.0_rp * x24 * x95 + 120.0_rp * x24 * x96 + 40.0_rp * x26 * x47 &
&+ 8.0_rp * x26 * (180.0_rp * ai * x8 + 30.0_rp * ai - x0 * x61 - 225.0_rp * x29 - 180.0_rp * x32 + 45.0_rp * x4 + 16.0_rp * x41 + 210&
&.0_rp * x42 + 72.0_rp * x44 * x7) + 120.0_rp * x28 * x82 + 4.0_rp * x49 * (120.0_rp * ai * x42 - 80.0_rp * x0 * x40 + 105.0_rp * x10 &
&* x7 - 180.0_rp * x23 * x6 - 120.0_rp * x30 * x5 + 72.0_rp * x30 * x8 - 150.0_rp * x37 + 32.0_rp * x4 * x41 + 16.0_rp * x52 * x7 + x&
&64) - 420.0_rp * x79 * (x77 - 5.0_rp) + 600.0_rp * x82 * (x84 - 3.0_rp) - x89 * (-70.0_rp * x5 + 63.0_rp * x8 + 15.0_rp))
            T(i, 3) = -xi * zi * (40.0_rp * x100 * x101 + 80.0_rp * x112 * x36 + 20.0_rp * x113 * x47 + 40.0_rp * x114 * x36 + 4.0&
&_rp * x121 * (120.0_rp * ai * x42 - 80.0_rp * x0 * x40 + 105.0_rp * x10 * x7 - 180.0_rp * x23 * x6 - 120.0_rp * x30 * x5 + 72.0_rp * &
&x30 * x8 - 150.0_rp * x37 + 32.0_rp * x4 * x41 + 16.0_rp * x52 * x7 + x64) + 180.0_rp * x14 * x22 + 30.0_rp * x14 * (-70.0_rp * x5 +&
& 63.0_rp * x8 + 15.0_rp) + 60.0_rp * x14 * (x107 * x7 - 8.0_rp * x23 - 12.0_rp * x5 + 12.0_rp * x62 * x7 + 15.0_rp * x8 + 1.0_rp) + 12&
&0.0_rp * x15 * x36 + 60.0_rp * x22 * x93 + 600.0_rp * x24 * x82 + 360.0_rp * x24 * x95 + 120.0_rp * x24 * x96 + 40.0_rp * x26 * x47 &
&+ 8.0_rp * x26 * (180.0_rp * ai * x8 + 30.0_rp * ai - x0 * x61 - 225.0_rp * x29 - 180.0_rp * x32 + 45.0_rp * x4 + 16.0_rp * x41 + 210&
&.0_rp * x42 + 72.0_rp * x44 * x7) + 120.0_rp * x28 * x82 + 4.0_rp * x49 * (120.0_rp * ai * x42 - 80.0_rp * x0 * x40 + 105.0_rp * x10 &
&* x7 - 180.0_rp * x23 * x6 - 120.0_rp * x30 * x5 + 72.0_rp * x30 * x8 - 150.0_rp * x37 + 32.0_rp * x4 * x41 + 16.0_rp * x52 * x7 + x&
&64) - 420.0_rp * x79 * (x77 - 5.0_rp) + 600.0_rp * x82 * (x84 - 3.0_rp) - x89 * (-70.0_rp * x5 + 63.0_rp * x8 + 15.0_rp))
            T(i, 4) = -x0 * x156 * x157 + 144.0_rp * x0 * x95 + x1 * x130 + x1 * x141 - x1 * x144 - x1 * x148 - x1 * x149 - x1 &
&* x150 - x1 * x168 - x1 * x169 - x1 * x173 - x1 * x175 - x1 * x177 - x1 * x178 - x1 * x181 - x1 * x182 - x1 * x185 - x1 * x187&
& - x1 * x188 - x1 * x190 - x1 * x191 + 24.0_rp * x114 * x24 + 4.0_rp * x121 * x47 - x123 * x166 * x47 - x126 * x165 * x195 + x12&
&9 * (x125 + 42.0_rp * x126 - x84 + 3.0_rp) - x13 * x152 * x159 + x13 * x195 * x37 - x140 * (120.0_rp * x1 * x42 - x107 * x7 + 4.0&
&_rp * x123 - 72.0_rp * x126 + 72.0_rp * x134 * x8 - 32.0_rp * x135 + 16.0_rp * x136 * x7 + 8.0_rp * x23 + 12.0_rp * x5 - 12.0_rp * x62&
& * x7 - 15.0_rp * x8 - 1.0_rp) - x146 * x159 * x161 + x146 * x200 * x29 + 72.0_rp * x15 * x24 + 24.0_rp * x15 * x39 - x151 * x153 &
&- 12.0_rp * x153 * x22 + x179 * x68 * x7 - x183 * x184 + 48.0_rp * x186 * x28 - x192 * x194 + x193 * x26 * x39 - x197 * x22 + 6.&
&0_rp * x22 * x26 + 4.0_rp * x47 * x49
            T(i, 5) = -x210 * (16.0_rp * x113 * x47 - x130 - x141 + x144 + x148 + x149 + x150 + x151 * x93 + x168 + x169 + x173&
& + x175 + x177 + x178 + x181 + x182 + x183 * x39 * x93 + x185 + x187 + x188 + x190 + x191 + x194 * x205 + x195 * x206 + x203 *&
& x82 + x207 * (x107 * x7 - 8.0_rp * x23 - 18.0_rp * x5 + 18.0_rp * x62 * x7 + 30.0_rp * x8 + 1.0_rp) + x209 * x47 + 12.0_rp * x22 * &
&x93 + 16.0_rp * x26 * x47 - 360.0_rp * x79 * (x84 - 2.0_rp))
            T(i, 6) = -x0 * x157 * x222 + 144.0_rp * x0 * x95 + 24.0_rp * x114 * x24 + 4.0_rp * x121 * x47 + x129 * (-12.0_rp * x2&
&11 + 42.0_rp * x213 - x84 + 3.0_rp) + x13 * x195 * x37 + x130 * x2 - x140 * (-x107 * x7 + 120.0_rp * x2 * x42 - 72.0_rp * x213 + x&
&214 + 72.0_rp * x216 * x8 - 32.0_rp * x217 + 16.0_rp * x218 * x7 + 8.0_rp * x23 + 12.0_rp * x5 - 12.0_rp * x62 * x7 - 15.0_rp * x8 - &
&1.0_rp) + x141 * x2 - x144 * x2 - x146 * x159 * x224 + x146 * x200 * x29 - x148 * x2 - x149 * x2 + 72.0_rp * x15 * x24 + 24.0_rp &
&* x15 * x39 - x150 * x2 - x151 * x220 - x159 * x223 - x165 * x195 * x213 - x166 * x211 * x47 - x168 * x2 - x169 * x2 - x173 * &
&x2 - x175 * x2 - x177 * x2 - x178 * x2 + x179 * x68 * x7 - x181 * x2 - x182 * x2 - x183 * x225 - x185 * x2 + 48.0_rp * x186 * x&
&28 - x187 * x2 - x188 * x2 - x190 * x2 - x191 * x2 + x193 * x26 * x39 - x194 * x226 - x197 * x22 - 12.0_rp * x22 * x220 + 6.0_rp&
& * x22 * x26 + 4.0_rp * x47 * x49
            T(i, 7) = -x270 * (3816.0_rp * x0 * x1 * x146 * x68 + x0 * x1 * x262 + x0 * x1 * x268 - 24.0_rp * x100 * x14 + 48.0_rp * x100 * x249 + 96.0_rp * x100 * x254 + 96.0_rp * x100 * x258 - 24.0_rp * x100 * x93 - x113 * x237 + x123 * x200 * x247 + 24.0_rp * x14 * x156 - 432.0_rp * x14 * x24 - 432.0_rp * x14 * x39 + 96.0_rp * x153 * x36 + 24.0_rp * x156 * x93 - 576.0_rp * x171 + x174&
& * x249 + x176 * x249 + x192 * x237 + x193 * x252 - 72.0_rp * x205 * x24 - 72.0_rp * x205 * x39 + x207 * (-6.0_rp * x123 + 24.0_rp&
& * x126 + x230 * x5 + x231) - x232 * x24 - x232 * x39 + x234 * x28 - x237 * x26 - x239 * x36 + 1152.0_rp * x24 * x254 + 864.0_rp&
& * x24 * x258 - x242 * x243 - x245 * x28 + 24.0_rp * x249 * x28 + x252 * x253 + 738.0_rp * x254 * x28 + 1152.0_rp * x254 * x39 + &
&96.0_rp * x256 * x36 + 204.0_rp * x258 * x28 + 864.0_rp * x258 * x39 + 1260.0_rp * x79 - 2088.0_rp * x82 - 198.0_rp * x95 - 36.0_rp *&
& x96)
            T(i, 8) = -x278 * (3816.0_rp * x0 * x1 * x146 * x68 + x0 * x1 * x262 + x0 * x1 * x268 - 8.0_rp * x100 * x14 + 48.0_rp&
& * x100 * x249 + 96.0_rp * x100 * x254 + 96.0_rp * x100 * x258 - 8.0_rp * x100 * x93 - x113 * x200 + x123 * x200 * x247 + 8.0_rp *&
& x14 * x156 - 144.0_rp * x14 * x24 - 144.0_rp * x14 * x39 + 96.0_rp * x153 * x36 + 8.0_rp * x156 * x93 + x157 * (-18.0_rp * x123 + &
&72.0_rp * x126 + 24.0_rp * x135 + x231) - 192.0_rp * x171 + x174 * x249 + x176 * x249 + x192 * x237 + x193 * x252 - x200 * x26 + &
&x203 * x254 + x203 * x258 - 24.0_rp * x205 * x24 - 24.0_rp * x205 * x39 - x209 * x36 + 1152.0_rp * x24 * x254 + 864.0_rp * x24 * x&
&258 - 96.0_rp * x24 * x93 - x242 * x243 - x245 * x28 + 24.0_rp * x249 * x28 + x252 * x253 + 738.0_rp * x254 * x28 + 1152.0_rp * x2&
&54 * x39 + 96.0_rp * x256 * x36 + 204.0_rp * x258 * x28 + 864.0_rp * x258 * x39 + x28 * x89 - 96.0_rp * x39 * x93 + 420.0_rp * x79 &
&- 696.0_rp * x82 - 66.0_rp * x95 - 12.0_rp * x96)
            T(i, 9) = -x270 * (3816.0_rp * x0 * x146 * x2 * x68 + x0 * x2 * x262 + x0 * x2 * x268 - 8.0_rp * x100 * x14 + 96.0_rp&
& * x100 * x280 + 96.0_rp * x100 * x281 + 48.0_rp * x100 * x285 - 8.0_rp * x100 * x93 - x113 * x200 + 8.0_rp * x14 * x222 - 144.0_rp&
& * x14 * x24 - 144.0_rp * x14 * x39 + x157 * (-18.0_rp * x211 + 72.0_rp * x213 + 24.0_rp * x217 + x231) - 192.0_rp * x171 + x174 * &
&x285 + x176 * x285 + x193 * x286 - 5670.0_rp * x2 * x242 + x200 * x211 * x247 - x200 * x26 + x203 * x280 + x203 * x281 - 24.0_rp&
& * x205 * x24 - 24.0_rp * x205 * x39 - x209 * x36 + 96.0_rp * x220 * x36 + 8.0_rp * x222 * x93 + x226 * x237 + 1152.0_rp * x24 * x&
&280 + 864.0_rp * x24 * x281 - 96.0_rp * x24 * x93 + x253 * x286 + 738.0_rp * x28 * x280 + 204.0_rp * x28 * x281 - x28 * x284 + 24.&
&0_rp * x28 * x285 + x28 * x89 + 1152.0_rp * x280 * x39 + 864.0_rp * x281 * x39 + 96.0_rp * x287 * x36 - 96.0_rp * x39 * x93 + 420.0&
&_rp * x79 - 696.0_rp * x82 - 66.0_rp * x95 - 12.0_rp * x96)
            T(i, 10) = -x278 * (3816.0_rp * x0 * x146 * x2 * x68 + x0 * x2 * x262 + x0 * x2 * x268 - 24.0_rp * x100 * x14 + 96.0&
&_rp * x100 * x280 + 96.0_rp * x100 * x281 + 48.0_rp * x100 * x285 - 24.0_rp * x100 * x93 - x113 * x237 + 24.0_rp * x14 * x222 - 432&
&.0_rp * x14 * x24 - 432.0_rp * x14 * x39 - 576.0_rp * x171 + x174 * x285 + x176 * x285 + x193 * x286 - 5670.0_rp * x2 * x242 + x20&
&0 * x211 * x247 - 72.0_rp * x205 * x24 - 72.0_rp * x205 * x39 + x207 * (-6.0_rp * x211 + 24.0_rp * x213 + 8.0_rp * x217 + x231) + 9&
&6.0_rp * x220 * x36 + 24.0_rp * x222 * x93 + x226 * x237 - x232 * x24 - x232 * x39 + x234 * x28 - x237 * x26 - x239 * x36 + 1152&
&.0_rp * x24 * x280 + 864.0_rp * x24 * x281 + x253 * x286 + 738.0_rp * x28 * x280 + 204.0_rp * x28 * x281 - x28 * x284 + 24.0_rp * x&
&28 * x285 + 1152.0_rp * x280 * x39 + 864.0_rp * x281 * x39 + 96.0_rp * x287 * x36 + 1260.0_rp * x79 - 2088.0_rp * x82 - 198.0_rp * x&
&95 - 36.0_rp * x96)
            T(i, 11) = -32.0_rp * ai**4.5_rp * x24 * x298 * x75 - 10080.0_rp * x0 * x146 * x303 - 3264.0_rp * x0 * x165 * x306 + 1&
&80.0_rp * x0 * x88 + 684.0_rp * x1 * x101 + 5040.0_rp * x1 * x171 + 1440.0_rp * x1 * x189 - 3780.0_rp * x1 * x79 + 8712.0_rp * x1 * &
&x82 - 42.0_rp * x112 - x113 * x193 - 12.0_rp * x114 + x123 * x24 * x332 + 192.0_rp * x126 * x247 - x146 * x306 * x310 - x146 * x3&
&06 * x315 - 456.0_rp * x15 + 576.0_rp * x153 * x24 - x165 * x316 * x317 - 120.0_rp * x165 * x316 * x39 + x176 * x192 + 264.0_rp * &
&x184 - 240.0_rp * x186 + 48.0_rp * x192 * x39 - x193 * x26 - 48.0_rp * x206 + x21 * x39 - x239 * x24 + 576.0_rp * x24 * x256 - x24&
&7 * x319 * x321 - x247 * x319 * x330 - 64.0_rp * x29 * x298 * x326 - 630.0_rp * x295 * x39 + x298 * x300 - x298 * x302 - x298 * &
&x309 * x39 - x298 * x311 - x298 * x314 + 945.0_rp * x298 * x39 * x78
            T(i, 12) = -x210 * (-48.0_rp * ai**3.5_rp * x24 * x75 + x0 * x1 * x301 + x0 * x1 * x351 - x1 * x300 + x1 * x314 + x1&
& * x352 - 342.0_rp * x101 + 64.0_rp * x126 * x326 - x14 * x176 - 2520.0_rp * x171 - 720.0_rp * x189 - 144.0_rp * x205 * x24 - 24.0_rp * x205 * x39 - x232 * x24 + x234 * x39 + 32.0_rp * x24 * x342 - x245 * x39 + x249 * x317 + 120.0_rp * x249 * x39 + x252 * x321&
& + x252 * x330 + x254 * x310 + 1122.0_rp * x254 * x39 + x258 * x310 + x258 * x315 - x29 * x332 - 132.0_rp * x39 * x93 + 1890.0_rp&
& * x79 - 4356.0_rp * x82)
            T(i, 13) = -x0 * x351 * x362 + 114.0_rp * x1 * x101 + 840.0_rp * x1 * x171 + 240.0_rp * x1 * x189 + 7560.0_rp * x1 * x&
&2 * x242 - x1 * x39 * x89 - 630.0_rp * x1 * x79 + 1452.0_rp * x1 * x82 + 114.0_rp * x101 * x2 - 14.0_rp * x112 - 16.0_rp * x113 * x&
&24 - 4.0_rp * x114 + 16.0_rp * x123 * x24 * x247 - 64.0_rp * x126 * x2 * x326 + 32.0_rp * x126 * x247 + x129 - 152.0_rp * x15 + 96.&
&0_rp * x153 * x24 + x166 * x224 * x39 + 840.0_rp * x171 * x2 + 44.0_rp * x184 - 80.0_rp * x186 + 240.0_rp * x189 * x2 + x192 * x193&
& + 8.0_rp * x192 * x39 + x193 * x226 + x197 * x39 - 32.0_rp * x2 * x24 * x342 + x2 * x245 * x39 - x2 * x249 * x317 - 120.0_rp * x&
&2 * x249 * x39 - x2 * x39 * x89 - 630.0_rp * x2 * x79 + 1452.0_rp * x2 * x82 - 16.0_rp * x206 - x209 * x24 + 16.0_rp * x211 * x24 &
&* x247 + 32.0_rp * x213 * x247 + 96.0_rp * x220 * x24 + 44.0_rp * x225 + 96.0_rp * x24 * x256 - 16.0_rp * x24 * x26 + 96.0_rp * x24 &
&* x287 - x247 * x321 * x367 - x247 * x330 * x367 - x302 * x362 - x309 * x362 * x39 - x310 * x366 - x311 * x362 - x314 * x362 -&
& x315 * x366 - x352 * x362
            T(i, 14) = -x210 * (-48.0_rp * ai**3.5_rp * x24 * x75 + x0 * x2 * x301 + x0 * x2 * x351 - 342.0_rp * x101 - x14 * x17&
&6 - 2520.0_rp * x171 - 720.0_rp * x189 - 7560.0_rp * x2 * x242 + x2 * x314 + x2 * x352 - 144.0_rp * x205 * x24 - 24.0_rp * x205 * x&
&39 + 64.0_rp * x213 * x326 - x232 * x24 + x234 * x39 + 32.0_rp * x24 * x371 + x280 * x310 + 1122.0_rp * x280 * x39 + x281 * x310 &
&+ x281 * x315 - x284 * x39 + x285 * x317 + 120.0_rp * x285 * x39 + x286 * x321 + x286 * x330 - x29 * x332 - 132.0_rp * x39 * x93&
& + 1890.0_rp * x79 - 4356.0_rp * x82)
            T(i, 15) = -32.0_rp * ai**4.5_rp * x24 * x374 * x75 - 10080.0_rp * x0 * x146 * x376 - 3264.0_rp * x0 * x165 * x377 + 1&
&80.0_rp * x0 * x88 + 684.0_rp * x101 * x2 - 42.0_rp * x112 - x113 * x193 - 12.0_rp * x114 - x146 * x310 * x377 - x146 * x315 * x37&
&7 - 456.0_rp * x15 - x165 * x317 * x380 - 120.0_rp * x165 * x380 * x39 + 5040.0_rp * x171 * x2 + x176 * x226 - 240.0_rp * x186 + 1&
&440.0_rp * x189 * x2 - x193 * x26 - 630.0_rp * x2 * x39 * x88 - 3780.0_rp * x2 * x79 + 8712.0_rp * x2 * x82 - 48.0_rp * x206 + x21 &
&* x39 + x211 * x24 * x332 + 192.0_rp * x213 * x247 + 576.0_rp * x220 * x24 + 264.0_rp * x225 + x226 * x253 - x239 * x24 + 576.0_rp&
& * x24 * x287 - x247 * x321 * x381 - x247 * x330 * x381 - x29 * x326 * x384 + x300 * x374 - x302 * x374 - x309 * x374 * x39 - &
&x311 * x374 - x314 * x374 + 945.0_rp * x374 * x39 * x78
            T(i, 16) = -x270 * (240.0_rp * ai**3.5_rp * x75 + 64.0_rp * ai**5.5_rp * x298 * x75 + 3150.0_rp * x14 + 840.0_rp * x205 &
&+ 9450.0_rp * x244 - 5040.0_rp * x249 - 1440.0_rp * x252 - 18900.0_rp * x254 - 12600.0_rp * x258 - x298 * x386 + x298 * x393 + x303&
& * x392 + x306 * x391 + x316 * x390 + x319 * x389 - 320.0_rp * x342 - 1575.0_rp * x88 + 2100.0_rp * x93)
            T(i, 17) = -x278 * (64.0_rp * ai**5.5_rp * x298 * x75 - 7560.0_rp * x1 * x170 + 5670.0_rp * x244 - 3024.0_rp * x249 - 8&
&64.0_rp * x252 - 11340.0_rp * x254 - x298 * x386 + x298 * x393 + x303 * x392 + x306 * x391 + x316 * x390 + x319 * x389 - 192.0_rp&
& * x342 + x396)
            T(i, 18) = -x270 * (64.0_rp * ai**5.5_rp * x12 * x123 * x2 + x10 * x362 * x391 + x152 * x2 * x390 - 3780.0_rp * x170 &
&* x2 + x245 - 504.0_rp * x249 - 144.0_rp * x252 - 1890.0_rp * x254 - 1260.0_rp * x258 - 5670.0_rp * x280 + 2835.0_rp * x283 - 1512.0&
&_rp * x285 - 432.0_rp * x286 - 32.0_rp * x342 - x362 * x386 + x362 * x393 + x367 * x389 - 96.0_rp * x371 + x392 * x398 + x396)
            T(i, 19) = -x278 * (64.0_rp * ai**5.5_rp * x12 * x123 * x2 + x10 * x362 * x391 + x152 * x2 * x390 - x243 * x81 + 283&
&5.0_rp * x244 - 1512.0_rp * x249 - 432.0_rp * x252 - 3780.0_rp * x258 - 1890.0_rp * x280 - 1260.0_rp * x281 + x284 - 504.0_rp * x285 &
&- 144.0_rp * x286 - 96.0_rp * x342 - x362 * x386 + x362 * x393 + x367 * x389 - 32.0_rp * x371 + x392 * x398 + x396)
            T(i, 20) = -x270 * (ai**5.5_rp * x384 * x75 - 11340.0_rp * x280 - 7560.0_rp * x281 + 5670.0_rp * x283 - 3024.0_rp * x28&
&5 - 864.0_rp * x286 - 192.0_rp * x371 - x374 * x386 + x374 * x393 + x376 * x392 + x377 * x391 + x380 * x390 + x381 * x389 + x396&
&)
            T(i, 21) = -x278 * (240.0_rp * ai**3.5_rp * x75 + ai**5.5_rp * x384 * x75 + 3150.0_rp * x14 + 840.0_rp * x205 - 18900.0&
&_rp * x280 - 12600.0_rp * x281 + 9450.0_rp * x283 - 5040.0_rp * x285 - 1440.0_rp * x286 - 320.0_rp * x371 - x374 * x386 + x374 * x39&
&3 + x376 * x392 + x377 * x391 + x380 * x390 + x381 * x389 - 1575.0_rp * x88 + 2100.0_rp * x93)
            T(i, 22) = -x1 * x403 * (-70.0_rp * x123 + 63.0_rp * x319 + 15.0_rp) - 120.0_rp * x13 * x152 * x407 * x409 - x161 * x6&
&6 * (120.0_rp * ai * x316 - 80.0_rp * x1 * x40 - 180.0_rp * x134 * x6 - 120.0_rp * x136 - 150.0_rp * x152 + 16.0_rp * x298 * x52 + 7&
&2.0_rp * x30 * x319 + 105.0_rp * x306 + 32.0_rp * x4 * x413 + x64) + x21 * (105.0_rp * x123 + x18 * yi**6 - 315.0_rp * x319 - 5.0_rp&
&) - x27 * x405 * x406 - x411 * x417 * x50 - x76 * (-1800.0_rp * ai * x316 + 360.0_rp * x1 * x40 + 810.0_rp * x134 * x6 + 540.0_rp &
&* x136 + 675.0_rp * x152 - 240.0_rp * x298 * x52 - 1080.0_rp * x30 * x319 - 1575.0_rp * x306 - 480.0_rp * x4 * x413 + x67 * yi**6 +&
& x69 * yi**6 + x70 * yi**6 + x71 * yi**6 + x72 * yi**6 + x73 * yi**6 + x74)
            T(i, 23) = -x210 * (20.0_rp * x113 * x417 + 4.0_rp * x121 * (120.0_rp * ai * x316 - 80.0_rp * x1 * x40 - 180.0_rp * x13&
&4 * x6 - 120.0_rp * x136 - 150.0_rp * x152 + 16.0_rp * x298 * x52 + 72.0_rp * x30 * x319 + 105.0_rp * x306 + 32.0_rp * x4 * x413 + x&
&64) + 360.0_rp * x14 * x406 * x407 + 40.0_rp * x14 * x411 * x432 + 30.0_rp * x14 * (-70.0_rp * x123 + 63.0_rp * x319 + 15.0_rp) + 60&
&.0_rp * x14 * (x107 * x298 + x125 - x230 + 12.0_rp * x298 * x62 + 15.0_rp * x319 + 1.0_rp) - 420.0_rp * x244 * (9.0_rp * x123 - 5.0_rp) + 600.0_rp * x254 * x406 + 120.0_rp * x254 * x407 + 600.0_rp * x254 * (7.0_rp * x123 - 3.0_rp) + 120.0_rp * x256 * x409 + 80.0_rp &
&* x26 * x409 * x411 + 40.0_rp * x26 * x417 + 8.0_rp * x26 * (-180.0_rp * ai * x123 + 180.0_rp * ai * x319 + 30.0_rp * ai - x1 * x61&
& - 225.0_rp * x161 + 72.0_rp * x298 * x44 + 210.0_rp * x316 + 45.0_rp * x4 + 16.0_rp * x413) + x403 * x405 + 60.0_rp * x405 * x93 + &
&120.0_rp * x406 * x407 * x93 + 40.0_rp * x409 * x438 + 4.0_rp * x49 * (120.0_rp * ai * x316 - 80.0_rp * x1 * x40 - 180.0_rp * x134 *&
& x6 - 120.0_rp * x136 - 150.0_rp * x152 + 16.0_rp * x298 * x52 + 72.0_rp * x30 * x319 + 105.0_rp * x306 + 32.0_rp * x4 * x413 + x64)&
& - x89 * (-70.0_rp * x123 + 63.0_rp * x319 + 15.0_rp))
            T(i, 24) = 144.0_rp * x1 * x14 * x407 - x1 * x157 * (-6.0_rp * x123 - 4.0_rp * x134 + x221 + 36.0_rp * x367 + 16.0_rp *&
& x450) + 4.0_rp * x121 * x417 + 32.0_rp * x13 * x152 * x409 - 66.0_rp * x14 * x2 * x405 - 192.0_rp * x14 * x362 * x409 - x140 * (-&
&x107 * x298 + 12.0_rp * x123 + 120.0_rp * x2 * x316 + x214 + 72.0_rp * x216 * x319 + 16.0_rp * x218 * x298 + x230 - 12.0_rp * x298 &
&* x62 - 15.0_rp * x319 - 72.0_rp * x367 - 32.0_rp * x450 - 1.0_rp) - 4128.0_rp * x142 * x2 * x298 + 16.0_rp * x146 * x161 * x409 - 9&
&60.0_rp * x146 * x2 * x303 - 16.0_rp * x146 * x224 * x417 - 128.0_rp * x153 * x2 * x409 + 48.0_rp * x153 * x407 - 32.0_rp * x165 * &
&x367 * x409 - x166 * x211 * x417 + x179 * x303 - x179 * x398 * (7.0_rp * x123 - 3.0_rp) - x197 * x405 - x2 * x207 * (x107 * x298&
& + x125 - x230 + 12.0_rp * x298 * x62 + 15.0_rp * x319 + 1.0_rp) + 600.0_rp * x2 * x244 * (7.0_rp * x123 - 3.0_rp) - 96.0_rp * x2 * x&
&249 * x407 + x2 * x405 * x89 - 12.0_rp * x220 * x405 - 192.0_rp * x220 * x406 * x411 - 48.0_rp * x220 * (x107 * x298 + x125 - x23&
&0 + 12.0_rp * x298 * x62 + 15.0_rp * x319 + 1.0_rp) - 16.0_rp * x223 * x417 - x226 * x411 * x444 + 72.0_rp * x256 * x406 + 24.0_rp *&
& x256 * x411 + 6.0_rp * x26 * x405 + x26 * x411 * x444 - 288.0_rp * x287 * x406 * x411 + 60.0_rp * x295 * (-7.0_rp * x123 + x221 +&
& 42.0_rp * x367) - 864.0_rp * x362 * x406 * x81 - 1152.0_rp * x362 * x407 * x81 - 288.0_rp * x362 * x411 * x81 - 128.0_rp * x362 * &
&x432 * x81 - 288.0_rp * x366 * x406 - 576.0_rp * x366 * x407 - 96.0_rp * x366 * x411 - 64.0_rp * x366 * x432 + 24.0_rp * x406 * x43&
&8 + 4.0_rp * x417 * x49
            T(i, 25) = -x210 * (-48.0_rp * x113 * x409 - 432.0_rp * x14 * x406 - 198.0_rp * x14 * x407 - 432.0_rp * x14 * x411 - 2&
&4.0_rp * x14 * x432 + 24.0_rp * x14 * (-6.0_rp * x123 - 4.0_rp * x134 + x221 + 36.0_rp * x367 + 16.0_rp * x450) + 3816.0_rp * x146 * &
&x398 - x2 * x241 * x243 - 72.0_rp * x205 * x406 - 72.0_rp * x205 * x411 + x207 * (-9.0_rp * x123 - 4.0_rp * x134 - 6.0_rp * x211 + &
&24.0_rp * x367 + 8.0_rp * x450 + 3.0_rp) + 16.0_rp * x211 * x247 * x409 + 96.0_rp * x220 * x409 + 48.0_rp * x226 * x409 - x232 * x40&
&6 - x232 * x411 + x234 * x407 - x239 * x409 + 1260.0_rp * x244 - 2088.0_rp * x254 - 576.0_rp * x258 - 48.0_rp * x26 * x409 + x262 &
&* x362 + x268 * x362 + 1152.0_rp * x280 * x406 + 738.0_rp * x280 * x407 + 1152.0_rp * x280 * x411 + 96.0_rp * x280 * x432 + 864.0_rp * x281 * x406 + 204.0_rp * x281 * x407 + 864.0_rp * x281 * x411 + 96.0_rp * x281 * x432 - x284 * x407 + 288.0_rp * x285 * x406 +&
& 24.0_rp * x285 * x407 + 288.0_rp * x285 * x411 + 48.0_rp * x285 * x432 + 48.0_rp * x286 * x411 + x286 * x444 + 96.0_rp * x287 * x4&
&09 - 36.0_rp * x407 * x93 - 24.0_rp * x432 * x93 + 24.0_rp * x93 * (-6.0_rp * x123 - 4.0_rp * x134 + x221 + 36.0_rp * x367 + 16.0_rp &
&* x450))
            T(i, 26) = -32.0_rp * ai**4.5_rp * x374 * x406 * x75 - 10080.0_rp * x1 * x146 * x376 - 3264.0_rp * x1 * x165 * x377 + &
&7560.0_rp * x1 * x241 * x374 - x1 * x301 * x374 - x113 * x444 - 768.0_rp * x146 * x377 * x406 - 492.0_rp * x146 * x377 * x411 - 6&
&40.0_rp * x152 * x247 * x374 - 240.0_rp * x153 - x161 * x326 * x384 - 384.0_rp * x165 * x380 * x406 - 120.0_rp * x165 * x380 * x41&
&1 - 48.0_rp * x192 - 3780.0_rp * x2 * x244 + 1440.0_rp * x2 * x249 - 630.0_rp * x2 * x411 * x88 + x21 * x411 + x211 * x332 * x406 &
&+ 576.0_rp * x220 * x406 + 264.0_rp * x220 * x411 + 288.0_rp * x226 * x406 + 48.0_rp * x226 * x411 - x239 * x406 + 192.0_rp * x247 &
&* x367 - 128.0_rp * x247 * x381 * x406 - 16.0_rp * x247 * x381 * x411 - 456.0_rp * x256 - 42.0_rp * x26 * x411 - x26 * x444 + 576.&
&0_rp * x287 * x406 + 684.0_rp * x287 * x411 + 180.0_rp * x295 - x309 * x374 * x411 + 8712.0_rp * x362 * x81 + 5040.0_rp * x366 - 76&
&8.0_rp * x374 * x406 * x81 + 945.0_rp * x374 * x411 * x78 - 12.0_rp * x438
            T(i, 27) = -x210 * (240.0_rp * ai**3.5_rp * x75 + ai**5.5_rp * x384 * x75 + 3150.0_rp * x14 + 840.0_rp * x205 - 18900.0&
&_rp * x280 - 12600.0_rp * x281 + 9450.0_rp * x283 - 5040.0_rp * x285 - 1440.0_rp * x286 - 320.0_rp * x371 - x374 * x386 + x374 * x39&
&3 + x376 * x392 + x377 * x391 + x380 * x390 + x381 * x389 - 1575.0_rp * x88 + 2100.0_rp * x93)
            T(i, 28) = -x2 * x403 * (-70.0_rp * x211 + 63.0_rp * x381 + 15.0_rp) + x21 * (x18 * zi**6 + 105.0_rp * x211 - 315.0_rp &
&* x381 - 5.0_rp) - 120.0_rp * x223 * (5.0_rp * x211 - 3.0_rp) * (ai * x214 + x107 * x2 + 3.0_rp * x224 + x35) - x224 * x66 * (120.0&
&_rp * ai * x380 - 80.0_rp * x2 * x40 - 180.0_rp * x216 * x6 - 120.0_rp * x218 - 150.0_rp * x219 + 72.0_rp * x30 * x381 + 32.0_rp * x3&
&74 * x4 * x40 + 16.0_rp * x374 * x52 + 105.0_rp * x377 + x64) - x27 * (-30.0_rp * x211 + 35.0_rp * x381 + 3.0_rp) * (x211 + 2.0_rp *&
& x216 - 1.0_rp) - x50 * (x2 * x34 - 1.0_rp) * (-24.0_rp * ai * x211 + 18.0_rp * ai * x381 - 24.0_rp * x2 * x30 - 18.0_rp * x224 + 8.&
&0_rp * x374 * x40 + 12.0_rp * x374 * x44 + 15.0_rp * x380 + x46) - x76 * (-1800.0_rp * ai * x380 + 360.0_rp * x2 * x40 + 810.0_rp * &
&x216 * x6 + 540.0_rp * x218 + 675.0_rp * x219 - 1080.0_rp * x30 * x381 - 480.0_rp * x374 * x4 * x40 - 240.0_rp * x374 * x52 - 1575.&
&0_rp * x377 + x67 * zi**6 + x69 * zi**6 + x70 * zi**6 + x71 * zi**6 + x72 * zi**6 + x73 * zi**6 + x74)
        end do
    end subroutine T6_damp_erf
end module T_tensor_damp_erf
