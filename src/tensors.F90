module T_tensor
    use tensors_recursive, only: Tn_recursive
    implicit none
    private

    public Tn
contains
    subroutine Tn(n, x, y, z, T)
        integer, intent(in) :: n
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), (n + 1)*(n + 2)*(n + 3)/6)
        select case (n)
        case (0)
            call T0(x, y, z, T(:, 1:1))
        case (1)
            call T0(x, y, z, T(:, 1:1))
            call T1(x, y, z, T(:, 2:4))
        case (2)
            call T0(x, y, z, T(:, 1:1))
            call T1(x, y, z, T(:, 2:4))
            call T2(x, y, z, T(:, 5:10))
        case (3)
            call T0(x, y, z, T(:, 1:1))
            call T1(x, y, z, T(:, 2:4))
            call T2(x, y, z, T(:, 5:10))
            call T3(x, y, z, T(:, 11:20))
        case (4)
            call T0(x, y, z, T(:, 1:1))
            call T1(x, y, z, T(:, 2:4))
            call T2(x, y, z, T(:, 5:10))
            call T3(x, y, z, T(:, 11:20))
            call T4(x, y, z, T(:, 21:35))
        case (5)
            call T0(x, y, z, T(:, 1:1))
            call T1(x, y, z, T(:, 2:4))
            call T2(x, y, z, T(:, 5:10))
            call T3(x, y, z, T(:, 11:20))
            call T4(x, y, z, T(:, 21:35))
            call T5(x, y, z, T(:, 36:56))
        case (6)
            call T0(x, y, z, T(:, 1:1))
            call T1(x, y, z, T(:, 2:4))
            call T2(x, y, z, T(:, 5:10))
            call T3(x, y, z, T(:, 11:20))
            call T4(x, y, z, T(:, 21:35))
            call T5(x, y, z, T(:, 36:56))
            call T6(x, y, z, T(:, 57:84))
        case (7)
            call T0(x, y, z, T(:, 1:1))
            call T1(x, y, z, T(:, 2:4))
            call T2(x, y, z, T(:, 5:10))
            call T3(x, y, z, T(:, 11:20))
            call T4(x, y, z, T(:, 21:35))
            call T5(x, y, z, T(:, 36:56))
            call T6(x, y, z, T(:, 57:84))
            call T7(x, y, z, T(:, 85:120))
        case (8)
            call T0(x, y, z, T(:, 1:1))
            call T1(x, y, z, T(:, 2:4))
            call T2(x, y, z, T(:, 5:10))
            call T3(x, y, z, T(:, 11:20))
            call T4(x, y, z, T(:, 21:35))
            call T5(x, y, z, T(:, 36:56))
            call T6(x, y, z, T(:, 57:84))
            call T7(x, y, z, T(:, 85:120))
            call T8(x, y, z, T(:, 121:165))
        case (9)
            call T0(x, y, z, T(:, 1:1))
            call T1(x, y, z, T(:, 2:4))
            call T2(x, y, z, T(:, 5:10))
            call T3(x, y, z, T(:, 11:20))
            call T4(x, y, z, T(:, 21:35))
            call T5(x, y, z, T(:, 36:56))
            call T6(x, y, z, T(:, 57:84))
            call T7(x, y, z, T(:, 85:120))
            call T8(x, y, z, T(:, 121:165))
            call T9(x, y, z, T(:, 166:220))
        case (10)
            call T0(x, y, z, T(:, 1:1))
            call T1(x, y, z, T(:, 2:4))
            call T2(x, y, z, T(:, 5:10))
            call T3(x, y, z, T(:, 11:20))
            call T4(x, y, z, T(:, 21:35))
            call T5(x, y, z, T(:, 36:56))
            call T6(x, y, z, T(:, 57:84))
            call T7(x, y, z, T(:, 85:120))
            call T8(x, y, z, T(:, 121:165))
            call T9(x, y, z, T(:, 166:220))
            call T10(x, y, z, T(:, 221:286))
        case (11)
            call T0(x, y, z, T(:, 1:1))
            call T1(x, y, z, T(:, 2:4))
            call T2(x, y, z, T(:, 5:10))
            call T3(x, y, z, T(:, 11:20))
            call T4(x, y, z, T(:, 21:35))
            call T5(x, y, z, T(:, 36:56))
            call T6(x, y, z, T(:, 57:84))
            call T7(x, y, z, T(:, 85:120))
            call T8(x, y, z, T(:, 121:165))
            call T9(x, y, z, T(:, 166:220))
            call T10(x, y, z, T(:, 221:286))
            call T11(x, y, z, T(:, 287:364))
        case (12)
            call T0(x, y, z, T(:, 1:1))
            call T1(x, y, z, T(:, 2:4))
            call T2(x, y, z, T(:, 5:10))
            call T3(x, y, z, T(:, 11:20))
            call T4(x, y, z, T(:, 21:35))
            call T5(x, y, z, T(:, 36:56))
            call T6(x, y, z, T(:, 57:84))
            call T7(x, y, z, T(:, 85:120))
            call T8(x, y, z, T(:, 121:165))
            call T9(x, y, z, T(:, 166:220))
            call T10(x, y, z, T(:, 221:286))
            call T11(x, y, z, T(:, 287:364))
            call T12(x, y, z, T(:, 365:455))
        case default
            call Tn_recursive(n, x, y, z, T)
        end select
    end subroutine Tn
    subroutine T0(x, y, z, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), 1)
        integer :: i
        real(8) :: xi, yi, zi
        do concurrent(i=1:size(x))
            xi = x(i)
            yi = y(i)
            zi = z(i)
            T(i, 1) = (xi**2 + yi**2 + zi**2)**(-0.5d0)
        end do
    end subroutine T0
    subroutine T1(x, y, z, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), 4)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: x0
        do concurrent(i=1:size(x))
            xi = x(i)
            yi = y(i)
            zi = z(i)
            x0 = (xi**2 + yi**2 + zi**2)**(-1.5d0)
            T(i, 1) = -x0*xi
            T(i, 2) = -x0*yi
            T(i, 3) = -x0*zi
        end do
    end subroutine T1
    subroutine T2(x, y, z, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), 10)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x6
        real(8) :: x7
        do concurrent(i=1:size(x))
            xi = x(i)
            yi = y(i)
            zi = z(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = x3**(-1.5d0)
            x5 = 3.0d0/x3
            x6 = x3**(-2.5d0)
            x7 = 3.0d0*x6*xi
            T(i, 1) = x4*(x0*x5 - 1.0d0)
            T(i, 2) = x7*yi
            T(i, 3) = x7*zi
            T(i, 4) = x4*(x1*x5 - 1.0d0)
            T(i, 5) = 3.0d0*x6*yi*zi
            T(i, 6) = x4*(x2*x5 - 1.0d0)
        end do
    end subroutine T2
    subroutine T3(x, y, z, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), 20)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x6
        real(8) :: x7
        real(8) :: x8
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x13
        real(8) :: x14
        do concurrent(i=1:size(x))
            xi = x(i)
            yi = y(i)
            zi = z(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 5.0d0/x3
            x5 = x0*x4
            x6 = 3.0d0*x3**(-2.5d0)
            x7 = x6*xi
            x8 = x6*(x5 - 1.0d0)
            x9 = x1*x4
            x10 = x9 - 1.0d0
            x11 = x2*x4
            x12 = x11 - 1.0d0
            x13 = x6*yi
            x14 = x6*zi
            T(i, 1) = -x7*(x5 - 3.0d0)
            T(i, 2) = -x8*yi
            T(i, 3) = -x8*zi
            T(i, 4) = -x10*x7
            T(i, 5) = -15.0d0*x3**(-3.5d0)*xi*yi*zi
            T(i, 6) = -x12*x7
            T(i, 7) = -x13*(x9 - 3.0d0)
            T(i, 8) = -x10*x14
            T(i, 9) = -x12*x13
            T(i, 10) = -x14*(x11 - 3.0d0)
        end do
    end subroutine T3
    subroutine T4(x, y, z, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), 35)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x6
        real(8) :: x7
        real(8) :: x8
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x13
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
        do concurrent(i=1:size(x))
            xi = x(i)
            yi = y(i)
            zi = z(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0/x3
            x5 = x0*x4
            x6 = 35.0d0/x3**2
            x7 = 3.0d0*x3**(-2.5d0)
            x8 = x3**(-3.5d0)
            x9 = 15.0d0*x8*yi
            x10 = 7.0d0*x5
            x11 = xi*(x10 - 3.0d0)
            x12 = 15.0d0*x8*zi
            x13 = 5.0d0*x4
            x14 = -x1*x13
            x15 = x0*x6
            x16 = 1.0d0 - 5.0d0*x5
            x17 = x9*zi
            x18 = -x13*x2
            x19 = 7.0d0*x4
            x20 = x1*x19
            x21 = x20 - 3.0d0
            x22 = x9*xi
            x23 = x12*xi
            x24 = x19*x2
            x25 = x24 - 3.0d0
            x26 = 30.0d0*x4
            T(i, 1) = x7*(-30.0d0*x5 + x6*xi**4 + 3.0d0)
            T(i, 2) = x11*x9
            T(i, 3) = x11*x12
            T(i, 4) = x7*(x1*x15 + x14 + x16)
            T(i, 5) = x17*(x10 - 1.0d0)
            T(i, 6) = x7*(x15*x2 + x16 + x18)
            T(i, 7) = x21*x22
            T(i, 8) = x23*(x20 - 1.0d0)
            T(i, 9) = x22*(x24 - 1.0d0)
            T(i, 10) = x23*x25
            T(i, 11) = x7*(-x1*x26 + x6*yi**4 + 3.0d0)
            T(i, 12) = x17*x21
            T(i, 13) = x7*(x1*x2*x6 + x14 + x18 + 1.0d0)
            T(i, 14) = x17*x25
            T(i, 15) = x7*(-x2*x26 + x6*zi**4 + 3.0d0)
        end do
    end subroutine T4
    subroutine T5(x, y, z, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), 56)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x6
        real(8) :: x7
        real(8) :: x8
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x13
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
        real(8) :: x27
        real(8) :: x28
        real(8) :: x29
        real(8) :: x30
        real(8) :: x31
        real(8) :: x32
        real(8) :: x33
        real(8) :: x34
        real(8) :: x35
        real(8) :: x36
        real(8) :: x37
        real(8) :: x38
        real(8) :: x39
        real(8) :: x40
        real(8) :: x41
        real(8) :: x42
        do concurrent(i=1:size(x))
            xi = x(i)
            yi = y(i)
            zi = z(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0/x3
            x5 = x0*x4
            x6 = x3**(-2)
            x7 = x6*xi**4
            x8 = x3**(-3.5d0)
            x9 = 15.0d0*x8
            x10 = x9*xi
            x11 = 45.0d0*x8
            x12 = x11*(-14.0d0*x5 + 21.0d0*x7 + 1.0d0)
            x13 = 21.0d0*x4
            x14 = -x1*x13
            x15 = 63.0d0*x6
            x16 = x0*x15
            x17 = x1*x16
            x18 = -7.0d0*x5
            x19 = x18 + 3.0d0
            x20 = 315.0d0*x3**(-4.5d0)*xi*yi*zi
            x21 = -x13*x2
            x22 = x16*x2
            x23 = 7.0d0*x4
            x24 = -x1*x23
            x25 = x17 + x24
            x26 = 3.0d0 - 21.0d0*x5
            x27 = x9*yi
            x28 = x18 + 1.0d0
            x29 = x9*zi
            x30 = x2*x23
            x31 = -x30
            x32 = x22 + x31
            x33 = x1*x4
            x34 = yi**4
            x35 = 21.0d0*x6
            x36 = -14.0d0*x33 + x34*x35 + 1.0d0
            x37 = x11*xi
            x38 = x1*x2
            x39 = x2*x4
            x40 = zi**4
            x41 = x35*x40 - 14.0d0*x39 + 1.0d0
            x42 = x15*x38 + 3.0d0
            T(i, 1) = -x10*(-70.0d0*x5 + 63.0d0*x7 + 15.0d0)
            T(i, 2) = -x12*yi
            T(i, 3) = -x12*zi
            T(i, 4) = -x10*(x14 + x17 + x19)
            T(i, 5) = -x20*(3.0d0*x5 - 1.0d0)
            T(i, 6) = -x10*(x19 + x21 + x22)
            T(i, 7) = -x27*(x25 + x26)
            T(i, 8) = -x29*(x25 + x28)
            T(i, 9) = -x27*(x28 + x32)
            T(i, 10) = -x29*(x26 + x32)
            T(i, 11) = -x36*x37
            T(i, 12) = -x20*(3.0d0*x33 - 1.0d0)
            T(i, 13) = -x10*(2.0d0*x33*(4.0d0*x39 - 1.0d0) + 20.0d0*x38*x6 + (x30 - 1.0d0)*(5.0d0*x33 - 1.0d0))
            T(i, 14) = -x20*(3.0d0*x39 - 1.0d0)
            T(i, 15) = -x37*x41
            T(i, 16) = -x27*(x15*x34 - 70.0d0*x33 + 15.0d0)
            T(i, 17) = -x11*x36*zi
            T(i, 18) = -x27*(x21 + x24 + x42)
            T(i, 19) = -x29*(x14 + x31 + x42)
            T(i, 20) = -x11*x41*yi
            T(i, 21) = -x29*(x15*x40 - 70.0d0*x39 + 15.0d0)
        end do
    end subroutine T5
    subroutine T6(x, y, z, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), 84)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x6
        real(8) :: x7
        real(8) :: x8
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x13
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
        real(8) :: x27
        real(8) :: x28
        real(8) :: x29
        real(8) :: x30
        real(8) :: x31
        real(8) :: x32
        real(8) :: x33
        real(8) :: x34
        real(8) :: x35
        real(8) :: x36
        real(8) :: x37
        real(8) :: x38
        real(8) :: x39
        real(8) :: x40
        real(8) :: x41
        real(8) :: x42
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
        real(8) :: x53
        real(8) :: x55
        real(8) :: x56
        real(8) :: x57
        real(8) :: x58
        real(8) :: x59
        real(8) :: x60
        real(8) :: x61
        real(8) :: x62
        real(8) :: x63
        do concurrent(i=1:size(x))
            xi = x(i)
            yi = y(i)
            zi = z(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0/x3
            x5 = x0*x4
            x6 = xi**4
            x7 = x3**(-2)
            x8 = x6*x7
            x9 = x3**(-3)
            x10 = 231.0d0*x9
            x11 = x3**(-3.5d0)
            x12 = 45.0d0*x11
            x13 = 315.0d0*yi
            x14 = 33.0d0*x8
            x15 = x3**(-4.5d0)
            x16 = x15*xi
            x17 = x16*(x14 - 30.0d0*x5 + 5.0d0)
            x18 = 315.0d0*zi
            x19 = 7.0d0*x4
            x20 = x1*x19
            x21 = x10*x6
            x22 = 14.0d0*x5 - 21.0d0*x8
            x23 = x0*x1
            x24 = 126.0d0*x7
            x25 = -x23*x24 - 1.0d0
            x26 = x15*zi
            x27 = x13*x26
            x28 = x0*x2
            x29 = -x24*x28
            x30 = x19*x2 - 1.0d0
            x31 = x1*x4
            x32 = 3.0d0*x31
            x33 = -x32
            x34 = 11.0d0*x7
            x35 = 1.0d0 - 3.0d0*x5
            x36 = 945.0d0*x16
            x37 = 33.0d0*x7
            x38 = x23*x37
            x39 = x16*x18
            x40 = x2*x4
            x41 = 9.0d0*x40
            x42 = x28*x37
            x43 = x13*x16
            x44 = 3.0d0*x40
            x45 = -x44
            x46 = 7.0d0*x5
            x47 = yi**4
            x48 = x0*x10
            x49 = 21.0d0*x7
            x50 = x47*x49
            x51 = 14.0d0*x31
            x52 = -x50 + x51
            x53 = 1.0d0 - 9.0d0*x5
            x55 = x1*x2
            x56 = zi**4
            x57 = 14.0d0*x40 - x49*x56 - 1.0d0
            x58 = -30.0d0*x31 + x37*x47 + 5.0d0
            x59 = x31*(4.0d0*x40 - 1.0d0)
            x60 = x37*x56
            x61 = -30.0d0*x40 + x60 + 5.0d0
            x62 = 315.0d0*x7
            x63 = -x24*x55
            T(i, 1) = x12*(x10*xi**6 + 105.0d0*x5 - 315.0d0*x8 - 5.0d0)
            T(i, 2) = x13*x17
            T(i, 3) = x17*x18
            T(i, 4) = x12*(x1*x21 + x20 + x22 + x25)
            T(i, 5) = x27*(x14 - 18.0d0*x5 + 1.0d0)
            T(i, 6) = x12*(x2*x21 + x22 + x29 + x30)
            T(i, 7) = x36*yi*(x23*x34 + x33 + x35)
            T(i, 8) = x39*(-9.0d0*x31 + x35 + x38)
            T(i, 9) = x43*(x35 - x41 + x42)
            T(i, 10) = x36*zi*(x28*x34 + x35 + x45)
            T(i, 11) = x12*(x25 + x46 + x47*x48 + x52)
            T(i, 12) = x27*(x33 + x38 + x53)
            T(i, 13) = 15.0d0*x11*(693.0d0*x2*x23*x9 + x20 - 63.0d0*x23*x7 - 63.0d0*x28*x7 + x30 + x46 - 63.0d0*x55*x7)
            T(i, 14) = x27*(x42 + x45 + x53)
            T(i, 15) = x12*(x29 + x46 + x48*x56 + x57)
            T(i, 16) = x43*x58
            T(i, 17) = x39*(4.0d0*x31*(x32 - 1.0d0) + x50 - x51 + 1.0d0)
            T(i, 18) = 105.0d0*x16*yi*(28.0d0*x55*x7 + 2.0d0*x59 + (x20 - 3.0d0)*(x41 - 1.0d0))
            T(i, 19) = 45.0d0*x16*zi*(10.0d0*x30*x31 + 8.0d0*x31*(2.0d0*x40 - 1.0d0) + 10.0d0*x59 + 7.0d0*(5.0d0*x31 - 1.0d0)*(x44 - 1.0d0))
            T(i, 20) = x43*(-18.0d0*x40 + x60 + 1.0d0)
            T(i, 21) = x39*x61
            T(i, 22) = x12*(x10*yi**6 + 105.0d0*x31 - x47*x62 - 5.0d0)
            T(i, 23) = x27*x58
            T(i, 24) = x12*(x10*x2*x47 + x30 + x52 + x63)
            T(i, 25) = 945.0d0*x26*yi*(x33 + x34*x55 + x45 + 1.0d0)
            T(i, 26) = x12*(x1*x10*x56 + x20 + x57 + x63)
            T(i, 27) = x27*x61
            T(i, 28) = x12*(x10*zi**6 + 105.0d0*x40 - x56*x62 - 5.0d0)
        end do
    end subroutine T6
    subroutine T7(x, y, z, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), 120)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x6
        real(8) :: x7
        real(8) :: x8
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x13
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
        real(8) :: x27
        real(8) :: x28
        real(8) :: x29
        real(8) :: x30
        real(8) :: x31
        real(8) :: x32
        real(8) :: x33
        real(8) :: x34
        real(8) :: x35
        real(8) :: x36
        real(8) :: x37
        real(8) :: x38
        real(8) :: x39
        real(8) :: x40
        real(8) :: x41
        real(8) :: x42
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
        real(8) :: x53
        real(8) :: x54
        real(8) :: x55
        real(8) :: x56
        real(8) :: x57
        real(8) :: x58
        real(8) :: x59
        real(8) :: x60
        real(8) :: x61
        real(8) :: x62
        real(8) :: x63
        real(8) :: x64
        real(8) :: x65
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
        real(8) :: x80
        real(8) :: x81
        real(8) :: x82
        real(8) :: x83
        real(8) :: x84
        real(8) :: x85
        real(8) :: x86
        real(8) :: x87
        real(8) :: x88
        real(8) :: x89
        real(8) :: x90
        real(8) :: x91
        real(8) :: x92
        real(8) :: x93
        real(8) :: x94
        real(8) :: x95
        do concurrent(i=1:size(x))
            xi = x(i)
            yi = y(i)
            zi = z(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0/x3
            x5 = x0*x4
            x6 = xi**4
            x7 = x3**(-2)
            x8 = x6*x7
            x9 = x3**(-3)
            x10 = 429.0d0*x9
            x11 = x10*xi**6
            x12 = x3**(-4.5d0)
            x13 = 315.0d0*x12
            x14 = x13*xi
            x15 = x13*(x11 + 135.0d0*x5 - 495.0d0*x8 - 5.0d0)
            x16 = 30.0d0*x5
            x17 = 45.0d0*x4
            x18 = x1*x17
            x19 = -33.0d0*x8
            x20 = x10*x6
            x21 = x1*x20 + x19
            x22 = x0*x7
            x23 = 330.0d0*x22
            x24 = -x1*x23 - 5.0d0
            x25 = x3**(-5.5d0)*xi*yi*zi
            x26 = 945.0d0*x25
            x27 = x17*x2
            x28 = x19 + x2*x20
            x29 = -x2*x23 - 5.0d0
            x30 = x1*x22
            x31 = -66.0d0*x30
            x32 = 143.0d0*x9
            x33 = x32*x6
            x34 = x1*x4
            x35 = 3.0d0*x34
            x36 = x35 - 1.0d0
            x37 = 18.0d0*x5
            x38 = x19 + x37
            x39 = 945.0d0*x12
            x40 = x39*yi
            x41 = 9.0d0*x34
            x42 = -198.0d0*x30 - 1.0d0
            x43 = x13*zi
            x44 = x2*x22
            x45 = -198.0d0*x44
            x46 = x2*x4
            x47 = 9.0d0*x46 - 1.0d0
            x48 = x13*yi
            x49 = -66.0d0*x44
            x50 = 3.0d0*x46 - 1.0d0
            x51 = x39*zi
            x52 = yi**4
            x53 = x0*x32
            x54 = 3.0d0*x5
            x55 = x54 - 1.0d0
            x56 = 33.0d0*x7
            x57 = -x52*x56
            x58 = 18.0d0*x34 + x57
            x59 = x39*xi
            x60 = 9.0d0 - 33.0d0*x5
            x61 = x1*x2
            x62 = x61*x7
            x63 = x0*x56
            x64 = x0*x10
            x65 = x61*x64
            x66 = -x1*x63 + x47 + x65
            x67 = -x2*x63 + x41
            x68 = zi**4
            x69 = x56*x68
            x70 = -x69
            x71 = 18.0d0*x46
            x72 = x70 + x71
            x73 = 45.0d0*x5
            x74 = x52*x64
            x75 = 30.0d0*x34 + x57
            x76 = 9.0d0*x5
            x77 = -x56*x61 + x76
            x78 = x64*x68
            x79 = 30.0d0*x46 + x70
            x80 = x52*x7
            x81 = x10*yi**6
            x82 = 135.0d0*x34 - 495.0d0*x80 + x81 - 5.0d0
            x83 = 4.0d0*x34
            x84 = 14.0d0*x34
            x85 = 4.0d0*x46
            x86 = 2.0d0*x46 - 1.0d0
            x87 = 8.0d0*x34
            x88 = x85 - 1.0d0
            x89 = x68*x7
            x90 = x10*zi**6
            x91 = 135.0d0*x46 - 495.0d0*x89 + x90 - 5.0d0
            x92 = x2*x52
            x93 = -330.0d0*x62 - 5.0d0
            x94 = -66.0d0*x62
            x95 = x1*x68
            T(i, 1) = -x14*(x11 + 315.0d0*x5 - 693.0d0*x8 - 35.0d0)
            T(i, 2) = -x15*yi
            T(i, 3) = -x15*zi
            T(i, 4) = -x14*(x16 + x18 + x21 + x24)
            T(i, 5) = -x26*(-110.0d0*x5 + 143.0d0*x8 + 15.0d0)
            T(i, 6) = -x14*(x16 + x27 + x28 + x29)
            T(i, 7) = -x40*(x1*x33 + x31 + x36 + x38)
            T(i, 8) = -x43*(x21 + x37 + x41 + x42)
            T(i, 9) = -x48*(x28 + x37 + x45 + x47)
            T(i, 10) = -x51*(x2*x33 + x38 + x49 + x50)
            T(i, 11) = -x59*(x31 + x52*x53 + x55 + x58)
            T(i, 12) = -x26*(143.0d0*x30 - 33.0d0*x34 + x60)
            T(i, 13) = -x14*(x54 - 99.0d0*x62 + x66 + x67)
            T(i, 14) = -x26*(143.0d0*x44 - 33.0d0*x46 + x60)
            T(i, 15) = -x59*(x49 + x53*x68 + x55 + x72)
            T(i, 16) = -x48*(x24 + x73 + x74 + x75)
            T(i, 17) = -x43*(x42 + x58 + x74 + x76)
            T(i, 18) = -x48*(x35 - 99.0d0*x44 + x66 + x77)
            T(i, 19) = -x43*(-99.0d0*x30 + x50 + x65 + x67 + x77)
            T(i, 20) = -x48*(x45 + x72 + x76 + x78 - 1.0d0)
            T(i, 21) = -x43*(x29 + x73 + x78 + x79)
            T(i, 22) = -x14*x82
            T(i, 23) = -x26*(-90.0d0*x34 + 99.0d0*x80 + x83*(11.0d0*x34 - 5.0d0) + 15.0d0)
            T(i, 24) = -x14*(56.0d0*x36*x62 + x47*(21.0d0*x80 - x84 + 1.0d0) + x83*(-x35 + 18.0d0*x62 - x85 + 1.0d0))
            T(i, 25) = -315.0d0*x25*(x47*x84 + x84*x88 + x86*x87 + 3.0d0*(7.0d0*x34 - 3.0d0)*(11.0d0*x46 - 3.0d0))
            T(i, 26) = -45.0d0*x12*xi*(20.0d0*x34*x88*(7.0d0*x46 - 1.0d0) + 280.0d0*x50*x62 + 160.0d0*x62*x86 + x87*(-12.0d0*x46 + 16.0d0*x89 + 1.0d0) + 7.0d0*(5.0d0*x34 - 1.0d0)*(x69 - x71 + 1.0d0))
            T(i, 27) = -x26*(-110.0d0*x46 + 143.0d0*x89 + 15.0d0)
            T(i, 28) = -x14*x91
            T(i, 29) = -x48*(315.0d0*x34 - 693.0d0*x80 + x81 - 35.0d0)
            T(i, 30) = -x43*x82
            T(i, 31) = -x48*(x10*x92 + x27 + x75 + x93)
            T(i, 32) = -x51*(x32*x92 + x50 + x58 + x94)
            T(i, 33) = -x40*(x32*x95 + x36 + x72 + x94)
            T(i, 34) = -x43*(x10*x95 + x18 + x79 + x93)
            T(i, 35) = -x48*x91
            T(i, 36) = -x43*(315.0d0*x46 - 693.0d0*x89 + x90 - 35.0d0)
        end do
    end subroutine T7
    subroutine T8(x, y, z, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), 165)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x6
        real(8) :: x7
        real(8) :: x8
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x13
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
        real(8) :: x27
        real(8) :: x28
        real(8) :: x29
        real(8) :: x30
        real(8) :: x31
        real(8) :: x32
        real(8) :: x33
        real(8) :: x34
        real(8) :: x35
        real(8) :: x36
        real(8) :: x37
        real(8) :: x38
        real(8) :: x39
        real(8) :: x40
        real(8) :: x41
        real(8) :: x42
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
        real(8) :: x53
        real(8) :: x54
        real(8) :: x55
        real(8) :: x56
        real(8) :: x57
        real(8) :: x58
        real(8) :: x59
        real(8) :: x60
        real(8) :: x61
        real(8) :: x62
        real(8) :: x63
        real(8) :: x64
        real(8) :: x65
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
        real(8) :: x80
        real(8) :: x81
        real(8) :: x82
        real(8) :: x83
        real(8) :: x84
        real(8) :: x85
        real(8) :: x86
        real(8) :: x87
        real(8) :: x88
        real(8) :: x89
        real(8) :: x90
        real(8) :: x91
        real(8) :: x92
        real(8) :: x93
        real(8) :: x94
        real(8) :: x95
        real(8) :: x96
        real(8) :: x97
        real(8) :: x98
        real(8) :: x99
        real(8) :: x100
        real(8) :: x101
        real(8) :: x102
        real(8) :: x103
        real(8) :: x104
        real(8) :: x105
        real(8) :: x106
        real(8) :: x107
        real(8) :: x108
        real(8) :: x109
        real(8) :: x110
        real(8) :: x111
        real(8) :: x112
        real(8) :: x113
        real(8) :: x114
        real(8) :: x115
        real(8) :: x116
        real(8) :: x117
        real(8) :: x118
        real(8) :: x119
        real(8) :: x120
        real(8) :: x121
        real(8) :: x122
        real(8) :: x123
        real(8) :: x124
        real(8) :: x125
        real(8) :: x126
        real(8) :: x127
        real(8) :: x128
        real(8) :: x129
        real(8) :: x130
        real(8) :: x131
        real(8) :: x132
        real(8) :: x134
        real(8) :: x135
        real(8) :: x136
        real(8) :: x138
        real(8) :: x139
        real(8) :: x140
        real(8) :: x141
        real(8) :: x142
        real(8) :: x143
        real(8) :: x144
        real(8) :: x145
        real(8) :: x146
        do concurrent(i=1:size(x))
            xi = x(i)
            yi = y(i)
            zi = z(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0/x3
            x5 = x0*x4
            x6 = xi**4
            x7 = x3**(-2)
            x8 = x6*x7
            x9 = xi**6
            x10 = x3**(-3)
            x11 = x10*x9
            x12 = x3**(-4)
            x13 = 6435.0d0*x12
            x14 = x3**(-4.5d0)
            x15 = 315.0d0*x14
            x16 = x3**(-5.5d0)
            x17 = x16*yi
            x18 = x17*xi
            x19 = 2027025.0d0*x11 + 1091475.0d0*x5 - 2837835.0d0*x8 - 99225.0d0
            x20 = x16*xi*zi
            x21 = 45.0d0*x4
            x22 = -x1*x21
            x23 = x10*x6
            x24 = 6435.0d0*x23
            x25 = x13*x9
            x26 = x0*x7
            x27 = 1485.0d0*x26
            x28 = x1*x27 + 5.0d0
            x29 = -429.0d0*x11 - 135.0d0*x5 + 495.0d0*x8
            x30 = 33.0d0*x5
            x31 = -143.0d0*x8
            x32 = x17*zi
            x33 = -x2*x21
            x34 = x2*x27 + 5.0d0
            x35 = 1430.0d0*x26
            x36 = -x1*x35
            x37 = x36 - 45.0d0
            x38 = 165.0d0*x4
            x39 = x1*x38
            x40 = 2145.0d0*x23
            x41 = x1*x40 + x39
            x42 = 330.0d0*x5 - 429.0d0*x8
            x43 = 945.0d0*x18
            x44 = x36 - 15.0d0
            x45 = x31 + 110.0d0*x5
            x46 = 945.0d0*x20
            x47 = -x2*x35
            x48 = x47 - 15.0d0
            x49 = x2*x38
            x50 = x2*x40 + x49
            x51 = x47 - 45.0d0
            x52 = yi**4
            x53 = x0*x52
            x54 = 858.0d0*x10
            x55 = 858.0d0*x23
            x56 = x1*x26
            x57 = 2145.0d0*x52
            x58 = x12*x6
            x59 = 33.0d0*x7
            x60 = x52*x59
            x61 = x1*x4
            x62 = x60 - 18.0d0*x61
            x63 = -18.0d0*x5 + 33.0d0*x8
            x64 = x63 + 1.0d0
            x65 = 945.0d0*x14
            x66 = 11.0d0*x61
            x67 = 715.0d0*x10
            x68 = x1*x6
            x69 = x31 + 66.0d0*x5
            x70 = -286.0d0*x56 - 3.0d0
            x71 = 2835.0d0*x32
            x72 = 429.0d0*x10
            x73 = x2*x6
            x74 = 99.0d0*x2
            x75 = x1*x7
            x76 = x13*x2
            x77 = x2*x4
            x78 = 9.0d0*x77
            x79 = x10*x2
            x80 = x0*x1
            x81 = -2574.0d0*x79*x80
            x82 = 198.0d0*x56 - x78 + x81
            x83 = -9.0d0*x61
            x84 = x2*x26
            x85 = x83 + 198.0d0*x84
            x86 = -286.0d0*x84
            x87 = 11.0d0*x77
            x88 = x87 - 3.0d0
            x89 = zi**4
            x90 = x0*x89
            x91 = 2145.0d0*x89
            x92 = x59*x89 - 18.0d0*x77 + 1.0d0
            x93 = 165.0d0*x5
            x94 = 2145.0d0*x10
            x95 = x53*x94 + x93
            x96 = x52*x7
            x97 = 330.0d0*x61 - 429.0d0*x96
            x98 = 143.0d0*x96
            x99 = -x98
            x100 = 11.0d0*x5
            x101 = 2835.0d0*x20
            x102 = x2*x75
            x103 = x2*x80*x94
            x104 = -429.0d0*x102 + x103 + x30 - 9.0d0
            x105 = 33.0d0*x61 - 429.0d0*x84
            x106 = 33.0d0*x77
            x107 = x106 - 429.0d0*x56
            x108 = 66.0d0*x77
            x109 = x7*x89
            x110 = 143.0d0*x109
            x111 = -x110
            x112 = 2835.0d0*x18
            x113 = x90*x94 + x93
            x114 = -429.0d0*x109 + 330.0d0*x77
            x115 = -45.0d0*x5
            x116 = 6435.0d0*x10
            x117 = yi**6
            x118 = x0*x13
            x119 = x117*x72
            x120 = 135.0d0*x61
            x121 = 495.0d0*x96
            x122 = -x119 - x120 + x121
            x123 = 110.0d0*x61
            x124 = 945.0d0*x32
            x125 = x2*x52
            x126 = 198.0d0*x102 - 9.0d0*x5
            x127 = x1*x89
            x128 = x1*x13
            x129 = 110.0d0*x77
            x130 = zi**6
            x131 = 495.0d0*x109 - x130*x72 - 135.0d0*x77
            x132 = x117*x67 + 385.0d0*x61 - 1001.0d0*x96 - 35.0d0
            x134 = 28.0d0*x61
            x135 = x134*(x78 - 1.0d0)
            x136 = 4.0d0*x77
            x138 = 8.0d0*x61
            x139 = 2.0d0*x77 - 1.0d0
            x140 = 16.0d0*x109
            x141 = x140 - 12.0d0*x77 + 1.0d0
            x142 = x136 - 1.0d0
            x143 = x10*x130
            x144 = -1001.0d0*x109 + x130*x67 + 385.0d0*x77 - 35.0d0
            x145 = 1485.0d0*x102 + 5.0d0
            x146 = -1430.0d0*x102 - 45.0d0
            T(i, 1) = x15*(-12012.0d0*x11 + x13*xi**8 - 1260.0d0*x5 + 6930.0d0*x8 + 35.0d0)
            T(i, 2) = x18*x19
            T(i, 3) = x19*x20
            T(i, 4) = x15*(-x1*x24 + x1*x25 + x22 + x28 + x29)
            T(i, 5) = 14175.0d0*x32*(143.0d0*x11 + x30 + x31 - 1.0d0)
            T(i, 6) = x15*(-x2*x24 + x2*x25 + x29 + x33 + x34)
            T(i, 7) = x43*(x37 + x41 + x42)
            T(i, 8) = x46*(x41 + x44 + x45)
            T(i, 9) = x43*(x45 + x48 + x50)
            T(i, 10) = x46*(x42 + x50 + x51)
            T(i, 11) = x65*(-x1*x55 - x53*x54 + 396.0d0*x56 + x57*x58 + x62 + x64)
            T(i, 12) = x71*(x66 + x67*x68 + x69 + x70)
            T(i, 13) = x15*(x64 - x68*x72 + x68*x76 - x72*x73 + x74*x75 + x82 + x85)
            T(i, 14) = x71*(x67*x73 + x69 + x86 + x88)
            T(i, 15) = x65*(-x2*x55 - x54*x90 + x58*x91 + x63 + 396.0d0*x84 + x92)
            T(i, 16) = x43*(x37 + x95 + x97)
            T(i, 17) = x101*(x100 + x53*x67 + 66.0d0*x61 + x70 + x99)
            T(i, 18) = x43*(x104 + x105 - 143.0d0*x56 + 99.0d0*x77)
            T(i, 19) = x46*(x104 + x107 + 99.0d0*x61 - 143.0d0*x84)
            T(i, 20) = x112*(x100 + x108 + x111 + x67*x90 + x86 - 3.0d0)
            T(i, 21) = x46*(x113 + x114 + x51)
            T(i, 22) = x15*(x115 - x116*x53 + x117*x118 + x122 + x28)
            T(i, 23) = x124*(x123 + x44 + x95 + x99)
            T(i, 24) = x15*(-x125*x72 + x126 + x26*x74 - x53*x72 + x53*x76 + x62 + x82 + 1.0d0)
            T(i, 25) = x124*(-143.0d0*x102 + x103 + x105 + x107 + 99.0d0*x5 - 9.0d0)
            T(i, 26) = x15*(x126 - x127*x72 + x128*x90 + 99.0d0*x56 - x72*x90 + x81 + x85 + x92)
            T(i, 27) = x124*(x111 + x113 + x129 + x48)
            T(i, 28) = x15*(x115 - x116*x90 + x118*x130 + x131 + x34)
            T(i, 29) = x112*x132
            T(i, 30) = x101*(x119 + x120 - x121 + 2.0d0*x61*(-x123 + x98 + 15.0d0) - 5.0d0)
            T(i, 31) = x43*(72.0d0*x102*(x66 - 5.0d0) + 4.0d0*x61*(66.0d0*x102 - x66 - 20.0d0*x77 + 5.0d0) + 3.0d0*(x87 - 1.0d0)*(x60 - 30.0d0*x61 + 5.0d0))
            T(i, 32) = x46*(x134*(18.0d0*x102 - x136 - 3.0d0*x61 + 1.0d0) + x135*(3.0d0*x61 - 1.0d0) + x138*(24.0d0*x102 - x136 + x83 + 2.0d0) + 3.0d0*x88*(-14.0d0*x61 + 21.0d0*x96 + 1.0d0))
            T(i, 33) = 315.0d0*x18*(224.0d0*x102*x139 + 168.0d0*x102*x88 + x135*x142 + x138*x141 + 3.0d0*(7.0d0*x61 - 3.0d0)*(-x108 + x110 + 3.0d0))
            T(i, 34) = 45.0d0*x20*(400.0d0*x139*x61*(7.0d0*x77 - 1.0d0) + 200.0d0*x141*x61 + 700.0d0*x142*x61*(3.0d0*x77 - 1.0d0) + 350.0d0*x61*x92 + 80.0d0*x61*(x140 - 16.0d0*x77 + 3.0d0) + 21.0d0*(5.0d0*x61 - 1.0d0)*(x110 - x129 + 15.0d0))
            T(i, 35) = 14175.0d0*x18*(x106 + x111 + 143.0d0*x143 - 1.0d0)
            T(i, 36) = x101*x144
            T(i, 37) = x15*(-12012.0d0*x10*x117 + x13*yi**8 - 1260.0d0*x61 + 6930.0d0*x96 + 35.0d0)
            T(i, 38) = x132*x71
            T(i, 39) = x15*(-x116*x125 + x117*x76 + x122 + x145 + x33)
            T(i, 40) = x124*(x146 + x49 + x57*x79 + x97)
            T(i, 41) = x65*(396.0d0*x102 + x12*x57*x89 - x125*x54 - x127*x54 + x62 + x92)
            T(i, 42) = x124*(x1*x10*x91 + x114 + x146 + x39)
            T(i, 43) = x15*(-x116*x127 + x128*x130 + x131 + x145 + x22)
            T(i, 44) = x144*x71
            T(i, 45) = x15*(6930.0d0*x109 + x13*zi**8 - 12012.0d0*x143 - 1260.0d0*x77 + 35.0d0)
        end do
    end subroutine T8
    subroutine T9(x, y, z, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), 220)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x6
        real(8) :: x7
        real(8) :: x8
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x13
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
        real(8) :: x27
        real(8) :: x28
        real(8) :: x29
        real(8) :: x30
        real(8) :: x31
        real(8) :: x32
        real(8) :: x33
        real(8) :: x34
        real(8) :: x35
        real(8) :: x36
        real(8) :: x37
        real(8) :: x38
        real(8) :: x39
        real(8) :: x40
        real(8) :: x41
        real(8) :: x42
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
        real(8) :: x53
        real(8) :: x54
        real(8) :: x55
        real(8) :: x56
        real(8) :: x57
        real(8) :: x58
        real(8) :: x59
        real(8) :: x60
        real(8) :: x61
        real(8) :: x62
        real(8) :: x63
        real(8) :: x64
        real(8) :: x65
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
        real(8) :: x80
        real(8) :: x81
        real(8) :: x82
        real(8) :: x83
        real(8) :: x84
        real(8) :: x85
        real(8) :: x86
        real(8) :: x87
        real(8) :: x88
        real(8) :: x89
        real(8) :: x90
        real(8) :: x91
        real(8) :: x92
        real(8) :: x93
        real(8) :: x94
        real(8) :: x95
        real(8) :: x96
        real(8) :: x97
        real(8) :: x98
        real(8) :: x99
        real(8) :: x100
        real(8) :: x101
        real(8) :: x102
        real(8) :: x103
        real(8) :: x104
        real(8) :: x105
        real(8) :: x106
        real(8) :: x107
        real(8) :: x108
        real(8) :: x109
        real(8) :: x110
        real(8) :: x111
        real(8) :: x112
        real(8) :: x113
        real(8) :: x114
        real(8) :: x115
        real(8) :: x116
        real(8) :: x117
        real(8) :: x118
        real(8) :: x119
        real(8) :: x120
        real(8) :: x121
        real(8) :: x122
        real(8) :: x123
        real(8) :: x124
        real(8) :: x125
        real(8) :: x126
        real(8) :: x127
        real(8) :: x128
        real(8) :: x129
        real(8) :: x130
        real(8) :: x131
        real(8) :: x132
        real(8) :: x133
        real(8) :: x134
        real(8) :: x135
        real(8) :: x136
        real(8) :: x137
        real(8) :: x138
        real(8) :: x139
        real(8) :: x140
        real(8) :: x141
        real(8) :: x142
        real(8) :: x143
        real(8) :: x144
        real(8) :: x145
        real(8) :: x146
        real(8) :: x147
        real(8) :: x148
        real(8) :: x149
        real(8) :: x150
        real(8) :: x151
        real(8) :: x152
        real(8) :: x153
        real(8) :: x154
        real(8) :: x155
        real(8) :: x156
        real(8) :: x157
        real(8) :: x158
        real(8) :: x159
        real(8) :: x160
        real(8) :: x161
        real(8) :: x162
        real(8) :: x163
        real(8) :: x164
        real(8) :: x165
        real(8) :: x166
        real(8) :: x169
        real(8) :: x170
        real(8) :: x173
        real(8) :: x175
        real(8) :: x176
        real(8) :: x177
        real(8) :: x178
        real(8) :: x179
        real(8) :: x180
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
        do concurrent(i=1:size(x))
            xi = x(i)
            yi = y(i)
            zi = z(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0/x3
            x5 = x0*x4
            x6 = xi**4
            x7 = x3**(-2)
            x8 = x6*x7
            x9 = xi**6
            x10 = x3**(-3)
            x11 = x10*x9
            x12 = x3**(-4)
            x13 = x12*xi**8
            x14 = x3**(-5.5d0)
            x15 = x14*xi
            x16 = 2835.0d0*x15
            x17 = 14175.0d0*x14
            x18 = x17*(-4004.0d0*x11 + 2431.0d0*x13 - 308.0d0*x5 + 2002.0d0*x8 + 7.0d0)
            x19 = x10*x6
            x20 = 15015.0d0*x19
            x21 = 12155.0d0*x12
            x22 = x21*x9
            x23 = 35.0d0 - 385.0d0*x5
            x24 = -715.0d0*x11 + x23 + 1001.0d0*x8
            x25 = 385.0d0*x4
            x26 = x1*x25
            x27 = -x26
            x28 = x0*x7
            x29 = 5005.0d0*x28
            x30 = x1*x29 + x27
            x31 = x3**(-6.5d0)*xi*yi*zi
            x32 = 155925.0d0*x31
            x33 = -x2*x25
            x34 = x2*x29 + x33
            x35 = 429.0d0*x28
            x36 = x1*x35
            x37 = x36 + 3.0d0
            x38 = 11.0d0*x4
            x39 = x1*x38
            x40 = -x39
            x41 = 2145.0d0*x19
            x42 = -x1*x41
            x43 = 2431.0d0*x12
            x44 = x43*x9
            x45 = x1*x44 + x40 + x42
            x46 = -429.0d0*x11 - 99.0d0*x5 + 429.0d0*x8
            x47 = x17*yi
            x48 = x36 + 1.0d0
            x49 = -33.0d0*x5
            x50 = 143.0d0*x8
            x51 = -143.0d0*x11 + x49 + x50
            x52 = x17*zi
            x53 = x2*x35
            x54 = x53 + 1.0d0
            x55 = x2*x38
            x56 = -x55
            x57 = -x2*x41
            x58 = x2*x44 + x56 + x57
            x59 = x53 + 3.0d0
            x60 = yi**4
            x61 = x0*x10
            x62 = 7150.0d0*x61
            x63 = x1*x19
            x64 = -4290.0d0*x63
            x65 = x1*x28
            x66 = x21*x6
            x67 = x60*x66
            x68 = 2860.0d0*x65 + x67
            x69 = x1*x4
            x70 = 715.0d0*x7
            x71 = x60*x70 - 330.0d0*x69
            x72 = -110.0d0*x5 + x50 + 15.0d0
            x73 = 26.0d0*x5 - 39.0d0*x8
            x74 = -130.0d0*x65 - 3.0d0
            x75 = x1*x2
            x76 = x7*x75
            x77 = 36465.0d0*x12
            x78 = x2*x4
            x79 = x61*x75
            x80 = -21450.0d0*x79
            x81 = 1430.0d0*x65 - 165.0d0*x78 + x80
            x82 = x2*x28
            x83 = -165.0d0*x69 + 1430.0d0*x82
            x84 = 945.0d0*x15
            x85 = -130.0d0*x82
            x86 = x19*x2
            x87 = 13.0d0*x78 - 3.0d0
            x88 = zi**4
            x89 = 2860.0d0*x82
            x90 = x66*x88
            x91 = -4290.0d0*x86 + x90
            x92 = x70*x88 - 330.0d0*x78
            x93 = 4290.0d0*x61
            x94 = -x60*x93
            x95 = 110.0d0*x69
            x96 = 143.0d0*x7
            x97 = x60*x96
            x98 = -x95 + x97 + 15.0d0
            x99 = -330.0d0*x5 + 715.0d0*x8
            x100 = 2835.0d0*x14
            x101 = x100*yi
            x102 = -66.0d0*x5 + x50
            x103 = x102 + 3.0d0
            x104 = -66.0d0*x69 + x97
            x105 = x100*zi
            x106 = 33.0d0*x78
            x107 = -x106
            x108 = -x75*x93
            x109 = x107 + x108 + 286.0d0*x65
            x110 = x40 + 858.0d0*x82
            x111 = x103 + x66*x75 + x75*x96
            x112 = -33.0d0*x69
            x113 = x108 + x112 + 286.0d0*x82
            x114 = x56 + 858.0d0*x65
            x115 = -x88*x93
            x116 = x88*x96
            x117 = x116 - 66.0d0*x78 + 3.0d0
            x118 = x116 - 110.0d0*x78 + 15.0d0
            x119 = -11.0d0*x5
            x120 = x60*x61
            x121 = -2145.0d0*x120
            x122 = yi**6
            x123 = x0*x43
            x124 = x119 + x121 + x122*x123
            x125 = 429.0d0*x10
            x126 = x122*x125
            x127 = x60*x7
            x128 = -x126 + 429.0d0*x127 - 99.0d0*x69
            x129 = 14175.0d0*x15
            x130 = 13.0d0*x5
            x131 = 2145.0d0*x10
            x132 = x2*x60
            x133 = -x131*x132
            x134 = x119 + 858.0d0*x76
            x135 = x0*x96
            x136 = x0*x21
            x137 = x104 + x132*x136 + x135*x2 + 3.0d0
            x138 = 31185.0d0*x31
            x139 = x1*x88
            x140 = -x131*x139
            x141 = x61*x88
            x142 = x1*x135 + x117 + x136*x139
            x143 = 26.0d0*x78
            x144 = x7*x88
            x145 = 39.0d0*x144
            x146 = -2145.0d0*x141
            x147 = zi**6
            x148 = x119 + x123*x147 + x146
            x149 = -x125*x147 + 429.0d0*x144 - 99.0d0*x78
            x150 = x10*x122
            x151 = 715.0d0*x150
            x152 = 1001.0d0*x127
            x153 = -x151 + x152
            x154 = x0*x77
            x155 = -165.0d0*x5 + 1430.0d0*x76
            x156 = 945.0d0*x14
            x157 = 715.0d0*x10
            x158 = x108 + x49 + 286.0d0*x76
            x159 = x10*x147
            x160 = 143.0d0*x159
            x161 = 1001.0d0*x144 - 715.0d0*x159
            x162 = yi**8
            x163 = 2002.0d0*x127 - 4004.0d0*x150 + x162*x43 - 308.0d0*x69 + 7.0d0
            x164 = 2.0d0*x69
            x165 = x55 - 1.0d0
            x166 = x10*x132
            x169 = 8.0d0*x69
            x170 = x55 - 3.0d0
            x173 = 4.0d0*x78
            x175 = 9.0d0*x78 - 1.0d0
            x176 = 56.0d0*x69
            x177 = 24.0d0*x78
            x178 = x10*x139
            x179 = 16.0d0*x144
            x180 = x179 - 16.0d0*x78 + 3.0d0
            x181 = 16.0d0*x69
            x182 = x179 - 12.0d0*x78 + 1.0d0
            x183 = 2.0d0*x78 - 1.0d0
            x184 = x69*(x173 - 1.0d0)
            x185 = zi**8
            x186 = 2002.0d0*x144 - 4004.0d0*x159 + x185*x43 - 308.0d0*x78 + 7.0d0
            x187 = x122*x2
            x188 = x27 + x33 + 5005.0d0*x76 + 35.0d0
            x189 = 429.0d0*x76 + 3.0d0
            x190 = x21*x60*x88 + 2860.0d0*x76
            x191 = x1*x147
            T(i, 1) = -x16*(-25740.0d0*x11 + 12155.0d0*x13 - 4620.0d0*x5 + 18018.0d0*x8 + 315.0d0)
            T(i, 2) = -x18*yi
            T(i, 3) = -x18*zi
            T(i, 4) = -x16*(-x1*x20 + x1*x22 + x24 + x30)
            T(i, 5) = -x32*(221.0d0*x11 + 91.0d0*x5 - 273.0d0*x8 - 7.0d0)
            T(i, 6) = -x16*(-x2*x20 + x2*x22 + x24 + x34)
            T(i, 7) = -x47*(x37 + x45 + x46)
            T(i, 8) = -x52*(x45 + x48 + x51)
            T(i, 9) = -x47*(x51 + x54 + x58)
            T(i, 10) = -x52*(x46 + x58 + x59)
            T(i, 11) = -x16*(-x60*x62 + x64 + x68 + x71 + x72)
            T(i, 12) = -x32*(221.0d0*x63 + 13.0d0*x69 + x73 + x74)
            T(i, 13) = -x84*(x42 + x57 + x6*x75*x77 + x72 + 2145.0d0*x76 + x81 + x83)
            T(i, 14) = -x32*(x73 + x85 + 221.0d0*x86 + x87)
            T(i, 15) = -x16*(-x62*x88 + x72 + x89 + x91 + x92)
            T(i, 16) = -x101*(-7150.0d0*x63 + x68 + x94 + x98 + x99)
            T(i, 17) = -x105*(x103 + x104 + x64 + 1716.0d0*x65 + x67 + x94)
            T(i, 18) = -x101*(x109 + x110 + x111 + x57 - 715.0d0*x63)
            T(i, 19) = -x105*(x111 + x113 + x114 + x42 - 715.0d0*x86)
            T(i, 20) = -x101*(x102 + x115 + x117 + 1716.0d0*x82 + x91)
            T(i, 21) = -x105*(x115 + x118 - 7150.0d0*x86 + x89 + x90 + x99)
            T(i, 22) = -x129*(x124 + x128 + x37)
            T(i, 23) = -x32*(221.0d0*x120 - 39.0d0*x127 + x130 + 26.0d0*x69 + x74)
            T(i, 24) = -x16*(x109 - 715.0d0*x120 + x133 + x134 + x137)
            T(i, 25) = -x138*(39.0d0*x5 - 195.0d0*x65 + 39.0d0*x69 - 195.0d0*x76 + 39.0d0*x78 + 1105.0d0*x79 - 195.0d0*x82 - 9.0d0)
            T(i, 26) = -x16*(x113 + x134 + x140 - 715.0d0*x141 + x142)
            T(i, 27) = -x32*(x130 + 221.0d0*x141 + x143 - x145 + x85 - 3.0d0)
            T(i, 28) = -x129*(x148 + x149 + x59)
            T(i, 29) = -x101*(-15015.0d0*x120 + x122*x136 + x153 + x23 + x30)
            T(i, 30) = -x52*(x112 + x124 - 143.0d0*x150 + x48 + x97)
            T(i, 31) = -x156*yi*(x121 + x132*x154 + x133 + x155 + x81 + 2145.0d0*x82 + x98)
            T(i, 32) = -x105*(x114 + x121 - x132*x157 + x137 + x158)
            T(i, 33) = -x101*(x110 - x139*x157 + x142 + x146 + x158)
            T(i, 34) = -x156*zi*(x118 + x139*x154 + x140 + x146 + x155 + 2145.0d0*x65 + x80 + x83)
            T(i, 35) = -x47*(x107 + x116 + x148 - x160 + x54)
            T(i, 36) = -x105*(x136*x147 - 15015.0d0*x141 + x161 + x23 + x34)
            T(i, 37) = -x129*x163
            T(i, 38) = -x138*(x151 - x152 + x164*(195.0d0*x127 - 182.0d0*x69 + 35.0d0) + x26 - 35.0d0)
            T(i, 39) = -x16*(x164*(1144.0d0*x166 - 660.0d0*x76 + 60.0d0*x78 + x95 - x97 - 15.0d0) + x165*(x126 - 495.0d0*x127 + 135.0d0*x69 - 5.0d0) + 36.0d0*x76*x98)
            T(i, 40) = -2835.0d0*x31*(36.0d0*x165*x69*(x39 - 5.0d0) + x169*(x112 + 88.0d0*x76 - 20.0d0*x78 + 10.0d0) + 36.0d0*x69*(x40 + 66.0d0*x76 - 20.0d0*x78 + 5.0d0) + 11.0d0*x87*(33.0d0*x127 - 30.0d0*x69 + 5.0d0))
            T(i, 41) = -x84*(3.0d0*x117*(21.0d0*x127 - 14.0d0*x69 + 1.0d0) + x169*(-32.0d0*x144 + x177 + 240.0d0*x178 + 9.0d0*x69 - 144.0d0*x76 - 2.0d0) + 336.0d0*x170*x76*(3.0d0*x69 - 1.0d0) + x175*x176*(-x173 - 3.0d0*x69 + 18.0d0*x76 + 1.0d0) + 224.0d0*x76*(-x173 - 9.0d0*x69 + 24.0d0*x76 + 2.0d0))
            T(i, 42) = -1575.0d0*x31*(42.0d0*x117*x69 + 84.0d0*x170*x184 + 112.0d0*x175*x183*x69 + x176*x182 + x180*x181 + 33.0d0*(7.0d0*x69 - 3.0d0)*(-x143 + x145 + 3.0d0))
            T(i, 43) = -675.0d0*x15*(84.0d0*x118*x76 + 160.0d0*x180*x76 + x181*(-80.0d0*x144 + 64.0d0*x159 + x177 - 1.0d0) + 40.0d0*x182*x69*(7.0d0*x78 - 1.0d0) + 1120.0d0*x183*x76*(3.0d0*x78 - 1.0d0) + 70.0d0*x184*(33.0d0*x144 - 18.0d0*x78 + 1.0d0) + 21.0d0*(5.0d0*x69 - 1.0d0)*(x106 - x116 + x160 - 1.0d0))
            T(i, 44) = -x32*(-273.0d0*x144 + 221.0d0*x159 + 91.0d0*x78 - 7.0d0)
            T(i, 45) = -x129*x186
            T(i, 46) = -x101*(18018.0d0*x127 - 25740.0d0*x150 + x162*x21 - 4620.0d0*x69 + 315.0d0)
            T(i, 47) = -x163*x52
            T(i, 48) = -x101*(x153 - 15015.0d0*x166 + x187*x21 + x188)
            T(i, 49) = -x52*(x128 + x133 + x187*x43 + x189 + x56)
            T(i, 50) = -x101*(-4290.0d0*x166 - 7150.0d0*x178 + x190 + x92 + x98)
            T(i, 51) = -x105*(x118 - 7150.0d0*x166 - 4290.0d0*x178 + x190 + x71)
            T(i, 52) = -x47*(x140 + x149 + x189 + x191*x43 + x40)
            T(i, 53) = -x105*(x161 - 15015.0d0*x178 + x188 + x191*x21)
            T(i, 54) = -x186*x47
            T(i, 55) = -x105*(18018.0d0*x144 - 25740.0d0*x159 + x185*x21 - 4620.0d0*x78 + 315.0d0)
        end do
    end subroutine T9
    subroutine T10(x, y, z, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), 286)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x6
        real(8) :: x7
        real(8) :: x8
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x13
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
        real(8) :: x27
        real(8) :: x28
        real(8) :: x29
        real(8) :: x30
        real(8) :: x31
        real(8) :: x32
        real(8) :: x33
        real(8) :: x34
        real(8) :: x35
        real(8) :: x36
        real(8) :: x37
        real(8) :: x38
        real(8) :: x39
        real(8) :: x40
        real(8) :: x41
        real(8) :: x42
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
        real(8) :: x53
        real(8) :: x54
        real(8) :: x55
        real(8) :: x56
        real(8) :: x57
        real(8) :: x58
        real(8) :: x59
        real(8) :: x60
        real(8) :: x61
        real(8) :: x62
        real(8) :: x63
        real(8) :: x64
        real(8) :: x65
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
        real(8) :: x80
        real(8) :: x81
        real(8) :: x82
        real(8) :: x83
        real(8) :: x84
        real(8) :: x85
        real(8) :: x86
        real(8) :: x87
        real(8) :: x88
        real(8) :: x89
        real(8) :: x90
        real(8) :: x91
        real(8) :: x92
        real(8) :: x93
        real(8) :: x94
        real(8) :: x95
        real(8) :: x96
        real(8) :: x97
        real(8) :: x98
        real(8) :: x99
        real(8) :: x100
        real(8) :: x101
        real(8) :: x102
        real(8) :: x103
        real(8) :: x104
        real(8) :: x105
        real(8) :: x106
        real(8) :: x107
        real(8) :: x108
        real(8) :: x109
        real(8) :: x110
        real(8) :: x111
        real(8) :: x112
        real(8) :: x113
        real(8) :: x114
        real(8) :: x115
        real(8) :: x116
        real(8) :: x117
        real(8) :: x118
        real(8) :: x119
        real(8) :: x120
        real(8) :: x121
        real(8) :: x122
        real(8) :: x123
        real(8) :: x124
        real(8) :: x125
        real(8) :: x126
        real(8) :: x127
        real(8) :: x128
        real(8) :: x129
        real(8) :: x130
        real(8) :: x131
        real(8) :: x132
        real(8) :: x133
        real(8) :: x134
        real(8) :: x135
        real(8) :: x136
        real(8) :: x137
        real(8) :: x138
        real(8) :: x139
        real(8) :: x140
        real(8) :: x141
        real(8) :: x142
        real(8) :: x143
        real(8) :: x144
        real(8) :: x145
        real(8) :: x146
        real(8) :: x147
        real(8) :: x148
        real(8) :: x149
        real(8) :: x150
        real(8) :: x151
        real(8) :: x152
        real(8) :: x153
        real(8) :: x154
        real(8) :: x155
        real(8) :: x156
        real(8) :: x157
        real(8) :: x158
        real(8) :: x159
        real(8) :: x160
        real(8) :: x161
        real(8) :: x162
        real(8) :: x163
        real(8) :: x164
        real(8) :: x165
        real(8) :: x166
        real(8) :: x167
        real(8) :: x168
        real(8) :: x169
        real(8) :: x170
        real(8) :: x171
        real(8) :: x172
        real(8) :: x173
        real(8) :: x174
        real(8) :: x175
        real(8) :: x176
        real(8) :: x177
        real(8) :: x178
        real(8) :: x179
        real(8) :: x180
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
        real(8) :: x196
        real(8) :: x197
        real(8) :: x198
        real(8) :: x199
        real(8) :: x200
        real(8) :: x201
        real(8) :: x202
        real(8) :: x203
        real(8) :: x204
        real(8) :: x205
        real(8) :: x206
        real(8) :: x207
        real(8) :: x208
        real(8) :: x209
        real(8) :: x210
        real(8) :: x211
        real(8) :: x212
        real(8) :: x213
        real(8) :: x214
        real(8) :: x215
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
        real(8) :: x227
        real(8) :: x228
        real(8) :: x229
        real(8) :: x230
        real(8) :: x231
        real(8) :: x232
        real(8) :: x233
        real(8) :: x234
        real(8) :: x235
        real(8) :: x236
        real(8) :: x237
        real(8) :: x239
        real(8) :: x240
        real(8) :: x241
        real(8) :: x242
        real(8) :: x243
        real(8) :: x244
        real(8) :: x245
        real(8) :: x246
        real(8) :: x247
        real(8) :: x248
        real(8) :: x252
        real(8) :: x255
        real(8) :: x257
        real(8) :: x258
        real(8) :: x260
        real(8) :: x262
        real(8) :: x263
        real(8) :: x264
        real(8) :: x265
        real(8) :: x266
        real(8) :: x267
        real(8) :: x268
        real(8) :: x269
        real(8) :: x270
        real(8) :: x271
        real(8) :: x272
        real(8) :: x274
        real(8) :: x275
        real(8) :: x276
        real(8) :: x277
        real(8) :: x278
        real(8) :: x279
        do concurrent(i=1:size(x))
            xi = x(i)
            yi = y(i)
            zi = z(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0/x3
            x5 = x0*x4
            x6 = xi**4
            x7 = x3**(-2)
            x8 = x6*x7
            x9 = xi**6
            x10 = x3**(-3)
            x11 = x10*x9
            x12 = xi**8
            x13 = x3**(-4)
            x14 = x12*x13
            x15 = x3**(-5)
            x16 = 46189.0d0*x15
            x17 = x3**(-5.5d0)
            x18 = 14175.0d0*x17
            x19 = x3**(-6.5d0)
            x20 = x19*yi
            x21 = 155925.0d0*x20
            x22 = 4199.0d0*x14
            x23 = xi*(-7956.0d0*x11 + x22 - 1092.0d0*x5 + 4914.0d0*x8 + 63.0d0)
            x24 = x19*zi
            x25 = 155925.0d0*x24
            x26 = 77.0d0*x4
            x27 = x1*x26
            x28 = x1*x9
            x29 = 68068.0d0*x13
            x30 = x10*x6
            x31 = 30030.0d0*x30
            x32 = x12*x16
            x33 = x0*x7
            x34 = 4004.0d0*x33
            x35 = -x1*x34 - 7.0d0
            x36 = 4004.0d0*x11 - 2431.0d0*x14 + 308.0d0*x5 - 2002.0d0*x8
            x37 = x21*zi
            x38 = x2*x26
            x39 = x2*x9
            x40 = -x2*x34 - 7.0d0
            x41 = -663.0d0*x11
            x42 = 4199.0d0*x13
            x43 = x28*x42
            x44 = x41 + x43
            x45 = 1365.0d0*x33
            x46 = x1*x45
            x47 = x46 + 21.0d0
            x48 = 91.0d0*x4
            x49 = x1*x48
            x50 = -x49
            x51 = 4641.0d0*x30
            x52 = -x1*x51 + x50
            x53 = -273.0d0*x5 + 819.0d0*x8
            x54 = x21*xi
            x55 = -91.0d0*x5
            x56 = x55 + 7.0d0
            x57 = x46 + x56
            x58 = -221.0d0*x11 + 273.0d0*x8
            x59 = x25*xi
            x60 = x2*x48
            x61 = -x60
            x62 = x2*x45
            x63 = x56 + x61 + x62
            x64 = x39*x42
            x65 = -x2*x51 + x64
            x66 = x62 + 21.0d0
            x67 = 14586.0d0*x13
            x68 = yi**4
            x69 = x0*x10
            x70 = 6435.0d0*x69
            x71 = x1*x30
            x72 = x16*x9
            x73 = x6*x68
            x74 = 36465.0d0*x13
            x75 = x1*x33
            x76 = -x73*x74 - 2574.0d0*x75 - 3.0d0
            x77 = 143.0d0*x7
            x78 = x68*x77
            x79 = -x78
            x80 = x1*x4
            x81 = x79 + 66.0d0*x80
            x82 = 429.0d0*x11 + 99.0d0*x5 - 429.0d0*x8
            x83 = -117.0d0*x5
            x84 = -13.0d0*x80
            x85 = 585.0d0*x75 + 3.0d0
            x86 = 585.0d0*x8
            x87 = -3315.0d0*x71 + x86
            x88 = -143.0d0*x8
            x89 = 33.0d0*x5
            x90 = 2431.0d0*x13
            x91 = x1*x2
            x92 = 2145.0d0*x71
            x93 = x2*x30
            x94 = 2145.0d0*x93
            x95 = x6*x91
            x96 = x16*x2
            x97 = x70*x91
            x98 = x2*x4
            x99 = 11.0d0*x98
            x100 = x99 - 1.0d0
            x101 = x100 - 429.0d0*x75 + x97
            x102 = 11.0d0*x80
            x103 = x2*x33
            x104 = x102 - 429.0d0*x103
            x105 = 13.0d0*x98
            x106 = -x105
            x107 = -3315.0d0*x93
            x108 = 585.0d0*x103 + 3.0d0
            x109 = zi**4
            x110 = x109*x6
            x111 = -2574.0d0*x103 - x110*x74 - 3.0d0
            x112 = x109*x77
            x113 = -x112
            x114 = 66.0d0*x98
            x115 = x113 + x114
            x116 = -130.0d0*x80
            x117 = -2210.0d0*x71
            x118 = 195.0d0*x7
            x119 = x118*x68
            x120 = 2210.0d0*x69
            x121 = x42*x73
            x122 = x119 - x120*x68 + x121
            x123 = 195.0d0*x8
            x124 = x123 - 130.0d0*x5 + 15.0d0
            x125 = 780.0d0*x75 + 3.0d0
            x126 = -26.0d0*x5 + 39.0d0*x8
            x127 = 663.0d0*x10
            x128 = x127*x6
            x129 = 39.0d0*x98
            x130 = -x129
            x131 = x118*x91
            x132 = x130 + x131
            x133 = -x120*x91
            x134 = x133 + 130.0d0*x75
            x135 = 390.0d0*x103 + x84
            x136 = x126 + 3.0d0
            x137 = x136 + x42*x95
            x138 = -39.0d0*x80
            x139 = 130.0d0*x103 + x133 + x138
            x140 = x106 + 390.0d0*x75
            x141 = x110*x42
            x142 = 780.0d0*x103 + x141
            x143 = x109*x118
            x144 = -x109*x120 + x143
            x145 = -130.0d0*x98
            x146 = -2210.0d0*x93
            x147 = yi**6
            x148 = x0*x67
            x149 = x68*x69
            x150 = x16*x6
            x151 = 66.0d0*x5 + x88
            x152 = 429.0d0*x10
            x153 = x147*x152
            x154 = x68*x7
            x155 = x153 - 429.0d0*x154 + 99.0d0*x80
            x156 = x123 - 78.0d0*x5
            x157 = 39.0d0*x154 - 26.0d0*x80
            x158 = 33.0d0*x98
            x159 = 2145.0d0*x10
            x160 = x2*x68
            x161 = x159*x160
            x162 = 12155.0d0*x13
            x163 = 230945.0d0*x15
            x164 = x81 - 3.0d0
            x165 = 72930.0d0*x13
            x166 = x0*x160
            x167 = x69*x91
            x168 = 25740.0d0*x167
            x169 = -858.0d0*x103 - x165*x166 + x168
            x170 = x7*x91
            x171 = x151 - x165*x95 - 858.0d0*x170
            x172 = 2835.0d0*x17
            x173 = 20995.0d0*x13
            x174 = -6630.0d0*x167
            x175 = x174 + 1170.0d0*x75 + 9.0d0
            x176 = 1170.0d0*x103 + x138
            x177 = 31185.0d0*x20
            x178 = 33.0d0*x80
            x179 = x1*x109
            x180 = x159*x179
            x181 = x109*x69
            x182 = x0*x179
            x183 = x115 - 3.0d0
            x184 = -x165*x182 + x183 - 858.0d0*x75
            x185 = -26.0d0*x98
            x186 = x109*x7
            x187 = x185 + 39.0d0*x186 + 3.0d0
            x188 = zi**6
            x189 = x152*x188 - 429.0d0*x186 + 99.0d0*x98
            x190 = -4641.0d0*x149
            x191 = -x127*x147
            x192 = x0*x42
            x193 = x147*x192
            x194 = x191 + x193
            x195 = 819.0d0*x154 - 273.0d0*x80
            x196 = -13.0d0*x5
            x197 = -3315.0d0*x149 + 585.0d0*x154
            x198 = x0*x118
            x199 = x198*x2
            x200 = x130 + x199
            x201 = 390.0d0*x170 + x196
            x202 = x157 + x160*x192 + 3.0d0
            x203 = 3315.0d0*x10
            x204 = -39.0d0*x5
            x205 = 1170.0d0*x170 + x204
            x206 = x24*xi
            x207 = x1*x198
            x208 = -3315.0d0*x181 + 585.0d0*x186
            x209 = x177*xi
            x210 = x179*x192 + x187 + x207
            x211 = -x127*x188
            x212 = x188*x192
            x213 = x211 + x212
            x214 = -4641.0d0*x181
            x215 = 819.0d0*x186 - 273.0d0*x98
            x216 = 77.0d0*x5
            x217 = x0*x147
            x218 = yi**8
            x219 = x0*x16
            x220 = x218*x90
            x221 = 2002.0d0*x154
            x222 = 308.0d0*x80
            x223 = x10*x147
            x224 = 4004.0d0*x223
            x225 = -x220 - x221 + x222 + x224
            x226 = 221.0d0*x223
            x227 = 273.0d0*x154
            x228 = x147*x2
            x229 = x0*x77
            x230 = 2145.0d0*x149
            x231 = -429.0d0*x170 + 11.0d0*x5
            x232 = x0*x127
            x233 = x10*x160
            x234 = x133 + 130.0d0*x170 + x204
            x235 = 2145.0d0*x181
            x236 = x109*x68
            x237 = x10*x179
            x239 = x1*x188
            x240 = x10*x188
            x241 = x113 + x158 + 143.0d0*x240 - 1.0d0
            x242 = 221.0d0*x240
            x243 = 273.0d0*x186
            x244 = x188*x29
            x245 = zi**8
            x246 = -2002.0d0*x186 + 4004.0d0*x240 - x245*x90 + 308.0d0*x98
            x247 = 4914.0d0*x154 + x218*x42 - 7956.0d0*x223 - 1092.0d0*x80 + 63.0d0
            x248 = 8.0d0*x80
            x252 = x105 - 3.0d0
            x255 = x20*xi
            x257 = x80*(x112 - x114 + 3.0d0)
            x258 = 9.0d0*x98 - 1.0d0
            x260 = 4.0d0*x98
            x262 = x99 - 3.0d0
            x263 = 16.0d0*x98
            x264 = 16.0d0*x186
            x265 = 24.0d0*x98
            x266 = 56.0d0*x80
            x267 = 4725.0d0*x206
            x268 = -x263 + x264 + 3.0d0
            x269 = 2.0d0*x98 - 1.0d0
            x270 = x80*(-80.0d0*x186 + 64.0d0*x240 + x265 - 1.0d0)
            x271 = 42.0d0*x260 - 42.0d0
            x272 = x264 - 12.0d0*x98 + 1.0d0
            x274 = x245*x42
            x275 = 4914.0d0*x186 - 7956.0d0*x240 + x274 - 1092.0d0*x98 + 63.0d0
            x276 = 109395.0d0*x13
            x277 = -4004.0d0*x170 - 7.0d0
            x278 = 1365.0d0*x170 + 21.0d0
            x279 = -2574.0d0*x170 - x236*x74
            T(i, 1) = x18*(90090.0d0*x11 - 109395.0d0*x14 + x16*xi**10 + 3465.0d0*x5 - 30030.0d0*x8 - 63.0d0)
            T(i, 2) = x21*x23
            T(i, 3) = x23*x25
            T(i, 4) = x18*(x1*x31 + x1*x32 + x27 - x28*x29 + x35 + x36)
            T(i, 5) = x37*(-6188.0d0*x11 + x22 - 364.0d0*x5 + 2730.0d0*x8 + 7.0d0)
            T(i, 6) = x18*(x2*x31 + x2*x32 - x29*x39 + x36 + x38 + x40)
            T(i, 7) = x54*(x44 + x47 + x52 + x53)
            T(i, 8) = x59*(x43 + x52 + x57 + x58)
            T(i, 9) = x54*(x58 + x63 + x65)
            T(i, 10) = x59*(x41 + x53 + x61 + x65 + x66)
            T(i, 11) = x18*(-x28*x67 + x68*x70 + x68*x72 + 12870.0d0*x71 + x76 + x81 + x82)
            T(i, 12) = x37*(x44 + x83 + x84 + x85 + x87)
            T(i, 13) = x18*(x101 + x104 + 143.0d0*x11 - x28*x90 + x28*x96 - x39*x90 - x74*x95 - x77*x91 + x88 + x89 + x92 + x94)
            T(i, 14) = x37*(x106 + x107 + x108 + x41 + x64 + x83 + x86)
            T(i, 15) = x18*(x109*x70 + x109*x72 + x111 + x115 - x39*x67 + x82 + 12870.0d0*x93)
            T(i, 16) = x54*(x116 + x117 + x122 + x124 + 1300.0d0*x75)
            T(i, 17) = x59*(x122 + x125 + x126 - 1326.0d0*x71 - 78.0d0*x80)
            T(i, 18) = x54*(-x128*x2 + x132 + x134 + x135 + x137 - 221.0d0*x71)
            T(i, 19) = x59*(-x1*x128 + x131 + x137 + x139 + x140 - 221.0d0*x93)
            T(i, 20) = x54*(x136 + x142 + x144 - 1326.0d0*x93 - 78.0d0*x98)
            T(i, 21) = x59*(1300.0d0*x103 + x124 + x141 + x144 + x145 + x146)
            T(i, 22) = x18*(-x147*x148 + x147*x150 + 12870.0d0*x149 + x151 + x155 + 6435.0d0*x71 + x76)
            T(i, 23) = x37*(x117 + x121 + x125 - 1326.0d0*x149 + x156 + x157)
            T(i, 24) = x172*(4290.0d0*x149 + x158 + x161 - x162*x73 + x163*x2*x73 + x164 + x169 + x171 + 4290.0d0*x71 - 1716.0d0*x75 + x94)
            T(i, 25) = x177*zi*(x107 + x132 + x173*x95 + x175 + x176 - 234.0d0*x5 + x87)
            T(i, 26) = x172*(x1*x110*x163 - 1716.0d0*x103 - x110*x162 + x168 + x171 + x178 + x180 + 4290.0d0*x181 + x184 + x92 + 4290.0d0*x93)
            T(i, 27) = x37*(x142 + x146 + x156 - 1326.0d0*x181 + x187)
            T(i, 28) = x18*(x111 - x148*x188 + x150*x188 + x151 + 12870.0d0*x181 + x189 + 6435.0d0*x93)
            T(i, 29) = x54*(x190 + x194 + x195 + x47 + x55)
            T(i, 30) = x59*(x194 + x196 + x197 - 117.0d0*x80 + x85)
            T(i, 31) = x54*(-x127*x160 + x134 - 221.0d0*x149 + x200 + x201 + x202)
            T(i, 32) = 31185.0d0*x206*(-x160*x203 + x166*x173 + x175 + x197 + x200 + x205 - 234.0d0*x80)
            T(i, 33) = x209*(x173*x182 + x174 + x176 - x179*x203 + x205 + x207 + x208 - 234.0d0*x98 + 9.0d0)
            T(i, 34) = x59*(-x127*x179 + x139 - 221.0d0*x181 + x201 + x210)
            T(i, 35) = x54*(x108 + x196 + x208 + x213 - 117.0d0*x98)
            T(i, 36) = x59*(x213 + x214 + x215 + x55 + x66)
            T(i, 37) = x18*(30030.0d0*x149 + x216 - x217*x29 + x218*x219 + x225 + x35)
            T(i, 38) = x37*(x190 + x193 - x226 + x227 + x50 + x57)
            T(i, 39) = x18*(x101 + x161 - x166*x74 + x178 - x2*x229 - x217*x90 + x217*x96 + 143.0d0*x223 - x228*x90 + x230 + x231 + x79)
            T(i, 40) = x37*(x140 + x199 + x202 - x232*x68 - 221.0d0*x233 + x234)
            T(i, 41) = x172*(x0*x163*x236 - x162*x236 + x169 - 1716.0d0*x170 + x184 + x230 + 4290.0d0*x233 + x235 + 4290.0d0*x237 + x81 + x89)
            T(i, 42) = x37*(-x109*x232 + x135 + x210 + x234 - 221.0d0*x237)
            T(i, 43) = x18*(-x0*x188*x90 - x1*x188*x90 - x1*x229 + x104 + x180 - x182*x74 + x219*x239 + x231 + x235 + x241 + x97)
            T(i, 44) = x37*(x212 + x214 - x242 + x243 + x63)
            T(i, 45) = x18*(-x0*x244 + 30030.0d0*x181 + x216 + x219*x245 + x246 + x40)
            T(i, 46) = x247*x54
            T(i, 47) = x59*(x220 + x221 - x222 - x224 + x248*(x226 - x227 + x49 - 7.0d0) + 7.0d0)
            T(i, 48) = x209*(44.0d0*x170*(x119 - 182.0d0*x80 + 35.0d0) + 2.0d0*x80*(-x119 - 1092.0d0*x170 + 1560.0d0*x233 + 182.0d0*x80 + 140.0d0*x98 - 35.0d0) + (x105 - 1.0d0)*(-1001.0d0*x154 + 715.0d0*x223 + 385.0d0*x80 - 35.0d0))
            T(i, 49) = 2835.0d0*x206*(54.0d0*x100*x80*(x78 - 110.0d0*x80 + 15.0d0) + x248*(-858.0d0*x154 - 1320.0d0*x170 + 2860.0d0*x233 + 495.0d0*x80 + 90.0d0*x98 - 45.0d0) + 11.0d0*x252*(x153 - 495.0d0*x154 + 135.0d0*x80 - 5.0d0) + 54.0d0*x80*(-660.0d0*x170 + 1144.0d0*x233 + x79 + 110.0d0*x80 + 60.0d0*x98 - 15.0d0))
            T(i, 50) = 2835.0d0*x255*(72.0d0*x100*x80*(-x102 + 66.0d0*x170 - 20.0d0*x98 + 5.0d0) + 528.0d0*x170*x252*(x102 - 5.0d0) + 288.0d0*x170*(88.0d0*x170 - x178 - 20.0d0*x98 + 10.0d0) + x248*(-528.0d0*x170 + x178 - 160.0d0*x186 + 880.0d0*x237 + 120.0d0*x98 - 10.0d0) + 33.0d0*(33.0d0*x154 - 30.0d0*x80 + 5.0d0)*(x185 + 65.0d0*x186 + 1.0d0))
            T(i, 51) = x267*(33.0d0*x187*(21.0d0*x154 - 14.0d0*x80 + 1.0d0) + 84.0d0*x257*(3.0d0*x80 - 1.0d0) + 112.0d0*x258*x80*(24.0d0*x170 - x260 - 9.0d0*x80 + 2.0d0) + 168.0d0*x262*x80*(18.0d0*x170 - x260 - 3.0d0*x80 + 1.0d0) + x266*(-144.0d0*x170 - 32.0d0*x186 + 240.0d0*x237 + x265 + 9.0d0*x80 - 2.0d0) + 32.0d0*x80*(-120.0d0*x170 + 144.0d0*x237 + x263 - x264 + 18.0d0*x80 - 3.0d0))
            T(i, 52) = 4725.0d0*x255*(924.0d0*x170*x187 + 672.0d0*x170*x262*x269 + 224.0d0*x170*x268 + x257*x271 + x258*x266*x272 + 16.0d0*x270 + 33.0d0*(7.0d0*x80 - 3.0d0)*(x129 - x143 + x242 - 1.0d0))
            T(i, 53) = x267*(210.0d0*x241*x80 + 80.0d0*x268*x80*(7.0d0*x98 - 1.0d0) + 280.0d0*x269*x80*(33.0d0*x186 - 18.0d0*x98 + 1.0d0) + 80.0d0*x270 + x271*x80*(x112 - 110.0d0*x98 + 15.0d0) + 280.0d0*x272*x80*(3.0d0*x98 - 1.0d0) + 128.0d0*x80*(-24.0d0*x186 + 16.0d0*x240 + 10.0d0*x98 - 1.0d0) + 33.0d0*(5.0d0*x80 - 1.0d0)*(x242 - x243 + x60 - 7.0d0))
            T(i, 54) = x54*(2730.0d0*x186 - 6188.0d0*x240 + x274 - 364.0d0*x98 + 7.0d0)
            T(i, 55) = x275*x59
            T(i, 56) = x18*(-30030.0d0*x154 + x16*yi**10 - x218*x276 + 90090.0d0*x223 + 3465.0d0*x80 - 63.0d0)
            T(i, 57) = x247*x37
            T(i, 58) = x18*(x218*x96 + x225 - x228*x29 + 30030.0d0*x233 + x277 + x38)
            T(i, 59) = x37*(x191 + x195 + x228*x42 - 4641.0d0*x233 + x278 + x61)
            T(i, 60) = x18*(x109*x147*x16 + x155 + x183 - x228*x67 + 12870.0d0*x233 + 6435.0d0*x237 + x279)
            T(i, 61) = x37*(x116 + x119 + x143 + x145 + 1300.0d0*x170 - 2210.0d0*x233 + x236*x42 - 2210.0d0*x237 + 15.0d0)
            T(i, 62) = x18*(x16*x188*x68 + x164 + x189 + 6435.0d0*x233 + 12870.0d0*x237 - x239*x67 + x279)
            T(i, 63) = x37*(x211 + x215 - 4641.0d0*x237 + x239*x42 + x278 + x50)
            T(i, 64) = x18*(x1*x16*x245 - x1*x244 + 30030.0d0*x237 + x246 + x27 + x277)
            T(i, 65) = x275*x37
            T(i, 66) = x18*(x16*zi**10 - 30030.0d0*x186 + 90090.0d0*x240 - x245*x276 + 3465.0d0*x98 - 63.0d0)
        end do
    end subroutine T10
    subroutine T11(x, y, z, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), 364)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x6
        real(8) :: x7
        real(8) :: x8
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x13
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
        real(8) :: x27
        real(8) :: x28
        real(8) :: x29
        real(8) :: x30
        real(8) :: x31
        real(8) :: x32
        real(8) :: x33
        real(8) :: x34
        real(8) :: x35
        real(8) :: x36
        real(8) :: x37
        real(8) :: x38
        real(8) :: x39
        real(8) :: x40
        real(8) :: x41
        real(8) :: x42
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
        real(8) :: x53
        real(8) :: x54
        real(8) :: x55
        real(8) :: x56
        real(8) :: x57
        real(8) :: x58
        real(8) :: x59
        real(8) :: x60
        real(8) :: x61
        real(8) :: x62
        real(8) :: x63
        real(8) :: x64
        real(8) :: x65
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
        real(8) :: x80
        real(8) :: x81
        real(8) :: x82
        real(8) :: x83
        real(8) :: x84
        real(8) :: x85
        real(8) :: x86
        real(8) :: x87
        real(8) :: x88
        real(8) :: x89
        real(8) :: x90
        real(8) :: x91
        real(8) :: x92
        real(8) :: x93
        real(8) :: x95
        real(8) :: x96
        real(8) :: x97
        real(8) :: x98
        real(8) :: x99
        real(8) :: x100
        real(8) :: x101
        real(8) :: x102
        real(8) :: x103
        real(8) :: x104
        real(8) :: x105
        real(8) :: x106
        real(8) :: x107
        real(8) :: x108
        real(8) :: x109
        real(8) :: x110
        real(8) :: x111
        real(8) :: x112
        real(8) :: x113
        real(8) :: x114
        real(8) :: x115
        real(8) :: x116
        real(8) :: x117
        real(8) :: x118
        real(8) :: x119
        real(8) :: x120
        real(8) :: x121
        real(8) :: x122
        real(8) :: x123
        real(8) :: x124
        real(8) :: x125
        real(8) :: x126
        real(8) :: x127
        real(8) :: x128
        real(8) :: x129
        real(8) :: x130
        real(8) :: x131
        real(8) :: x132
        real(8) :: x133
        real(8) :: x134
        real(8) :: x135
        real(8) :: x136
        real(8) :: x137
        real(8) :: x138
        real(8) :: x139
        real(8) :: x140
        real(8) :: x141
        real(8) :: x142
        real(8) :: x143
        real(8) :: x144
        real(8) :: x145
        real(8) :: x146
        real(8) :: x147
        real(8) :: x148
        real(8) :: x149
        real(8) :: x150
        real(8) :: x151
        real(8) :: x152
        real(8) :: x153
        real(8) :: x154
        real(8) :: x155
        real(8) :: x156
        real(8) :: x157
        real(8) :: x158
        real(8) :: x159
        real(8) :: x160
        real(8) :: x161
        real(8) :: x162
        real(8) :: x163
        real(8) :: x164
        real(8) :: x165
        real(8) :: x166
        real(8) :: x167
        real(8) :: x168
        real(8) :: x169
        real(8) :: x170
        real(8) :: x171
        real(8) :: x172
        real(8) :: x173
        real(8) :: x174
        real(8) :: x175
        real(8) :: x176
        real(8) :: x177
        real(8) :: x178
        real(8) :: x179
        real(8) :: x180
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
        real(8) :: x196
        real(8) :: x197
        real(8) :: x198
        real(8) :: x199
        real(8) :: x200
        real(8) :: x201
        real(8) :: x202
        real(8) :: x203
        real(8) :: x204
        real(8) :: x205
        real(8) :: x206
        real(8) :: x207
        real(8) :: x208
        real(8) :: x209
        real(8) :: x210
        real(8) :: x211
        real(8) :: x212
        real(8) :: x213
        real(8) :: x214
        real(8) :: x215
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
        real(8) :: x227
        real(8) :: x228
        real(8) :: x229
        real(8) :: x230
        real(8) :: x231
        real(8) :: x232
        real(8) :: x233
        real(8) :: x234
        real(8) :: x235
        real(8) :: x236
        real(8) :: x237
        real(8) :: x238
        real(8) :: x239
        real(8) :: x240
        real(8) :: x241
        real(8) :: x242
        real(8) :: x243
        real(8) :: x244
        real(8) :: x245
        real(8) :: x246
        real(8) :: x247
        real(8) :: x248
        real(8) :: x249
        real(8) :: x250
        real(8) :: x251
        real(8) :: x252
        real(8) :: x253
        real(8) :: x254
        real(8) :: x255
        real(8) :: x256
        real(8) :: x257
        real(8) :: x258
        real(8) :: x259
        real(8) :: x260
        real(8) :: x261
        real(8) :: x262
        real(8) :: x263
        real(8) :: x264
        real(8) :: x265
        real(8) :: x266
        real(8) :: x267
        real(8) :: x268
        real(8) :: x269
        real(8) :: x270
        real(8) :: x271
        real(8) :: x272
        real(8) :: x273
        real(8) :: x274
        real(8) :: x275
        real(8) :: x276
        real(8) :: x277
        real(8) :: x278
        real(8) :: x279
        real(8) :: x280
        real(8) :: x281
        real(8) :: x282
        real(8) :: x283
        real(8) :: x284
        real(8) :: x285
        real(8) :: x286
        real(8) :: x287
        real(8) :: x288
        real(8) :: x289
        real(8) :: x290
        real(8) :: x294
        real(8) :: x295
        real(8) :: x296
        real(8) :: x297
        real(8) :: x298
        real(8) :: x299
        real(8) :: x300
        real(8) :: x301
        real(8) :: x303
        real(8) :: x306
        real(8) :: x307
        real(8) :: x308
        real(8) :: x309
        real(8) :: x311
        real(8) :: x313
        real(8) :: x315
        real(8) :: x316
        real(8) :: x317
        real(8) :: x318
        real(8) :: x319
        real(8) :: x320
        real(8) :: x321
        real(8) :: x322
        real(8) :: x323
        real(8) :: x324
        real(8) :: x327
        real(8) :: x328
        real(8) :: x329
        real(8) :: x330
        real(8) :: x331
        real(8) :: x332
        real(8) :: x333
        real(8) :: x334
        real(8) :: x335
        real(8) :: x336
        real(8) :: x337
        real(8) :: x338
        real(8) :: x339
        real(8) :: x340
        real(8) :: x341
        do concurrent(i=1:size(x))
            xi = x(i)
            yi = y(i)
            zi = z(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0/x3
            x5 = x0*x4
            x6 = xi**4
            x7 = x3**(-2)
            x8 = x6*x7
            x9 = xi**6
            x10 = x3**(-3)
            x11 = x10*x9
            x12 = xi**8
            x13 = x3**(-4)
            x14 = x12*x13
            x15 = x3**(-5)
            x16 = x15*xi**10
            x17 = x3**(-6.5d0)
            x18 = 155925.0d0*x17
            x19 = x18*xi
            x20 = 467775.0d0*x17
            x21 = x20*(46410.0d0*x11 - 62985.0d0*x14 + 29393.0d0*x16 + 1365.0d0*x5 - 13650.0d0*x8 - 21.0d0)
            x22 = 819.0d0*x4
            x23 = x1*x22
            x24 = x13*x9
            x25 = 151164.0d0*x24
            x26 = x10*x6
            x27 = 83538.0d0*x26
            x28 = 1092.0d0*x5
            x29 = 88179.0d0*x15
            x30 = x12*x29
            x31 = x1*x30
            x32 = x28 + x31
            x33 = x0*x7
            x34 = 16380.0d0*x33
            x35 = -x1*x34 - 63.0d0
            x36 = -4199.0d0*x14
            x37 = 7956.0d0*x11 + x36 - 4914.0d0*x8
            x38 = x3**(-7.5d0)*xi*yi*zi
            x39 = 6081075.0d0*x38
            x40 = x2*x22
            x41 = x2*x30
            x42 = x28 + x41
            x43 = -x2*x34 - 63.0d0
            x44 = 91.0d0*x4
            x45 = x1*x44
            x46 = 5460.0d0*x33
            x47 = -x1*x46
            x48 = x47 - 21.0d0
            x49 = 117572.0d0*x24
            x50 = 46410.0d0*x26
            x51 = -x1*x49 + x1*x50
            x52 = 18564.0d0*x11 - 12597.0d0*x14 - 8190.0d0*x8
            x53 = x18*yi
            x54 = x45 - 7.0d0
            x55 = 6188.0d0*x11 + x36 + 364.0d0*x5 - 2730.0d0*x8
            x56 = x18*zi
            x57 = x2*x44
            x58 = x57 - 7.0d0
            x59 = -x2*x46
            x60 = -x2*x49 + x2*x50 + x59
            x61 = x57 - 21.0d0
            x62 = yi**4
            x63 = x0*x62
            x64 = 7735.0d0*x10
            x65 = x1*x26
            x66 = x13*x6
            x67 = 29393.0d0*x66
            x68 = x1*x33
            x69 = -x62*x67 - 2730.0d0*x68
            x70 = x1*x4
            x71 = 182.0d0*x70
            x72 = 455.0d0*x7
            x73 = -x62*x72 + x71
            x74 = 221.0d0*x11
            x75 = x1*x24
            x76 = 29393.0d0*x15
            x77 = x76*x9
            x78 = x62*x77 + x74 - 8398.0d0*x75
            x79 = 91.0d0*x5
            x80 = x79 - 273.0d0*x8
            x81 = x80 - 7.0d0
            x82 = x20*xi
            x83 = 595.0d0*x68 + 7.0d0
            x84 = -323.0d0*x11 - 105.0d0*x5 + 357.0d0*x8
            x85 = 4199.0d0*x13
            x86 = x85*x9
            x87 = -x2*x86
            x88 = x2*x26
            x89 = x1*x2
            x90 = x66*x89
            x91 = x29*x9
            x92 = x89*x91
            x93 = -x1*x86 + x92
            x95 = x1*x7
            x96 = x0*x89
            x97 = x10*x96
            x98 = -1365.0d0*x2*x33 - 1365.0d0*x2*x95 + x57 - 1365.0d0*x68 + 23205.0d0*x97
            x99 = x54 + x98
            x100 = x2*x4
            x101 = x2*x24
            x102 = x2*x33
            x103 = 595.0d0*x102 + 7.0d0
            x104 = zi**4
            x105 = x0*x104
            x106 = -2730.0d0*x102 - x104*x67
            x107 = 182.0d0*x100 - x104*x72
            x108 = -8398.0d0*x101 + x104*x77 + x74
            x109 = 9945.0d0*x10
            x110 = x109*x63
            x111 = x62*x66
            x112 = -62985.0d0*x111 - 5850.0d0*x68 - 15.0d0
            x113 = x62*x7
            x114 = 195.0d0*x113
            x115 = -x114
            x116 = x115 + 130.0d0*x70
            x117 = 3315.0d0*x11 + 585.0d0*x5 - 2925.0d0*x8
            x118 = 3315.0d0*x10
            x119 = x118*x63
            x120 = 6630.0d0*x65
            x121 = -195.0d0*x8
            x122 = 26.0d0*x70
            x123 = x121 + x122
            x124 = 39.0d0*x5
            x125 = x124 - 1.0d0
            x126 = -20995.0d0*x111
            x127 = -1170.0d0*x68
            x128 = x126 + x127
            x129 = x20*zi
            x130 = 9945.0d0*x88
            x131 = x118*x6
            x132 = x1*x131
            x133 = 12597.0d0*x13
            x134 = x133*x9
            x135 = 39.0d0*x100
            x136 = x109*x96
            x137 = x136 - 3.0d0
            x138 = x135 + x137 - 585.0d0*x68
            x139 = 13.0d0*x70
            x140 = -1755.0d0*x102 + x139
            x141 = x2*x95
            x142 = 663.0d0*x11 - 195.0d0*x141 + 117.0d0*x5 - 585.0d0*x8 - 62985.0d0*x90
            x143 = 9945.0d0*x65
            x144 = x131*x2
            x145 = 39.0d0*x70
            x146 = -585.0d0*x102 + x145
            x147 = 13.0d0*x100
            x148 = x147 - 3.0d0
            x149 = x136 + x148 - 1755.0d0*x68
            x150 = x104*x7
            x151 = 65.0d0*x150
            x152 = x105*x118
            x153 = 6630.0d0*x88
            x154 = x104*x66
            x155 = -20995.0d0*x154
            x156 = -1170.0d0*x102
            x157 = x155 + x156
            x158 = 26.0d0*x100
            x159 = x121 + x158
            x160 = x20*yi
            x161 = x105*x109
            x162 = -5850.0d0*x102 - 62985.0d0*x154 - 15.0d0
            x163 = -195.0d0*x150
            x164 = 130.0d0*x100 + x163
            x165 = yi**6
            x166 = x0*x13
            x167 = 41990.0d0*x166
            x168 = x10*x63
            x169 = x29*x6
            x170 = x121 + 130.0d0*x5
            x171 = -2925.0d0*x113 + x118*x165 + 585.0d0*x70
            x172 = -150.0d0*x5 + 255.0d0*x8 + 15.0d0
            x173 = 2027025.0d0*x38
            x174 = x2*x62
            x175 = x118*x174
            x176 = 26.0d0*x5
            x177 = 13260.0d0*x97
            x178 = x177 - 3.0d0
            x179 = -1170.0d0*x141 + x176 + x178 - 39.0d0*x8 - 25194.0d0*x90
            x180 = x6*x85
            x181 = x135 + x169*x174 - x180*x62 - 780.0d0*x68
            x182 = 41990.0d0*x13
            x183 = x2*x63
            x184 = -390.0d0*x102 - x182*x183 + 78.0d0*x70
            x185 = x115 + x184
            x186 = -3230.0d0*x97
            x187 = -45.0d0*x100 + x186 + 510.0d0*x68 + 9.0d0
            x188 = 510.0d0*x102 - 45.0d0*x70
            x189 = x1*x104
            x190 = x118*x189
            x191 = x10*x105
            x192 = -780.0d0*x102 - x104*x180 + x169*x189
            x193 = x1*x105
            x194 = 78.0d0*x100 + x163 - x182*x193 - 390.0d0*x68
            x195 = zi**6
            x196 = 585.0d0*x100 + x118*x195 - 2925.0d0*x150
            x197 = 273.0d0*x113
            x198 = 221.0d0*x10
            x199 = x165*x198
            x200 = -x197 + x199 + x54
            x201 = 8398.0d0*x166
            x202 = x6*x76
            x203 = -x165*x201 + x165*x202
            x204 = 182.0d0*x5 - 455.0d0*x8
            x205 = x176 - 65.0d0*x8
            x206 = x115 + 6630.0d0*x168
            x207 = x132 + x145
            x208 = x10*x174
            x209 = -390.0d0*x141 + 78.0d0*x5 - 41990.0d0*x90
            x210 = x178 + x209
            x211 = 25194.0d0*x13
            x212 = -39.0d0*x113 + x156 - x183*x211
            x213 = 146965.0d0*x15
            x214 = x213*x6
            x215 = 93555.0d0*x17
            x216 = 6630.0d0*x191
            x217 = x10*x189
            x218 = 39.0d0*x150
            x219 = x127 - x193*x211 - x218
            x220 = x195*x198
            x221 = x135 + x163 + x220 - 1.0d0
            x222 = -x195*x201 + x195*x202
            x223 = -273.0d0*x150 + x220
            x224 = x223 + x58
            x225 = 1092.0d0*x70
            x226 = yi**8
            x227 = x0*x29
            x228 = x226*x227
            x229 = x225 + x228
            x230 = x165*x166
            x231 = 46410.0d0*x168 - 117572.0d0*x230 + x79
            x232 = x10*x165
            x233 = -8190.0d0*x113 - x133*x226 + 18564.0d0*x232
            x234 = 323.0d0*x232
            x235 = -35.0d0*x5
            x236 = x109*x174
            x237 = x165*x2
            x238 = 13.0d0*x5
            x239 = -1755.0d0*x141 + x238
            x240 = x0*x85
            x241 = x227*x237
            x242 = -x165*x240 + x241
            x243 = x13*x183
            x244 = -195.0d0*x102 - 585.0d0*x113 + 663.0d0*x232 - 62985.0d0*x243 + 117.0d0*x70
            x245 = 510.0d0*x141 - 45.0d0*x5
            x246 = x104*x62
            x247 = x13*x246
            x248 = x104*x63
            x249 = x178 + x194
            x250 = 90.0d0*x100
            x251 = x13*x193
            x252 = x109*x189
            x253 = x133*x195
            x254 = x1*x195
            x255 = x227*x254
            x256 = -x195*x240 + x255
            x257 = x10*x195
            x258 = 117.0d0*x100 + x137 - 585.0d0*x150 - 62985.0d0*x251 + 663.0d0*x257 - 195.0d0*x68
            x259 = 323.0d0*x257
            x260 = 105.0d0*x100
            x261 = 357.0d0*x150
            x262 = x166*x195
            x263 = 1092.0d0*x100
            x264 = zi**8
            x265 = x227*x264
            x266 = x263 + x265
            x267 = 46410.0d0*x191 - 117572.0d0*x262 + x59 + x79
            x268 = -x133*x264 - 8190.0d0*x150 + 18564.0d0*x257 - 21.0d0
            x269 = 819.0d0*x5
            x270 = 4914.0d0*x113
            x271 = x226*x85
            x272 = -x271
            x273 = 7956.0d0*x232
            x274 = -x270 + x272 + x273
            x275 = -x237*x85
            x276 = x124 - 585.0d0*x141
            x277 = x124 - 780.0d0*x141 - x246*x85 + x248*x29
            x278 = -x254*x85
            x279 = x264*x85
            x280 = -x279
            x281 = 2730.0d0*x150
            x282 = 364.0d0*x100
            x283 = 6188.0d0*x257
            x284 = -4914.0d0*x150 + 7956.0d0*x257 + x280
            x285 = x13*x226
            x286 = yi**10
            x287 = -13650.0d0*x113 + 46410.0d0*x232 - 62985.0d0*x285 + x286*x76 + 1365.0d0*x70 - 21.0d0
            x288 = 8.0d0*x70
            x289 = x147 - 1.0d0
            x290 = x13*x237
            x294 = 264.0d0*x148
            x295 = -x158
            x296 = x151 + x295 + 1.0d0
            x297 = 11.0d0*x100
            x298 = x70*(x297 - 1.0d0)
            x299 = 240.0d0*x150
            x300 = 2640.0d0*x141
            x301 = x17*xi
            x303 = 10.0d0*x100
            x306 = 80.0d0*x150
            x307 = -x306
            x308 = 32.0d0*x70
            x309 = x218 + x295 + 3.0d0
            x311 = x297 - 3.0d0
            x313 = 4.0d0*x100
            x315 = 143.0d0*x150
            x316 = x70*(-66.0d0*x100 + x315 + 3.0d0)
            x317 = 16.0d0*x100
            x318 = 16.0d0*x150
            x319 = 9.0d0*x100 - 1.0d0
            x320 = 24.0d0*x100
            x321 = 64.0d0*x257
            x322 = x13*x254
            x323 = -24.0d0*x150 + 16.0d0*x257 + x303 - 1.0d0
            x324 = x307 + x320 + x321 - 1.0d0
            x327 = 2.0d0*x100 - 1.0d0
            x328 = x313 - 1.0d0
            x329 = -x317 + x318 + 3.0d0
            x330 = x70*(-12.0d0*x100 + x318 + 1.0d0)
            x331 = x13*x264
            x332 = zi**10
            x333 = 1365.0d0*x100 - 13650.0d0*x150 + 46410.0d0*x257 - 62985.0d0*x331 + x332*x76 - 21.0d0
            x334 = x2*x226*x29 + x225
            x335 = -16380.0d0*x141 - 63.0d0
            x336 = -5460.0d0*x141
            x337 = x104*x165
            x338 = -2730.0d0*x141 - 29393.0d0*x247
            x339 = -5850.0d0*x141 - 62985.0d0*x247 - 15.0d0
            x340 = x195*x62
            x341 = x1*x264*x29 + x263
            T(i, 1) = -x19*(218790.0d0*x11 - 230945.0d0*x14 + 88179.0d0*x16 + 15015.0d0*x5 - 90090.0d0*x8 - 693.0d0)
            T(i, 2) = -x21*yi
            T(i, 3) = -x21*zi
            T(i, 4) = -x19*(-x1*x25 + x1*x27 + x23 + x32 + x35 + x37)
            T(i, 5) = -x39*(-3876.0d0*x11 + 2261.0d0*x14 - 420.0d0*x5 + 2142.0d0*x8 + 21.0d0)
            T(i, 6) = -x19*(-x2*x25 + x2*x27 + x37 + x40 + x42 + x43)
            T(i, 7) = -x53*(x32 + x45 + x48 + x51 + x52)
            T(i, 8) = -x56*(x31 + x47 + x51 + x54 + x55)
            T(i, 9) = -x53*(x41 + x55 + x58 + x60)
            T(i, 10) = -x56*(x42 + x52 + x60 + x61)
            T(i, 11) = -x82*(x63*x64 + 9282.0d0*x65 + x69 + x73 + x78 + x81)
            T(i, 12) = -x39*(-2261.0d0*x65 - 35.0d0*x70 + 2261.0d0*x75 + x83 + x84)
            T(i, 13) = -x19*(4641.0d0*x65 + x74 + x80 + x87 + 4641.0d0*x88 - 88179.0d0*x90 + x93 + x99)
            T(i, 14) = -x39*(-35.0d0*x100 + 2261.0d0*x101 + x103 + x84 - 2261.0d0*x88)
            T(i, 15) = -x82*(x105*x64 + x106 + x107 + x108 + x81 + 9282.0d0*x88)
            T(i, 16) = -x53*(x110 + x112 + x116 + x117 + x62*x91 + 33150.0d0*x65 - 41990.0d0*x75)
            T(i, 17) = -x129*(-65.0d0*x113 + x119 + x120 + x123 + x125 + x128 + x78)
            T(i, 18) = -x53*(x130 + x132 - x134*x2 + x138 + x140 + x142 + x93)
            T(i, 19) = -x56*(-x1*x134 + x142 + x143 + x144 + x146 + x149 + x87 + x92)
            T(i, 20) = -x160*(x108 + x125 - x151 + x152 + x153 + x157 + x159)
            T(i, 21) = -x56*(-41990.0d0*x101 + x104*x91 + x117 + x161 + x162 + x164 + 33150.0d0*x88)
            T(i, 22) = -x19*(x112 + x143 - x165*x167 + x165*x169 + 33150.0d0*x168 + x170 + x171)
            T(i, 23) = -x173*(6783.0d0*x111 + 255.0d0*x113 - 3230.0d0*x168 + x172 - 3230.0d0*x65 + 1700.0d0*x68 - 150.0d0*x70)
            T(i, 24) = -x19*(2210.0d0*x168 + x175 + x179 + x181 + x185 + 1326.0d0*x65 + 663.0d0*x88)
            T(i, 25) = -x173*(255.0d0*x141 + x187 + x188 - 90.0d0*x5 - 969.0d0*x65 + 153.0d0*x8 - 969.0d0*x88 + 6783.0d0*x90)
            T(i, 26) = -x19*(x145 + x179 + x190 + 2210.0d0*x191 + x192 + x194 + 663.0d0*x65 + 1326.0d0*x88)
            T(i, 27) = -x173*(-150.0d0*x100 + 1700.0d0*x102 + 255.0d0*x150 + 6783.0d0*x154 + x172 - 3230.0d0*x191 - 3230.0d0*x88)
            T(i, 28) = -x19*(x130 + x162 - x167*x195 + x169*x195 + x170 + 33150.0d0*x191 + x196)
            T(i, 29) = -x160*(9282.0d0*x168 + x200 + x203 + x204 + 7735.0d0*x65 + x69)
            T(i, 30) = -x129*(x128 + x199 + x203 + x205 + x206 + x207 - 1.0d0)
            T(i, 31) = -x53*(x123 + x144 + 1326.0d0*x168 + x181 + 663.0d0*x208 + x210 + x212 + 2210.0d0*x65)
            T(i, 32) = -x215*zi*(x120 + x121 + x126 + x148 + x174*x214 + x177 + x184 + x206 + 1105.0d0*x208 + x209 - 2340.0d0*x68 + 1105.0d0*x88)
            T(i, 33) = -x215*yi*(-2340.0d0*x102 + x121 + x139 + x153 + x155 + x189*x214 + x194 + x210 + x216 + 1105.0d0*x217 + 1105.0d0*x65)
            T(i, 34) = -x56*(x159 + 1326.0d0*x191 + x192 + x207 + x210 + 663.0d0*x217 + x219 + 2210.0d0*x88)
            T(i, 35) = -x160*(x144 + x157 + x205 + x216 + x221 + x222)
            T(i, 36) = -x129*(x106 + 9282.0d0*x191 + x204 + x222 + x224 + 7735.0d0*x88)
            T(i, 37) = -x19*(x229 + x231 + x233 + x48)
            T(i, 38) = -x39*(357.0d0*x113 - 2261.0d0*x168 + 2261.0d0*x230 - x234 + x235 - 105.0d0*x70 + x83)
            T(i, 39) = -x19*(x119 - x133*x237 + x138 + x236 + x239 + x242 + x244)
            T(i, 40) = -x173*(255.0d0*x102 + 153.0d0*x113 - 969.0d0*x168 + x187 - 969.0d0*x208 + 6783.0d0*x243 + x245 - 90.0d0*x70)
            T(i, 41) = -x215*xi*(-2340.0d0*x141 + 1105.0d0*x168 + x185 + 1105.0d0*x191 + 6630.0d0*x208 + x213*x248 + 6630.0d0*x217 + x238 - 20995.0d0*x247 + x249)
            T(i, 42) = -x173*(153.0d0*x150 + x186 + x188 - 969.0d0*x191 - 969.0d0*x217 + x245 - x250 + 6783.0d0*x251 + 255.0d0*x68 + 9.0d0)
            T(i, 43) = -x19*(-x1*x253 + x146 + x152 + x239 + x252 + x256 + x258)
            T(i, 44) = -x39*(x103 - 2261.0d0*x191 + x235 - x259 - x260 + x261 + 2261.0d0*x262)
            T(i, 45) = -x19*(x266 + x267 + x268)
            T(i, 46) = -x53*(83538.0d0*x168 + x229 - 151164.0d0*x230 + x269 + x274 + x35)
            T(i, 47) = -x56*(-2730.0d0*x113 + x228 + x231 + 6188.0d0*x232 + x272 + x47 + 364.0d0*x70 - 7.0d0)
            T(i, 48) = -x53*(4641.0d0*x168 + x200 + 4641.0d0*x208 + x242 - 88179.0d0*x243 + x275 + x79 + x98)
            T(i, 49) = -x56*(-x0*x133*x165 + x110 + x149 + x175 + x241 + x244 + x275 + x276)
            T(i, 50) = -x53*(x122 + x152 + 663.0d0*x168 + 1326.0d0*x208 + x212 + 2210.0d0*x217 + x249 + x277)
            T(i, 51) = -x56*(x119 + x158 + x178 + x185 + 663.0d0*x191 + 2210.0d0*x208 + 1326.0d0*x217 + x219 + x277)
            T(i, 52) = -x53*(-x0*x253 + x140 + x161 + x190 + x255 + x258 + x276 + x278)
            T(i, 53) = -x56*(4641.0d0*x191 + 4641.0d0*x217 + x223 - 88179.0d0*x251 + x256 + x278 + x79 + x99)
            T(i, 54) = -x53*(x265 + x267 + x280 - x281 + x282 + x283 - 7.0d0)
            T(i, 55) = -x56*(83538.0d0*x191 - 151164.0d0*x262 + x266 + x269 + x284 + x43)
            T(i, 56) = -x287*x82
            T(i, 57) = -x173*(-x225 + x270 + x271 - x273 + x288*(-459.0d0*x113 + x234 + 189.0d0*x70 - 21.0d0) + 63.0d0)
            T(i, 58) = -x19*(176.0d0*x141*x200 + x288*(-28.0d0*x100 + 546.0d0*x141 + x197 - x199 - 2184.0d0*x208 + 2210.0d0*x290 - x45 + 7.0d0) + x289*(2002.0d0*x113 - 4004.0d0*x232 + 2431.0d0*x285 - 308.0d0*x70 + 7.0d0))
            T(i, 59) = -93555.0d0*x38*(x288*(70.0d0*x100 - 390.0d0*x113 - 728.0d0*x141 + 1300.0d0*x208 + 273.0d0*x70 - 35.0d0) + 22.0d0*x289*x70*(x114 - x71 + 35.0d0) + 22.0d0*x70*(140.0d0*x100 + x115 - 1092.0d0*x141 + 1560.0d0*x208 + x71 - 35.0d0) + 13.0d0*(5.0d0*x100 - 1.0d0)*(-1001.0d0*x113 + 715.0d0*x232 + 385.0d0*x70 - 35.0d0))
            T(i, 60) = -8505.0d0*x301*(x141*x294*(143.0d0*x113 - 110.0d0*x70 + 15.0d0) + 96.0d0*x141*(-858.0d0*x113 - 1320.0d0*x141 + 2860.0d0*x208 + x250 + 495.0d0*x70 - 45.0d0) + x288*(-180.0d0*x100 + 286.0d0*x113 - 5720.0d0*x208 - 4400.0d0*x217 + 11440.0d0*x247 + x299 + x300 - 165.0d0*x70 + 15.0d0) + 11.0d0*x296*(-495.0d0*x113 + 429.0d0*x232 + 135.0d0*x70 - 5.0d0) + 36.0d0*x298*(60.0d0*x100 - 143.0d0*x113 - 660.0d0*x141 + 1144.0d0*x208 + 110.0d0*x70 - 15.0d0))
            T(i, 61) = -14175.0d0*x38*(x294*x70*(-20.0d0*x100 + 66.0d0*x141 - 11.0d0*x70 + 5.0d0) + 396.0d0*x296*x70*(11.0d0*x70 - 5.0d0) + 144.0d0*x298*(-20.0d0*x100 + 88.0d0*x141 - 33.0d0*x70 + 10.0d0) + x308*(80.0d0*x100 - 440.0d0*x141 + 528.0d0*x217 + x307 + 66.0d0*x70 - 15.0d0) + 72.0d0*x70*(120.0d0*x100 - 528.0d0*x141 - 160.0d0*x150 + 880.0d0*x217 + 33.0d0*x70 - 10.0d0) + 429.0d0*(33.0d0*x113 - 30.0d0*x70 + 5.0d0)*(17.0d0*x150 - x303 + 1.0d0))
            T(i, 62) = -14175.0d0*x301*(1848.0d0*x141*x309*(3.0d0*x70 - 1.0d0) + 672.0d0*x141*x311*(24.0d0*x141 - x313 - 9.0d0*x70 + 2.0d0) + 448.0d0*x141*(-120.0d0*x141 + 144.0d0*x217 + x317 - x318 + 18.0d0*x70 - 3.0d0) + 33.0d0*x221*(21.0d0*x113 - 14.0d0*x70 + 1.0d0) + x308*(180.0d0*x141 - 720.0d0*x217 + x306 - x320 - x321 + 672.0d0*x322 - 6.0d0*x70 + 1.0d0) + 84.0d0*x316*(18.0d0*x141 - x313 - 3.0d0*x70 + 1.0d0) + 56.0d0*x319*x70*(-144.0d0*x141 - 32.0d0*x150 + 240.0d0*x217 + x320 + 9.0d0*x70 - 2.0d0))
            T(i, 63) = -4725.0d0*x38*(3234.0d0*x221*x70 + 3234.0d0*x309*x328*x70 + 1176.0d0*x311*x330 + 1176.0d0*x316*x327 + 784.0d0*x319*x329*x70 + 896.0d0*x323*x70 + 784.0d0*x324*x70 + 429.0d0*(7.0d0*x70 - 3.0d0)*(x259 + x260 - x261 - 7.0d0))
            T(i, 64) = -4725.0d0*x301*(5120.0d0*x141*x323 + 1344.0d0*x141*x327*(-110.0d0*x100 + x315 + 15.0d0) + 4480.0d0*x141*x329*(3.0d0*x100 - 1.0d0) + x224*x300 + 320.0d0*x324*x70*(7.0d0*x100 - 1.0d0) + 840.0d0*x328*x70*(33.0d0*x100 + 143.0d0*x257 - x315 - 1.0d0) + 560.0d0*x330*(-18.0d0*x100 + 33.0d0*x150 + 1.0d0) + 128.0d0*x70*(-40.0d0*x100 - 448.0d0*x257 + x299 + 256.0d0*x331 + 1.0d0) + 33.0d0*(5.0d0*x70 - 1.0d0)*(x279 + x281 - x282 - x283 + 7.0d0))
            T(i, 65) = -x39*(-420.0d0*x100 + 2142.0d0*x150 - 3876.0d0*x257 + 2261.0d0*x331 + 21.0d0)
            T(i, 66) = -x333*x82
            T(i, 67) = -x53*(-90090.0d0*x113 + 218790.0d0*x232 - 230945.0d0*x285 + x286*x29 + 15015.0d0*x70 - 693.0d0)
            T(i, 68) = -x129*x287
            T(i, 69) = -x53*(83538.0d0*x208 + x274 - 151164.0d0*x290 + x334 + x335 + x40)
            T(i, 70) = -x56*(46410.0d0*x208 + x233 - 117572.0d0*x290 + x334 + x336 + x61)
            T(i, 71) = -x160*(x107 + x189*x64 + x200 + 9282.0d0*x208 - 8398.0d0*x290 + x337*x76 + x338)
            T(i, 72) = -x56*(x164 + x171 + 33150.0d0*x208 + x252 + x29*x337 - 41990.0d0*x290 + x339)
            T(i, 73) = -x53*(x116 + x196 + 33150.0d0*x217 + x236 + x29*x340 - 41990.0d0*x322 + x339)
            T(i, 74) = -x129*(x174*x64 + 9282.0d0*x217 + x224 - 8398.0d0*x322 + x338 + x340*x76 + x73)
            T(i, 75) = -x53*(46410.0d0*x217 + x268 - 117572.0d0*x322 + x336 + x341 + x45)
            T(i, 76) = -x56*(83538.0d0*x217 + x23 + x284 - 151164.0d0*x322 + x335 + x341)
            T(i, 77) = -x160*x333
            T(i, 78) = -x56*(15015.0d0*x100 - 90090.0d0*x150 + 218790.0d0*x257 + x29*x332 - 230945.0d0*x331 - 693.0d0)
        end do
    end subroutine T11
    subroutine T12(x, y, z, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), 455)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: x0
        real(8) :: x1
        real(8) :: x2
        real(8) :: x3
        real(8) :: x4
        real(8) :: x5
        real(8) :: x6
        real(8) :: x7
        real(8) :: x8
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x13
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
        real(8) :: x27
        real(8) :: x28
        real(8) :: x29
        real(8) :: x30
        real(8) :: x31
        real(8) :: x32
        real(8) :: x33
        real(8) :: x34
        real(8) :: x35
        real(8) :: x36
        real(8) :: x37
        real(8) :: x38
        real(8) :: x39
        real(8) :: x40
        real(8) :: x41
        real(8) :: x42
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
        real(8) :: x53
        real(8) :: x54
        real(8) :: x55
        real(8) :: x56
        real(8) :: x57
        real(8) :: x58
        real(8) :: x59
        real(8) :: x60
        real(8) :: x61
        real(8) :: x62
        real(8) :: x63
        real(8) :: x64
        real(8) :: x65
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
        real(8) :: x80
        real(8) :: x81
        real(8) :: x82
        real(8) :: x83
        real(8) :: x84
        real(8) :: x85
        real(8) :: x86
        real(8) :: x87
        real(8) :: x88
        real(8) :: x89
        real(8) :: x90
        real(8) :: x91
        real(8) :: x92
        real(8) :: x93
        real(8) :: x94
        real(8) :: x95
        real(8) :: x96
        real(8) :: x97
        real(8) :: x98
        real(8) :: x99
        real(8) :: x100
        real(8) :: x101
        real(8) :: x102
        real(8) :: x103
        real(8) :: x104
        real(8) :: x105
        real(8) :: x106
        real(8) :: x107
        real(8) :: x108
        real(8) :: x109
        real(8) :: x110
        real(8) :: x111
        real(8) :: x112
        real(8) :: x113
        real(8) :: x114
        real(8) :: x115
        real(8) :: x116
        real(8) :: x117
        real(8) :: x118
        real(8) :: x119
        real(8) :: x120
        real(8) :: x121
        real(8) :: x122
        real(8) :: x123
        real(8) :: x124
        real(8) :: x125
        real(8) :: x126
        real(8) :: x127
        real(8) :: x128
        real(8) :: x129
        real(8) :: x130
        real(8) :: x131
        real(8) :: x132
        real(8) :: x133
        real(8) :: x134
        real(8) :: x135
        real(8) :: x136
        real(8) :: x137
        real(8) :: x138
        real(8) :: x139
        real(8) :: x140
        real(8) :: x141
        real(8) :: x142
        real(8) :: x143
        real(8) :: x144
        real(8) :: x145
        real(8) :: x146
        real(8) :: x147
        real(8) :: x148
        real(8) :: x149
        real(8) :: x150
        real(8) :: x151
        real(8) :: x152
        real(8) :: x153
        real(8) :: x154
        real(8) :: x155
        real(8) :: x156
        real(8) :: x157
        real(8) :: x158
        real(8) :: x159
        real(8) :: x160
        real(8) :: x161
        real(8) :: x162
        real(8) :: x163
        real(8) :: x164
        real(8) :: x165
        real(8) :: x166
        real(8) :: x167
        real(8) :: x168
        real(8) :: x169
        real(8) :: x170
        real(8) :: x171
        real(8) :: x172
        real(8) :: x173
        real(8) :: x174
        real(8) :: x175
        real(8) :: x176
        real(8) :: x177
        real(8) :: x178
        real(8) :: x179
        real(8) :: x180
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
        real(8) :: x196
        real(8) :: x197
        real(8) :: x198
        real(8) :: x199
        real(8) :: x200
        real(8) :: x201
        real(8) :: x202
        real(8) :: x203
        real(8) :: x204
        real(8) :: x205
        real(8) :: x206
        real(8) :: x207
        real(8) :: x208
        real(8) :: x209
        real(8) :: x210
        real(8) :: x211
        real(8) :: x212
        real(8) :: x213
        real(8) :: x214
        real(8) :: x215
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
        real(8) :: x227
        real(8) :: x228
        real(8) :: x229
        real(8) :: x230
        real(8) :: x231
        real(8) :: x232
        real(8) :: x233
        real(8) :: x234
        real(8) :: x235
        real(8) :: x236
        real(8) :: x237
        real(8) :: x238
        real(8) :: x239
        real(8) :: x240
        real(8) :: x241
        real(8) :: x242
        real(8) :: x243
        real(8) :: x244
        real(8) :: x245
        real(8) :: x246
        real(8) :: x247
        real(8) :: x248
        real(8) :: x249
        real(8) :: x250
        real(8) :: x251
        real(8) :: x252
        real(8) :: x253
        real(8) :: x254
        real(8) :: x255
        real(8) :: x256
        real(8) :: x257
        real(8) :: x258
        real(8) :: x259
        real(8) :: x260
        real(8) :: x261
        real(8) :: x262
        real(8) :: x263
        real(8) :: x264
        real(8) :: x265
        real(8) :: x266
        real(8) :: x267
        real(8) :: x268
        real(8) :: x269
        real(8) :: x270
        real(8) :: x271
        real(8) :: x272
        real(8) :: x273
        real(8) :: x274
        real(8) :: x275
        real(8) :: x276
        real(8) :: x277
        real(8) :: x278
        real(8) :: x279
        real(8) :: x280
        real(8) :: x281
        real(8) :: x282
        real(8) :: x283
        real(8) :: x284
        real(8) :: x285
        real(8) :: x286
        real(8) :: x287
        real(8) :: x288
        real(8) :: x289
        real(8) :: x290
        real(8) :: x291
        real(8) :: x292
        real(8) :: x293
        real(8) :: x294
        real(8) :: x295
        real(8) :: x296
        real(8) :: x298
        real(8) :: x299
        real(8) :: x300
        real(8) :: x301
        real(8) :: x302
        real(8) :: x303
        real(8) :: x304
        real(8) :: x305
        real(8) :: x306
        real(8) :: x307
        real(8) :: x308
        real(8) :: x309
        real(8) :: x310
        real(8) :: x311
        real(8) :: x312
        real(8) :: x313
        real(8) :: x314
        real(8) :: x315
        real(8) :: x316
        real(8) :: x317
        real(8) :: x318
        real(8) :: x319
        real(8) :: x320
        real(8) :: x321
        real(8) :: x322
        real(8) :: x323
        real(8) :: x324
        real(8) :: x325
        real(8) :: x326
        real(8) :: x327
        real(8) :: x328
        real(8) :: x329
        real(8) :: x330
        real(8) :: x331
        real(8) :: x332
        real(8) :: x333
        real(8) :: x334
        real(8) :: x335
        real(8) :: x336
        real(8) :: x337
        real(8) :: x338
        real(8) :: x339
        real(8) :: x340
        real(8) :: x341
        real(8) :: x342
        real(8) :: x343
        real(8) :: x344
        real(8) :: x345
        real(8) :: x346
        real(8) :: x347
        real(8) :: x348
        real(8) :: x349
        real(8) :: x350
        real(8) :: x351
        real(8) :: x352
        real(8) :: x353
        real(8) :: x354
        real(8) :: x355
        real(8) :: x356
        real(8) :: x357
        real(8) :: x358
        real(8) :: x359
        real(8) :: x360
        real(8) :: x361
        real(8) :: x362
        real(8) :: x363
        real(8) :: x364
        real(8) :: x365
        real(8) :: x366
        real(8) :: x367
        real(8) :: x368
        real(8) :: x369
        real(8) :: x370
        real(8) :: x371
        real(8) :: x372
        real(8) :: x373
        real(8) :: x374
        real(8) :: x375
        real(8) :: x376
        real(8) :: x377
        real(8) :: x378
        real(8) :: x379
        real(8) :: x380
        real(8) :: x381
        real(8) :: x382
        real(8) :: x383
        real(8) :: x384
        real(8) :: x385
        real(8) :: x386
        real(8) :: x387
        real(8) :: x390
        real(8) :: x391
        real(8) :: x395
        real(8) :: x396
        real(8) :: x397
        real(8) :: x399
        real(8) :: x400
        real(8) :: x403
        real(8) :: x404
        real(8) :: x405
        real(8) :: x406
        real(8) :: x407
        real(8) :: x408
        real(8) :: x409
        real(8) :: x410
        real(8) :: x411
        real(8) :: x412
        real(8) :: x415
        real(8) :: x416
        real(8) :: x417
        real(8) :: x418
        real(8) :: x421
        real(8) :: x423
        real(8) :: x424
        real(8) :: x426
        real(8) :: x427
        real(8) :: x430
        real(8) :: x431
        real(8) :: x432
        real(8) :: x434
        real(8) :: x435
        real(8) :: x437
        real(8) :: x438
        real(8) :: x439
        real(8) :: x440
        real(8) :: x441
        real(8) :: x442
        real(8) :: x443
        real(8) :: x444
        real(8) :: x445
        real(8) :: x446
        real(8) :: x447
        real(8) :: x448
        real(8) :: x449
        real(8) :: x450
        real(8) :: x451
        do concurrent(i=1:size(x))
            xi = x(i)
            yi = y(i)
            zi = z(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0/x3
            x5 = x0*x4
            x6 = xi**4
            x7 = x3**(-2)
            x8 = x6*x7
            x9 = xi**6
            x10 = x3**(-3)
            x11 = x10*x9
            x12 = xi**8
            x13 = x3**(-4)
            x14 = x12*x13
            x15 = xi**10
            x16 = x3**(-5)
            x17 = x15*x16
            x18 = x3**(-6)
            x19 = 676039.0d0*x18
            x20 = x3**(-6.5d0)
            x21 = 467775.0d0*x20
            x22 = x3**(-7.5d0)
            x23 = x22*yi
            x24 = 6081075.0d0*xi
            x25 = x24*(106590.0d0*x11 - 124355.0d0*x14 + 52003.0d0*x17 + 5775.0d0*x5 - 39270.0d0*x8 - 231.0d0)
            x26 = x22*zi
            x27 = 273.0d0*x4
            x28 = x1*x27
            x29 = -x28
            x30 = x1*x12
            x31 = 1322685.0d0*x16
            x32 = x10*x6
            x33 = 232050.0d0*x32
            x34 = x15*x19
            x35 = x13*x9
            x36 = 881790.0d0*x1
            x37 = x0*x7
            x38 = 20475.0d0*x37
            x39 = x1*x38 + 21.0d0
            x40 = -46410.0d0*x11 + 62985.0d0*x14 - 29393.0d0*x17 - 1365.0d0*x5 + 13650.0d0*x8
            x41 = 225.0d0*x5
            x42 = x23*zi
            x43 = 42567525.0d0*x42
            x44 = -x2*x27
            x45 = x12*x2
            x46 = x2*x35
            x47 = x2*x38 + 21.0d0
            x48 = 7140.0d0*x37
            x49 = -x1*x48
            x50 = x49 - 63.0d0
            x51 = 315.0d0*x4
            x52 = x1*x51
            x53 = x1*x35
            x54 = 40698.0d0*x32
            x55 = 52003.0d0*x16
            x56 = x1*x54 + x30*x55 + x52 - 81396.0d0*x53
            x57 = 11628.0d0*x11 - 6783.0d0*x14 + 1260.0d0*x5 - 6426.0d0*x8
            x58 = x23*x24
            x59 = x49 - 21.0d0
            x60 = 3876.0d0*x11 - 2261.0d0*x14 + 420.0d0*x5 - 2142.0d0*x8
            x61 = x24*x26
            x62 = -x2*x48
            x63 = x62 - 21.0d0
            x64 = x2*x51
            x65 = x2*x54 + x45*x55 - 81396.0d0*x46 + x64
            x66 = x62 - 63.0d0
            x67 = yi**4
            x68 = x67*x9
            x69 = 823004.0d0*x16
            x70 = x16*x30
            x71 = x1*x32
            x72 = x0*x10
            x73 = 30940.0d0*x72
            x74 = x12*x19
            x75 = x1*x4
            x76 = 182.0d0*x75
            x77 = -x76
            x78 = 455.0d0*x7
            x79 = x67*x78 + x77
            x80 = x1*x37
            x81 = x13*x6
            x82 = 293930.0d0*x81
            x83 = x67*x82 + 10920.0d0*x80 + 7.0d0
            x84 = -6188.0d0*x11 + 4199.0d0*x14 - 364.0d0*x5 + 2730.0d0*x8
            x85 = -340.0d0*x80
            x86 = 3230.0d0*x71
            x87 = 5.0d0*x75 - 1.0d0
            x88 = 1292.0d0*x11 - 969.0d0*x14 + 60.0d0*x5 - 510.0d0*x8
            x89 = x16*x45
            x90 = x2*x32
            x91 = x1*x2
            x92 = x7*x91
            x93 = x9*x91
            x94 = 2469012.0d0*x16
            x95 = 2028117.0d0*x18
            x96 = x2*x95
            x97 = 91.0d0*x75
            x98 = -x97
            x99 = x98 + 7.0d0
            x100 = x2*x37
            x101 = x72*x91
            x102 = -92820.0d0*x101
            x103 = 5460.0d0*x100 + x102
            x104 = x2*x4
            x105 = 91.0d0*x104
            x106 = -x105 + 5460.0d0*x80
            x107 = 155925.0d0*x20
            x108 = -340.0d0*x100
            x109 = 3230.0d0*x90
            x110 = 5.0d0*x104 - 1.0d0
            x111 = zi**4
            x112 = x111*x9
            x113 = 10920.0d0*x100 + x111*x82
            x114 = -182.0d0*x104 + x111*x78 + 7.0d0
            x115 = 11305.0d0*x72
            x116 = x115*x67
            x117 = 1615.0d0*x11
            x118 = x55*x68
            x119 = x117 + x118 - 22610.0d0*x53
            x120 = 47481.0d0*x81
            x121 = -x120*x67
            x122 = x121 - 5950.0d0*x80 - 35.0d0
            x123 = 595.0d0*x7
            x124 = -x123*x67
            x125 = x124 + 350.0d0*x75
            x126 = 525.0d0*x5 - 1785.0d0*x8
            x127 = x121 - 3570.0d0*x80
            x128 = 105.0d0*x5
            x129 = x128 - 7.0d0
            x130 = 323.0d0*x11 - 357.0d0*x8
            x131 = x129 + x130
            x132 = x55*x93
            x133 = x132 - 6783.0d0*x46
            x134 = 105.0d0*x104 - 7.0d0
            x135 = x0*x123
            x136 = x115*x91
            x137 = -x1*x135 + x136
            x138 = x134 + x137
            x139 = -x123*x91
            x140 = -x120*x91 + x139
            x141 = -1785.0d0*x100 + x128 + 35.0d0*x75
            x142 = -6783.0d0*x53
            x143 = 105.0d0*x75
            x144 = -x135*x2
            x145 = x143 + x144
            x146 = 35.0d0*x104 + x136 - 1785.0d0*x80
            x147 = -x111*x120
            x148 = -3570.0d0*x100 + x147
            x149 = -x111*x123
            x150 = x112*x55
            x151 = x111*x115 + x149 + x150
            x152 = x117 - 22610.0d0*x46
            x153 = 350.0d0*x104 - 35.0d0
            x154 = -5950.0d0*x100 + x147
            x155 = 440895.0d0*x16
            x156 = yi**6
            x157 = x156*x6
            x158 = x67*x72
            x159 = x0*x13
            x160 = 62985.0d0*x159
            x161 = x67*x81
            x162 = x19*x9
            x163 = -1105.0d0*x11 - 195.0d0*x5 + 975.0d0*x8 + 5.0d0
            x164 = 1105.0d0*x10
            x165 = -x156*x164
            x166 = x67*x7
            x167 = x165 + 975.0d0*x166 - 195.0d0*x75
            x168 = 4845.0d0*x158
            x169 = x41 - 1275.0d0*x8 - 5.0d0
            x170 = -2550.0d0*x80
            x171 = -33915.0d0*x161 + x170
            x172 = 6081075.0d0*x42
            x173 = 29393.0d0*x16
            x174 = x2*x67
            x175 = 4199.0d0*x13
            x176 = x175*x9
            x177 = x19*x2
            x178 = -19890.0d0*x101
            x179 = -3315.0d0*x158 + x178
            x180 = 1.0d0 - 3315.0d0*x90
            x181 = 195.0d0*x8
            x182 = 20995.0d0*x161
            x183 = x181 + x182 - 6630.0d0*x71
            x184 = 13.0d0*x104
            x185 = x155*x6
            x186 = -x174*x185 - x184 + 1170.0d0*x80
            x187 = 585.0d0*x100 + x160*x174 + 65.0d0*x166 - 26.0d0*x75
            x188 = 176358.0d0*x16
            x189 = x81*x91
            x190 = -221.0d0*x11 - x188*x93 + 125970.0d0*x189 - 39.0d0*x5 + 390.0d0*x92
            x191 = 4845.0d0*x71
            x192 = 4845.0d0*x90
            x193 = 4845.0d0*x101
            x194 = 15.0d0*x104
            x195 = x194 - 3.0d0
            x196 = x193 + x195 - 765.0d0*x80
            x197 = 15.0d0*x75
            x198 = -765.0d0*x100 + x197
            x199 = x111*x72
            x200 = -3315.0d0*x199
            x201 = x1*x111
            x202 = x1*x19
            x203 = x178 - 3315.0d0*x71
            x204 = x111*x81
            x205 = 20995.0d0*x204
            x206 = x205 - 6630.0d0*x90
            x207 = -26.0d0*x104
            x208 = x111*x7
            x209 = x207 + 65.0d0*x208 + 1.0d0
            x210 = x160*x201 + x209 + 585.0d0*x80
            x211 = 1170.0d0*x100 - x185*x201 - 13.0d0*x75
            x212 = 85.0d0*x208
            x213 = 4845.0d0*x199
            x214 = -2550.0d0*x100
            x215 = -33915.0d0*x204 + x214
            x216 = zi**6
            x217 = x216*x6
            x218 = -195.0d0*x104 - x164*x216 + 975.0d0*x208
            x219 = 11305.0d0*x71
            x220 = -595.0d0*x8
            x221 = x220 + 350.0d0*x5
            x222 = 1615.0d0*x10
            x223 = x156*x222
            x224 = 22610.0d0*x159
            x225 = x157*x55
            x226 = -x156*x224 + x223 + x225
            x227 = -1785.0d0*x166 + 525.0d0*x75
            x228 = x191 + 225.0d0*x75
            x229 = 50.0d0*x5 - 85.0d0*x8 - 5.0d0
            x230 = 4845.0d0*x10
            x231 = x174*x230
            x232 = x174*x6
            x233 = 156009.0d0*x16
            x234 = 225.0d0*x104
            x235 = x192 + x234
            x236 = 3230.0d0*x158
            x237 = -255.0d0*x166
            x238 = -6783.0d0*x161
            x239 = x236 + x237 + x238
            x240 = -255.0d0*x8
            x241 = x240 + x86
            x242 = 32300.0d0*x101
            x243 = -67830.0d0*x189 + x242 + 150.0d0*x5 - 2550.0d0*x92 - 15.0d0
            x244 = 67830.0d0*x159
            x245 = -x174*x244 + x214 + 150.0d0*x75
            x246 = 2027025.0d0*xi
            x247 = x23*x246
            x248 = 6460.0d0*x101
            x249 = x195 + x232*x55 + x248 - 1020.0d0*x80
            x250 = -170.0d0*x100 - x174*x224 + 90.0d0*x75
            x251 = -13566.0d0*x189 + 30.0d0*x5 - 51.0d0*x8 - 510.0d0*x92
            x252 = -255.0d0*x208
            x253 = 3230.0d0*x199
            x254 = -6783.0d0*x204
            x255 = x252 + x253 + x254
            x256 = x201*x6
            x257 = x248 - 3.0d0
            x258 = -1020.0d0*x100 + x197 + x256*x55 + x257
            x259 = 90.0d0*x104
            x260 = -x201*x224 + x259 - 170.0d0*x80
            x261 = x201*x230
            x262 = x109 + x240
            x263 = 150.0d0*x104 + x170 - x201*x244
            x264 = x216*x222
            x265 = x217*x55
            x266 = -x216*x224 + x264 + x265
            x267 = 11305.0d0*x90
            x268 = 525.0d0*x104 - 1785.0d0*x208 - 35.0d0
            x269 = yi**8
            x270 = x0*x269
            x271 = x156*x159
            x272 = x19*x6
            x273 = -182.0d0*x5 + 455.0d0*x8
            x274 = x175*x269
            x275 = x10*x156
            x276 = 2730.0d0*x166 + x274 - 6188.0d0*x275 - 364.0d0*x75
            x277 = x220 + 210.0d0*x5
            x278 = 323.0d0*x275
            x279 = x143 - 357.0d0*x166 + x278
            x280 = x2*x6
            x281 = x156*x2
            x282 = x10*x174
            x283 = 1.0d0 - 3315.0d0*x282
            x284 = 195.0d0*x166
            x285 = -6630.0d0*x158 + x284
            x286 = 62985.0d0*x189 - 26.0d0*x5 + 65.0d0*x8 + 585.0d0*x92
            x287 = 221.0d0*x275
            x288 = -x287
            x289 = x0*x188
            x290 = x159*x174
            x291 = 390.0d0*x100 - x281*x289 + x288 + 125970.0d0*x290 - 39.0d0*x75
            x292 = -22610.0d0*x189 + 90.0d0*x5 - 170.0d0*x92
            x293 = 30.0d0*x75
            x294 = -510.0d0*x100 - 51.0d0*x166 - 13566.0d0*x290 + x293
            x295 = x10*x201
            x296 = -6630.0d0*x295
            x298 = x111*x67
            x299 = x0*x298
            x300 = x111*x36
            x301 = x159*x201
            x302 = 195.0d0*x208
            x303 = -6630.0d0*x199 + x302
            x304 = x13*x298
            x305 = 20995.0d0*x304
            x306 = -6630.0d0*x282 + x305
            x307 = x1*x6
            x308 = 30.0d0*x104
            x309 = -51.0d0*x208 - 13566.0d0*x301 + x308 - 510.0d0*x80
            x310 = -3315.0d0*x295
            x311 = x1*x216
            x312 = x159*x216
            x313 = x10*x216
            x314 = 221.0d0*x313
            x315 = 39.0d0*x104
            x316 = x178 - x289*x311 + 125970.0d0*x301 - x314 - x315 + 390.0d0*x80
            x317 = 323.0d0*x313
            x318 = x134 - 357.0d0*x208 + x317
            x319 = zi**8
            x320 = x0*x319
            x321 = -364.0d0*x104 + x175*x319 + 2730.0d0*x208 - 6188.0d0*x313 + 7.0d0
            x322 = 315.0d0*x5
            x323 = 40698.0d0*x158 + x270*x55 - 81396.0d0*x271 + x322
            x324 = x13*x269
            x325 = -6426.0d0*x166 + 11628.0d0*x275 - 6783.0d0*x324 + 1260.0d0*x75
            x326 = 7429.0d0*x16
            x327 = 5.0d0*x5 - 1.0d0
            x328 = 42567525.0d0*xi
            x329 = 6783.0d0*x13
            x330 = x0*x55
            x331 = x281*x330
            x332 = -x281*x329 + x331
            x333 = 35.0d0*x5 - 1785.0d0*x92
            x334 = x144 + x279 - 47481.0d0*x290
            x335 = 135.0d0*x75
            x336 = -6783.0d0*x271
            x337 = 15.0d0*x5
            x338 = x337 - 765.0d0*x92
            x339 = x0*x222
            x340 = -6783.0d0*x304
            x341 = x252 + 3230.0d0*x295 + x340
            x342 = x257 + x299*x55 + x337 - 1020.0d0*x92
            x343 = x237 + 3230.0d0*x282
            x344 = -6783.0d0*x312
            x345 = x311*x330
            x346 = -x311*x329 + x345
            x347 = x137 - 47481.0d0*x301 + x318
            x348 = 60.0d0*x104
            x349 = x13*x319
            x350 = 969.0d0*x349
            x351 = 510.0d0*x208
            x352 = 1292.0d0*x313
            x353 = x23*x328
            x354 = 40698.0d0*x199 - 81396.0d0*x312 + x320*x55 + x322
            x355 = 1260.0d0*x104 - 6426.0d0*x208 + 11628.0d0*x313 - 6783.0d0*x349
            x356 = -273.0d0*x5
            x357 = yi**10
            x358 = x0*x19
            x359 = 46410.0d0*x275
            x360 = x173*x357
            x361 = 1365.0d0*x75
            x362 = 13650.0d0*x166
            x363 = 62985.0d0*x324
            x364 = -x359 - x360 - x361 + x362 + x363
            x365 = 2261.0d0*x324
            x366 = 2142.0d0*x166
            x367 = 420.0d0*x75
            x368 = 3876.0d0*x275
            x369 = 88179.0d0*x16
            x370 = x2*x269
            x371 = x13*x281
            x372 = x0*x94
            x373 = -91.0d0*x5 + 5460.0d0*x92
            x374 = x111*x156
            x375 = x0*x164
            x376 = x0*x175
            x377 = -x155*x299 - 13.0d0*x5 + 1170.0d0*x92
            x378 = x216*x67
            x379 = x13*x311
            x380 = x1*x319
            x381 = 2261.0d0*x349
            x382 = 2142.0d0*x208
            x383 = 420.0d0*x104
            x384 = 3876.0d0*x313
            x385 = zi**10
            x386 = -1365.0d0*x104 - x173*x385 + 13650.0d0*x208 - 46410.0d0*x313 + 62985.0d0*x349
            x387 = -39270.0d0*x166 + 106590.0d0*x275 - 124355.0d0*x324 + x357*x55 + 5775.0d0*x75 - 231.0d0
            x390 = 8.0d0*x75
            x391 = x184 - 1.0d0
            x395 = 728.0d0*x92
            x396 = 16.0d0*x75
            x397 = x26*xi
            x399 = -x383
            x400 = x23*xi
            x403 = x209*x75
            x404 = 10.0d0*x104
            x405 = 17.0d0*x208 - x404 + 1.0d0
            x406 = 11.0d0*x104
            x407 = x406 - 1.0d0
            x408 = x184 - 3.0d0
            x409 = 660.0d0*x92
            x410 = 240.0d0*x208
            x411 = 72.0d0*x75
            x412 = 42525.0d0*x397
            x415 = -20.0d0*x104
            x416 = 80.0d0*x104
            x417 = 80.0d0*x208
            x418 = -x417
            x421 = -x302 + x314 + x315 - 1.0d0
            x423 = x207 + 39.0d0*x208 + 3.0d0
            x424 = 4.0d0*x104
            x426 = 143.0d0*x208
            x427 = -66.0d0*x104 + x426 + 3.0d0
            x430 = 9.0d0*x104 - 1.0d0
            x431 = 16.0d0*x104
            x432 = 16.0d0*x208
            x434 = x406 - 3.0d0
            x435 = 24.0d0*x104
            x437 = 64.0d0*x313
            x438 = -24.0d0*x208 + 16.0d0*x313 + x404 - 1.0d0
            x439 = 2.0d0*x104 - 1.0d0
            x440 = -x431 + x432 + 3.0d0
            x441 = 256.0d0*x349
            x442 = x75*(-40.0d0*x104 - 448.0d0*x313 + x410 + x441 + 1.0d0)
            x443 = x75*(x424 - 1.0d0)
            x444 = x75*(x418 + x435 + x437 - 1.0d0)
            x445 = 336.0d0*x75*(-12.0d0*x104 + x432 + 1.0d0)
            x446 = x16*x385
            x447 = 5775.0d0*x104 - 39270.0d0*x208 + 106590.0d0*x313 - 124355.0d0*x349 + x385*x55 - 231.0d0
            x448 = 20475.0d0*x92 + 21.0d0
            x449 = -7140.0d0*x92 - 63.0d0
            x450 = 293930.0d0*x304 + 10920.0d0*x92
            x451 = -47481.0d0*x304 - 5950.0d0*x92
            T(i, 1) = x21*(-1021020.0d0*x11 + 2078505.0d0*x14 - 1939938.0d0*x17 + x19*xi**12 - 18018.0d0*x5 + 225225.0d0*x8 + 231.0d0)
            T(i, 2) = x23*x25
            T(i, 3) = x25*x26
            T(i, 4) = x21*(-x1*x33 + x1*x34 + x29 - x30*x31 + x35*x36 + x39 + x40)
            T(i, 5) = x43*(9690.0d0*x11 - 14535.0d0*x14 + 7429.0d0*x17 + x41 - 2550.0d0*x8 - 3.0d0)
            T(i, 6) = x21*(-x2*x33 + x2*x34 - x31*x45 + x40 + x44 + 881790.0d0*x46 + x47)
            T(i, 7) = x58*(x50 + x56 + x57)
            T(i, 8) = x61*(x56 + x59 + x60)
            T(i, 9) = x58*(x60 + x63 + x65)
            T(i, 10) = x61*(x57 + x65 + x66)
            T(i, 11) = x21*(235144.0d0*x53 - x67*x73 + x67*x74 - x68*x69 - 176358.0d0*x70 - 92820.0d0*x71 + x79 + x83 + x84)
            T(i, 12) = x43*(-9044.0d0*x53 + 7429.0d0*x70 + x85 + x86 + x87 + x88)
            T(i, 13) = x107*(x103 + x106 + x2*x36*x81 + x30*x96 + 117572.0d0*x46 + 117572.0d0*x53 - 88179.0d0*x70 - 46410.0d0*x71 + x84 - 88179.0d0*x89 - 46410.0d0*x90 + 1365.0d0*x92 - x93*x94 + x99)
            T(i, 14) = x43*(x108 + x109 + x110 - 9044.0d0*x46 + x88 + 7429.0d0*x89)
            T(i, 15) = x21*(-x111*x73 + x111*x74 - x112*x69 + x113 + x114 + 235144.0d0*x46 + x84 - 176358.0d0*x89 - 92820.0d0*x90)
            T(i, 16) = x58*(x116 + x119 + x122 + x125 + x126 + 22610.0d0*x71)
            T(i, 17) = x61*(x116 + x118 + x124 + x127 + x131 - 13566.0d0*x53 + 13566.0d0*x71 + 210.0d0*x75)
            T(i, 18) = x58*(x130 + x133 + x138 + x140 + x141 - 2261.0d0*x53 + 2261.0d0*x71 + 6783.0d0*x90)
            T(i, 19) = x61*(x131 + x132 + x140 + x142 + x145 + x146 - 2261.0d0*x46 + 6783.0d0*x71 + 2261.0d0*x90)
            T(i, 20) = x58*(210.0d0*x104 + x131 + x148 + x151 - 13566.0d0*x46 + 13566.0d0*x90)
            T(i, 21) = x61*(x126 + x151 + x152 + x153 + x154 + 22610.0d0*x90)
            T(i, 22) = x21*(-x155*x157 - x155*x68 + x156*x160 + x156*x162 - 49725.0d0*x158 + 314925.0d0*x161 + x163 + x167 + 62985.0d0*x53 - 49725.0d0*x71 + 8775.0d0*x80)
            T(i, 23) = x172*(x119 - 85.0d0*x166 + x168 + x169 + x171 + 16150.0d0*x71 + 50.0d0*x75)
            T(i, 24) = x21*(-x164*x174 - x173*x68 + x176*x2 + x177*x68 + x179 + x180 + x183 + x186 + x187 + x190 + 8398.0d0*x53)
            T(i, 25) = x172*(969.0d0*x11 + x133 + x142 - 33915.0d0*x189 + x191 + x192 + x196 + x198 + 135.0d0*x5 - 765.0d0*x8 - 85.0d0*x92)
            T(i, 26) = x21*(x1*x176 - x112*x173 + x112*x202 - x164*x201 + x181 + x190 + x200 + x203 + x206 + x210 + x211 + 8398.0d0*x46)
            T(i, 27) = x172*(50.0d0*x104 + x150 + x152 + x169 - x212 + x213 + x215 + 16150.0d0*x90)
            T(i, 28) = x21*(8775.0d0*x100 - x112*x155 - x155*x217 + x160*x216 + x162*x216 + x163 - 49725.0d0*x199 + 314925.0d0*x204 + x218 + 62985.0d0*x46 - 49725.0d0*x90)
            T(i, 29) = x58*(x122 + 22610.0d0*x158 + x219 + x221 + x226 + x227)
            T(i, 30) = x61*(16150.0d0*x158 - 1275.0d0*x166 + x171 + x226 + x228 + x229)
            T(i, 31) = x247*(x231 + x232*x233 + x235 + x239 + x241 + x243 + x245 - 1700.0d0*x80)
            T(i, 32) = x61*(x174*x222 + x239 + x249 + x250 + x251 + 1938.0d0*x71 + 323.0d0*x90)
            T(i, 33) = x58*(x201*x222 + x251 + x255 + x258 + x260 + 323.0d0*x71 + 1938.0d0*x90)
            T(i, 34) = x246*x26*(-1700.0d0*x100 + x228 + x233*x256 + x243 + x255 + x261 + x262 + x263)
            T(i, 35) = x58*(16150.0d0*x199 - 1275.0d0*x208 + x215 + x229 + x235 + x266)
            T(i, 36) = x61*(x154 + 22610.0d0*x199 + x221 + x266 + x267 + x268)
            T(i, 37) = x21*(-x157*x69 - 92820.0d0*x158 - x188*x270 + x269*x272 + 235144.0d0*x271 + x273 + x276 - 30940.0d0*x71 + x83)
            T(i, 38) = x172*(x127 + 13566.0d0*x158 + x219 + x225 - 13566.0d0*x271 + x277 + x279 - 7.0d0)
            T(i, 39) = x21*(-x157*x173 + x157*x177 - x164*x280 + x175*x281 + x182 + x186 + x203 + 8398.0d0*x271 + x283 + x285 + x286 + x291)
            T(i, 40) = x172*(1938.0d0*x158 + x222*x280 + x238 + x241 + x249 + 323.0d0*x282 + x292 + x294)
            T(i, 41) = 93555.0d0*x20*(2340.0d0*x100 - 79560.0d0*x101 - 78.0d0*x104 - 881790.0d0*x16*x232 - 881790.0d0*x16*x299 - x16*x300*x6 + 3380195.0d0*x18*x298*x6 + x183 + 251940.0d0*x189 + x206 + x285 + 251940.0d0*x290 + x296 + 251940.0d0*x301 + x303 + x306 - 78.0d0*x5 - 78.0d0*x75 + 2340.0d0*x80 + 2340.0d0*x92 + 3.0d0)
            T(i, 42) = x172*(1938.0d0*x199 + x222*x307 + x254 + x258 + x262 + x292 + 323.0d0*x295 + x309)
            T(i, 43) = x21*(-x164*x307 - x173*x217 + x175*x311 + x180 + x202*x217 + x205 + x211 + x286 + x303 + x310 + 8398.0d0*x312 + x316)
            T(i, 44) = x172*(x148 + 13566.0d0*x199 + x265 + x267 + x277 - 13566.0d0*x312 + x318)
            T(i, 45) = x21*(x113 - x188*x320 - 92820.0d0*x199 - x217*x69 + x272*x319 + x273 + 235144.0d0*x312 + x321 - 30940.0d0*x90)
            T(i, 46) = x58*(x323 + x325 + x50)
            T(i, 47) = x26*x328*(-510.0d0*x166 + x236 + x270*x326 - 9044.0d0*x271 + 1292.0d0*x275 - 969.0d0*x324 + x327 + 60.0d0*x75 + x85)
            T(i, 48) = x58*(x138 + 2261.0d0*x158 - 2261.0d0*x271 + 6783.0d0*x282 + x332 + x333 + x334)
            T(i, 49) = x61*(-85.0d0*x100 - 765.0d0*x166 + x168 + x196 + x231 + 969.0d0*x275 - 33915.0d0*x290 + x332 + x335 + x336 + x338)
            T(i, 50) = x58*(x111*x339 + 323.0d0*x158 + x260 + 1938.0d0*x282 + x294 + x341 + x342)
            T(i, 51) = x61*(323.0d0*x199 + x250 + 1938.0d0*x295 + x309 + x339*x67 + x340 + x342 + x343)
            T(i, 52) = x58*(135.0d0*x104 + x193 + x198 - 765.0d0*x208 + x213 + x261 - 33915.0d0*x301 + 969.0d0*x313 + x338 + x344 + x346 - 85.0d0*x80 - 3.0d0)
            T(i, 53) = x61*(x145 + 2261.0d0*x199 + 6783.0d0*x295 - 2261.0d0*x312 + x333 + x346 + x347)
            T(i, 54) = x353*(x108 + x253 - 9044.0d0*x312 + x320*x326 + x327 + x348 - x350 - x351 + x352)
            T(i, 55) = x61*(x354 + x355 + x66)
            T(i, 56) = x21*(-232050.0d0*x158 - x270*x31 + 881790.0d0*x271 + x356 + x357*x358 + x364 + x39)
            T(i, 57) = x172*(x323 - x365 - x366 + x367 + x368 + x59)
            T(i, 58) = x107*(1365.0d0*x100 + x102 + x106 - 46410.0d0*x158 - x270*x369 + x270*x96 + 117572.0d0*x271 + x276 - x281*x372 - 46410.0d0*x282 + 881790.0d0*x290 - x369*x370 + 117572.0d0*x371 + x373 + 7.0d0)
            T(i, 59) = x172*(x129 + x139 + x146 + 6783.0d0*x158 + 2261.0d0*x282 + x331 + x334 + x336 - 2261.0d0*x371)
            T(i, 60) = x21*(-x111*x375 + x156*x376 - x173*x374 + x179 + x210 + x284 + x291 + x306 + x310 + x358*x374 + 8398.0d0*x371 + x377)
            T(i, 61) = 2027025.0d0*x42*(x168 + x213 + x233*x299 + x242 + x245 + x263 + x341 + x343 + x41 - 1700.0d0*x92 - 15.0d0)
            T(i, 62) = x21*(-x173*x378 + x187 + x200 + x216*x376 + x283 + x296 + x302 + x305 + x316 + x358*x378 - x375*x67 + x377 + 8398.0d0*x379)
            T(i, 63) = x172*(x139 + x141 + 6783.0d0*x199 + 2261.0d0*x295 + x344 + x345 + x347 - 2261.0d0*x379)
            T(i, 64) = x107*(x1*x320*x95 + x103 + x159*x300 - 46410.0d0*x199 - 46410.0d0*x295 - x311*x372 + 117572.0d0*x312 - x320*x369 + x321 - x369*x380 + x373 + 117572.0d0*x379 + 1365.0d0*x80 + x98)
            T(i, 65) = x172*(x354 - x381 - x382 + x383 + x384 + x63)
            T(i, 66) = x21*(-232050.0d0*x199 - x31*x320 + 881790.0d0*x312 + x356 + x358*x385 + x386 + x47)
            T(i, 67) = x387*x58
            T(i, 68) = x61*(x359 + x360 + x361 - x362 - x363 + 10.0d0*x75*(x365 + x366 - x367 - x368 + 21.0d0) - 21.0d0)
            T(i, 69) = x247*(x390*(-84.0d0*x104 + 459.0d0*x166 - x278 - 3672.0d0*x282 + 3230.0d0*x371 - 189.0d0*x75 + 1134.0d0*x92 + 21.0d0) + 208.0d0*x92*(-459.0d0*x166 + x278 + 189.0d0*x75 - 21.0d0) + (x194 - 1.0d0)*(4914.0d0*x166 + x274 - 7956.0d0*x275 - 1092.0d0*x75 + 63.0d0))
            T(i, 70) = 467775.0d0*x397*(13.0d0*x110*(2002.0d0*x166 - 4004.0d0*x275 + 2431.0d0*x324 - 308.0d0*x75 + 7.0d0) + 88.0d0*x391*x75*(-273.0d0*x166 + x287 + x97 - 7.0d0) + x396*(-28.0d0*x104 + x165 + 1092.0d0*x166 - 3640.0d0*x282 + x29 + 4420.0d0*x371 + x395 + 14.0d0) + 88.0d0*x75*(-28.0d0*x104 + 273.0d0*x166 - 2184.0d0*x282 + x288 + 2210.0d0*x371 + 546.0d0*x92 + x99))
            T(i, 71) = 93555.0d0*x400*(1144.0d0*x110*x92*(x284 + x77 + 35.0d0) + x390*(390.0d0*x166 + 560.0d0*x208 - 7800.0d0*x282 + x29 - 7280.0d0*x295 + 15600.0d0*x304 + x399 + 4368.0d0*x92 + 35.0d0) + 44.0d0*x391*x75*(140.0d0*x104 + 1560.0d0*x282 - x284 + x76 - 1092.0d0*x92 - 35.0d0) + 352.0d0*x92*(70.0d0*x104 - 390.0d0*x166 + x28 + 1300.0d0*x282 - x395 - 35.0d0) + 13.0d0*(x212 - x308 + 1.0d0)*(-1001.0d0*x166 + 715.0d0*x275 + 385.0d0*x75 - 35.0d0))
            T(i, 72) = x412*(x396*(-240.0d0*x104 + 1430.0d0*x166 - 11440.0d0*x282 - 5280.0d0*x295 + 16016.0d0*x304 + x410 - 660.0d0*x75 + 4400.0d0*x92 + 45.0d0) + 198.0d0*x403*(143.0d0*x166 - 110.0d0*x75 + 15.0d0) + 143.0d0*x405*(-495.0d0*x166 + 429.0d0*x275 + x335 - 5.0d0) + 48.0d0*x407*x75*(-858.0d0*x166 + x259 + 2860.0d0*x282 + 495.0d0*x75 - 1320.0d0*x92 - 45.0d0) + 132.0d0*x408*x75*(-143.0d0*x166 + 1144.0d0*x282 + x348 - x409 + 110.0d0*x75 - 15.0d0) + x411*(-180.0d0*x104 + 286.0d0*x166 - 5720.0d0*x282 - 4400.0d0*x295 + 11440.0d0*x304 + x410 - 165.0d0*x75 + 2640.0d0*x92 + 15.0d0))
            T(i, 73) = 42525.0d0*x400*(396.0d0*x403*(x415 - 11.0d0*x75 + 66.0d0*x92 + 5.0d0) + 10296.0d0*x405*x92*(11.0d0*x75 - 5.0d0) + x407*x411*(120.0d0*x104 - 160.0d0*x208 + 880.0d0*x295 + 33.0d0*x75 - 528.0d0*x92 - 10.0d0) + 1056.0d0*x408*x92*(x415 - 33.0d0*x75 + 88.0d0*x92 + 10.0d0) + 32.0d0*x75*(-120.0d0*x104 + 400.0d0*x208 - 2640.0d0*x295 - 320.0d0*x313 + 2464.0d0*x379 + x409 - 22.0d0*x75 + 5.0d0) + 576.0d0*x92*(528.0d0*x295 + x416 + x418 + 66.0d0*x75 - 440.0d0*x92 - 15.0d0) + 143.0d0*(33.0d0*x166 - x293 + 5.0d0)*(45.0d0*x104 + x252 + x317 - 1.0d0))
            T(i, 74) = 14175.0d0*x397*(429.0d0*x318*(21.0d0*x166 - 14.0d0*x75 + 1.0d0) + 6468.0d0*x421*x75*(3.0d0*x75 - 1.0d0) + 6468.0d0*x423*x75*(-x424 - 3.0d0*x75 + 18.0d0*x92 + 1.0d0) + 1176.0d0*x427*x75*(-x424 - 9.0d0*x75 + 24.0d0*x92 + 2.0d0) + 1568.0d0*x430*x75*(144.0d0*x295 + x431 - x432 + 18.0d0*x75 - 120.0d0*x92 - 3.0d0) + 1176.0d0*x434*x75*(-32.0d0*x208 + 240.0d0*x295 + x435 + 9.0d0*x75 - 144.0d0*x92 - 2.0d0) + 896.0d0*x75*(-x197 + 48.0d0*x208 - 504.0d0*x295 - 32.0d0*x313 + 384.0d0*x379 + x415 + 180.0d0*x92 + 2.0d0) + 1568.0d0*x75*(-720.0d0*x295 + 672.0d0*x379 + x417 - x435 - x437 - 6.0d0*x75 + 180.0d0*x92 + 1.0d0))
            T(i, 75) = 33075.0d0*x400*(6864.0d0*x318*x92 + 1848.0d0*x421*x443 + 14784.0d0*x423*x439*x92 + x427*x445 + 448.0d0*x430*x444 + 2688.0d0*x434*x440*x92 + 7168.0d0*x438*x92 + 128.0d0*x442 + 429.0d0*(7.0d0*x75 - 3.0d0)*(-x348 + x350 + x351 - x352 + 1.0d0))
            T(i, 76) = x412*(330.0d0*x321*x75 + 2560.0d0*x438*x75*(7.0d0*x104 - 1.0d0) + 3360.0d0*x439*x75*(33.0d0*x104 + 143.0d0*x313 - x426 - 1.0d0) + 1120.0d0*x440*x75*(-18.0d0*x104 + 33.0d0*x208 + 1.0d0) + 640.0d0*x442 + 1320.0d0*x443*(x105 - 273.0d0*x208 + x314 - 7.0d0) + 2240.0d0*x444*(3.0d0*x104 - 1.0d0) + x445*(-110.0d0*x104 + x426 + 15.0d0) + 256.0d0*x75*(336.0d0*x208 - 512.0d0*x313 - x416 + x441 + 5.0d0) + 143.0d0*x87*(x381 + x382 - x384 + x399 + 21.0d0))
            T(i, 77) = x353*(-2550.0d0*x208 + x234 + 9690.0d0*x313 - 14535.0d0*x349 + 7429.0d0*x446 - 3.0d0)
            T(i, 78) = x447*x61
            T(i, 79) = x21*(-1939938.0d0*x16*x357 + 225225.0d0*x166 + x19*yi**12 - 1021020.0d0*x275 + 2078505.0d0*x324 - 18018.0d0*x75 + 231.0d0)
            T(i, 80) = x172*x387
            T(i, 81) = x21*(x177*x357 - 232050.0d0*x282 - x31*x370 + x364 + 881790.0d0*x371 + x44 + x448)
            T(i, 82) = x172*(40698.0d0*x282 + x325 + x370*x55 - 81396.0d0*x371 + x449 + x64)
            T(i, 83) = x21*(x111*x19*x269 + x114 - x188*x370 + x276 - 92820.0d0*x282 - 30940.0d0*x295 + 235144.0d0*x371 - x374*x69 + x450)
            T(i, 84) = x172*(x149 + x153 + x223 + x227 + 22610.0d0*x282 + 11305.0d0*x295 - 22610.0d0*x371 + x374*x55 + x451)
            T(i, 85) = x21*(-x155*x374 - x155*x378 + x156*x19*x216 + x167 + x218 - 49725.0d0*x282 - 49725.0d0*x295 + 314925.0d0*x304 + 62985.0d0*x371 + 62985.0d0*x379 + 8775.0d0*x92 + 5.0d0)
            T(i, 86) = x172*(x125 + x264 + x268 + 11305.0d0*x282 + 22610.0d0*x295 + x378*x55 - 22610.0d0*x379 + x451)
            T(i, 87) = x21*(-x188*x380 + x19*x319*x67 - 30940.0d0*x282 - 92820.0d0*x295 + x321 - x378*x69 + 235144.0d0*x379 + x450 + x79)
            T(i, 88) = x172*(40698.0d0*x295 + x355 - 81396.0d0*x379 + x380*x55 + x449 + x52)
            T(i, 89) = x21*(x13*x216*x36 + x202*x385 + x29 - 232050.0d0*x295 - x31*x380 + x386 + x448)
            T(i, 90) = x172*x447
            T(i, 91) = x21*(-18018.0d0*x104 + x19*zi**12 + 225225.0d0*x208 - 1021020.0d0*x313 + 2078505.0d0*x349 - 1939938.0d0*x446 + 231.0d0)
        end do
    end subroutine T12
end module T_tensor
