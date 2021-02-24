module T_tensor_damp_amoeba
    implicit none
    private

    public Tn_damp_amoeba
    public T0_damp_amoeba
    public T1_damp_amoeba
    public T2_damp_amoeba
    public T3_damp_amoeba
    public T4_damp_amoeba
    public T5_damp_amoeba
    public T6_damp_amoeba
contains
    pure subroutine Tn_damp_amoeba(n, x, y, z, a, T)
        integer, intent(in) :: n
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(in) :: a(:)
        real(8), intent(inout) :: T(size(x), (n + 1)*(n + 2)*(n + 3)/6)
        select case (n)
        case (0)
            call T0_damp_amoeba(x, y, z, a, T(:, 1:1))
        case (1)
            call T0_damp_amoeba(x, y, z, a, T(:, 1:1))
            call T1_damp_amoeba(x, y, z, a, T(:, 2:4))
        case (2)
            call T0_damp_amoeba(x, y, z, a, T(:, 1:1))
            call T1_damp_amoeba(x, y, z, a, T(:, 2:4))
            call T2_damp_amoeba(x, y, z, a, T(:, 5:10))
        case (3)
            call T0_damp_amoeba(x, y, z, a, T(:, 1:1))
            call T1_damp_amoeba(x, y, z, a, T(:, 2:4))
            call T2_damp_amoeba(x, y, z, a, T(:, 5:10))
            call T3_damp_amoeba(x, y, z, a, T(:, 11:20))
        case (4)
            call T0_damp_amoeba(x, y, z, a, T(:, 1:1))
            call T1_damp_amoeba(x, y, z, a, T(:, 2:4))
            call T2_damp_amoeba(x, y, z, a, T(:, 5:10))
            call T3_damp_amoeba(x, y, z, a, T(:, 11:20))
            call T4_damp_amoeba(x, y, z, a, T(:, 21:35))
        case (5)
            call T0_damp_amoeba(x, y, z, a, T(:, 1:1))
            call T1_damp_amoeba(x, y, z, a, T(:, 2:4))
            call T2_damp_amoeba(x, y, z, a, T(:, 5:10))
            call T3_damp_amoeba(x, y, z, a, T(:, 11:20))
            call T4_damp_amoeba(x, y, z, a, T(:, 21:35))
            call T5_damp_amoeba(x, y, z, a, T(:, 36:56))
        case (6)
            call T0_damp_amoeba(x, y, z, a, T(:, 1:1))
            call T1_damp_amoeba(x, y, z, a, T(:, 2:4))
            call T2_damp_amoeba(x, y, z, a, T(:, 5:10))
            call T3_damp_amoeba(x, y, z, a, T(:, 11:20))
            call T4_damp_amoeba(x, y, z, a, T(:, 21:35))
            call T5_damp_amoeba(x, y, z, a, T(:, 36:56))
            call T6_damp_amoeba(x, y, z, a, T(:, 57:84))
        case default
            call T0_damp_amoeba(x, y, z, a, T(:, 1:1))
            call T1_damp_amoeba(x, y, z, a, T(:, 2:4))
            call T2_damp_amoeba(x, y, z, a, T(:, 5:10))
            call T3_damp_amoeba(x, y, z, a, T(:, 11:20))
            call T4_damp_amoeba(x, y, z, a, T(:, 21:35))
            call T5_damp_amoeba(x, y, z, a, T(:, 36:56))
            call T6_damp_amoeba(x, y, z, a, T(:, 57:84))
        end select
    end subroutine Tn_damp_amoeba
    pure subroutine T0_damp_amoeba(x, y, z, a, T)
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
            T(i, 1) = 0
        end do
    end subroutine T0_damp_amoeba
    pure subroutine T1_damp_amoeba(x, y, z, a, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(in) :: a(:)
        real(8), intent(inout) :: T(size(x), 4)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: ai
        real(8) :: x1
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x1 = (1.0d0 - exp(-ai**3*(xi**2 + yi**2 + zi**2)**1.5d0))*(xi**2 + yi**2 + zi**2)**(-1.5d0)
            T(i, 1) = -x1*xi
            T(i, 2) = -x1*yi
            T(i, 3) = -x1*zi
        end do
    end subroutine T1_damp_amoeba
    pure subroutine T2_damp_amoeba(x, y, z, a, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(in) :: a(:)
        real(8), intent(inout) :: T(size(x), 10)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: ai
        real(8) :: x3
        real(8) :: x8
        real(8) :: x10
        real(8) :: x12
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x3 = xi**2 + yi**2 + zi**2
            x8 = -x3**(-1.5d0)*(1.0d0 - exp(-ai**3*x3**1.5d0))
            x10 = 3.0d0*ai**3*exp(-ai**3*x3**1.5d0)/x3
            x12 = 3.0d0*x3**(-2.5d0)*(1.0d0 - exp(-ai**3*x3**1.5d0))
            T(i, 1) = -x10*xi**2 + x12*xi**2 + x8
            T(i, 2) = 3.0d0*xi*yi*(-ai**3*exp(-ai**3*x3**1.5d0)/x3 + x3**(-2.5d0)*(1.0d0 - exp(-ai**3*x3**1.5d0)))
            T(i, 3) = 3.0d0*xi*zi*(-ai**3*exp(-ai**3*x3**1.5d0)/x3 + x3**(-2.5d0)*(1.0d0 - exp(-ai**3*x3**1.5d0)))
            T(i, 4) = -x10*yi**2 + x12*yi**2 + x8
            T(i, 5) = 3.0d0*yi*zi*(-ai**3*exp(-ai**3*x3**1.5d0)/x3 + x3**(-2.5d0)*(1.0d0 - exp(-ai**3*x3**1.5d0)))
            T(i, 6) = -x10*zi**2 + x12*zi**2 + x8
        end do
    end subroutine T2_damp_amoeba
    pure subroutine T3_damp_amoeba(x, y, z, a, T)
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(in) :: a(:)
        real(8), intent(inout) :: T(size(x), 20)
        integer :: i
        real(8) :: xi, yi, zi
        real(8) :: ai
        real(8) :: x0
        real(8) :: x1
        real(8) :: x3
        real(8) :: x5
        real(8) :: x6
        real(8) :: x8
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x14
        real(8) :: x15
        real(8) :: x16
        real(8) :: x18
        real(8) :: x20
        real(8) :: x21
        real(8) :: x22
        real(8) :: x23
        real(8) :: x24
        real(8) :: x29
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x3 = x0 + x1 + zi**2
            x5 = 5.0d0/x3
            x6 = ai**3
            x8 = exp(-x3**1.5d0*x6)
            x10 = x3**(-2.5d0)*(1.0d0 - x8)
            x11 = sqrt(x3)
            x12 = 1d0/x11
            x14 = 3.0d0*x3*x6
            x15 = x6*x8
            x16 = x15*x3**(-1.5d0)
            x18 = 6.0d0*x15/x3**2
            x20 = -2.0d0*x10 + 2.0d0*x15/x3
            x21 = 3.0d0*xi
            x22 = 5.0d0*x3**(-3.5d0)*(1.0d0 - x8)
            x23 = 3.0d0*ai**6*x8
            x24 = 5.0d0*x15/x3**2
            x29 = x10*(x5*zi**2 - 1.0d0) + x16*(x11 + x12*zi**2 - x14*zi**2) - x18*zi**2
            T(i, 1) = -x21*(-x0*x18 + x10*(x0*x5 - 1.0d0) + x16*(x0*x12 - x0*x14 + x11) + x20)
            T(i, 2) = yi*(3.0d0*x0*x12*x23 - 3.0d0*x0*x22 + 3.0d0*x0*x24 + 3.0d0*x10 - 3.0d0*x15/x3)
            T(i, 3) = zi*(3.0d0*x0*x12*x23 - 3.0d0*x0*x22 + 3.0d0*x0*x24 + 3.0d0*x10 - 3.0d0*x15/x3)
            T(i, 4) = -x21*(-x1*x18 + x10*(x1*x5 - 1.0d0) + x16*(x1*x12 - x1*x14 + x11))
            T(i, 5) = x21*yi*zi*(x12*x23 - x22 + x24)
            T(i, 6) = -x21*x29
            T(i, 7) = -3.0d0*yi*(-x1*x18 + x10*(x1*x5 - 1.0d0) + x16*(x1*x12 - x1*x14 + x11) + x20)
            T(i, 8) = 3.0d0*zi*(x1*x12*x23 - x1*x22 + x1*x24 + x10 - x15/x3)
            T(i, 9) = -3.0d0*x29*yi
            T(i, 10) = -3.0d0*zi*(x20 + x29)
        end do
    end subroutine T3_damp_amoeba
    pure subroutine T4_damp_amoeba(x, y, z, a, T)
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
        real(8) :: x6
        real(8) :: x7
        real(8) :: x8
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x14
        real(8) :: x15
        real(8) :: x16
        real(8) :: x17
        real(8) :: x18
        real(8) :: x22
        real(8) :: x24
        real(8) :: x25
        real(8) :: x26
        real(8) :: x27
        real(8) :: x29
        real(8) :: x32
        real(8) :: x37
        real(8) :: x38
        real(8) :: x39
        real(8) :: x43
        real(8) :: x44
        real(8) :: x45
        real(8) :: x46
        real(8) :: x47
        real(8) :: x55
        real(8) :: x56
        real(8) :: x59
        real(8) :: x62
        real(8) :: x65
        real(8) :: x66
        real(8) :: x70
        real(8) :: x74
        real(8) :: x76
        real(8) :: x88
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = x3**(-2)
            x5 = ai**3
            x6 = x3**1.5d0
            x7 = exp(-x5*x6)
            x8 = x5*x7
            x9 = x4*x8
            x10 = 18.0d0*x9
            x11 = 1d0/x3
            x14 = x3**(-2.5d0)
            x15 = 1.0d0 - x7
            x16 = x14*x15
            x17 = 3.0d0*x16
            x18 = 9.0d0*x5
            x22 = 5.0d0*x15*x3**(-3.5d0)
            x24 = 1d0/x6
            x25 = x24*x8
            x26 = sqrt(x3)
            x27 = 1d0/x26
            x29 = 3.0d0*x3*x5
            x32 = x0*x24
            x37 = x18*x3 - 3.0d0*x27
            x38 = 3.0d0*yi
            x39 = x15*x3**(-4.5d0)
            x43 = x8/x3**3
            x44 = x0*x43
            x45 = ai**6*x7
            x46 = 18.0d0*x45
            x47 = x32*x46
            x55 = 9.0d0*ai**9*x0*x7
            x56 = 35.0d0*x0*x39
            x59 = 5.0d0*x9
            x62 = 3.0d0*x27*x45
            x65 = x38*zi
            x66 = 35.0d0*x2
            x70 = 5.0d0*x1*x11 - 1.0d0
            x74 = x1*x27 - x1*x29 + x26
            x76 = x1*x24
            x88 = x2*x24
            T(i, 1) = 3.0d0*x0*x10 + 3.0d0*x0*x14*x18*x7*(x0*x27 - x0*x29 + x26) - 3.0d0*x0*x18*x4*x7*(5.0d0*x0*x11 - 1.0d0) + 3.0d0*x0*x22*(7.0d0*x0*x11 - 3.0d0) - 3.0d0*x17*(5.0d0*x0*x11 - 1.0d0) - 9.0d0*x25*(x0*x27 - x0*x29 + x26) + 3.0d0*x32*x8*(-9.0d0*ai**6*x0*x6 + x0*x18 + x32 + x37)
            T(i, 2) = x38*xi*(10.0d0*x0*x39 + 3.0d0*x11*x45*(x0*x27 - x0*x29 + x26) + 3.0d0*x14*x8*(x0*x27 - x0*x29 + x26) - 10.0d0*x15*x3**(-3.5d0) + x22*(5.0d0*x0*x11 - 1.0d0) + x25*(6.0d0*x0*x5 - x27 + x32) + 6.0d0*x27*x45 - 24.0d0*x44 - x47 - 3.0d0*x9*(5.0d0*x0*x11 - 1.0d0) + 10.0d0*x9)
            T(i, 3) = 3.0d0*xi*zi*(10.0d0*x0*x39 + 3.0d0*x11*x45*(x0*x27 - x0*x29 + x26) + 3.0d0*x14*x8*(x0*x27 - x0*x29 + x26) - 10.0d0*x15*x3**(-3.5d0) + x22*(5.0d0*x0*x11 - 1.0d0) + x25*(6.0d0*x0*x5 - x27 + x32) + 6.0d0*x27*x45 - 24.0d0*x44 - x47 - 3.0d0*x9*(5.0d0*x0*x11 - 1.0d0) + 10.0d0*x9)
            T(i, 4) = -3.0d0*x0*x22 + 9.0d0*x0*x27*x45 + 3.0d0*x0*x59 - 3.0d0*x1*x22 - 105.0d0*x1*x44 - 3.0d0*x1*x47 - 3.0d0*x1*x55 + 3.0d0*x1*x56 + 3.0d0*x1*x59 + 3.0d0*x1*x62 - 3.0d0*x11*x8 + 3.0d0*x16
            T(i, 5) = -x65*(x22 + 35.0d0*x44 + x47 + x55 - x56 - x59 - x62)
            T(i, 6) = -3.0d0*x0*x22 + 9.0d0*x0*x27*x45 + 3.0d0*x0*x59 - 3.0d0*x11*x8 + 3.0d0*x16 - 3.0d0*x2*x22 - 3.0d0*x2*x47 - 3.0d0*x2*x55 + 3.0d0*x2*x56 + 3.0d0*x2*x59 + 3.0d0*x2*x62 - 3.0d0*x44*x66
            T(i, 7) = x38*xi*(x14*x18*x7*x74 - x18*x4*x7*x70 + x22*(7.0d0*x1*x11 - 3.0d0) + x25*(-9.0d0*ai**6*x1*x6 + x1*x18 + x37 + x76))
            T(i, 8) = 3.0d0*xi*zi*(10.0d0*x1*x39 - 24.0d0*x1*x43 + 3.0d0*x11*x45*x74 + 3.0d0*x14*x74*x8 + x22*x70 + x25*(6.0d0*x1*x5 - x27 + x76) - x46*x76 - 3.0d0*x70*x9)
            T(i, 9) = -x38*xi*(9.0d0*ai**9*x2*x7 + x22 - x39*x66 + x43*x66 + x46*x88 - x59 - x62)
            T(i, 10) = 3.0d0*xi*zi*(x14*x18*x7*(x2*x27 - x2*x29 + x26) - x18*x4*x7*(5.0d0*x11*x2 - 1.0d0) + x22*(7.0d0*x11*x2 - 3.0d0) + x25*(-9.0d0*ai**6*x2*x6 + x18*x2 + x37 + x88))
            T(i, 11) = 3.0d0*x1*x10 + 3.0d0*x1*x14*x18*x7*x74 - 3.0d0*x1*x18*x4*x7*x70 + 3.0d0*x1*x22*(7.0d0*x1*x11 - 3.0d0) - 3.0d0*x17*x70 - 9.0d0*x25*x74 + 3.0d0*x76*x8*(-9.0d0*ai**6*x1*x6 + x1*x18 + x37 + x76)
            T(i, 12) = x65*(10.0d0*x1*x39 - 24.0d0*x1*x43 + 3.0d0*x11*x45*x74 + 3.0d0*x14*x74*x8 - 10.0d0*x15*x3**(-3.5d0) + x22*x70 + x25*(6.0d0*x1*x5 - x27 + x76) + 6.0d0*x27*x45 - x46*x76 - 3.0d0*x70*x9 + 10.0d0*x9)
            T(i, 13) = -27.0d0*ai**9*x1*x2*x7 - 3.0d0*x1*x22 + 3.0d0*x1*x39*x66 - 3.0d0*x1*x43*x66 + 3.0d0*x1*x59 + 3.0d0*x1*x62 - 3.0d0*x11*x8 + 3.0d0*x16 - 3.0d0*x2*x22 - 3.0d0*x2*x46*x76 + 3.0d0*x2*x59 + 3.0d0*x2*x62
            T(i, 14) = x65*(x14*x18*x7*(x2*x27 - x2*x29 + x26) - x18*x4*x7*(5.0d0*x11*x2 - 1.0d0) + x22*(7.0d0*x11*x2 - 3.0d0) + x25*(-9.0d0*ai**6*x2*x6 + x18*x2 + x37 + x88))
            T(i, 15) = 3.0d0*x10*x2 + 3.0d0*x14*x18*x2*x7*(x2*x27 - x2*x29 + x26) - 3.0d0*x17*(5.0d0*x11*x2 - 1.0d0) - 3.0d0*x18*x2*x4*x7*(5.0d0*x11*x2 - 1.0d0) + 3.0d0*x2*x22*(7.0d0*x11*x2 - 3.0d0) - 9.0d0*x25*(x2*x27 - x2*x29 + x26) + 3.0d0*x8*x88*(-9.0d0*ai**6*x2*x6 + x18*x2 + x37 + x88)
        end do
    end subroutine T4_damp_amoeba
    pure subroutine T5_damp_amoeba(x, y, z, a, T)
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
        real(8) :: x8
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x13
        real(8) :: x14
        real(8) :: x15
        real(8) :: x17
        real(8) :: x18
        real(8) :: x19
        real(8) :: x20
        real(8) :: x21
        real(8) :: x24
        real(8) :: x25
        real(8) :: x26
        real(8) :: x27
        real(8) :: x28
        real(8) :: x29
        real(8) :: x30
        real(8) :: x31
        real(8) :: x32
        real(8) :: x34
        real(8) :: x35
        real(8) :: x36
        real(8) :: x37
        real(8) :: x38
        real(8) :: x39
        real(8) :: x42
        real(8) :: x43
        real(8) :: x44
        real(8) :: x45
        real(8) :: x46
        real(8) :: x48
        real(8) :: x49
        real(8) :: x51
        real(8) :: x52
        real(8) :: x53
        real(8) :: x54
        real(8) :: x56
        real(8) :: x57
        real(8) :: x58
        real(8) :: x59
        real(8) :: x60
        real(8) :: x63
        real(8) :: x64
        real(8) :: x65
        real(8) :: x66
        real(8) :: x67
        real(8) :: x70
        real(8) :: x72
        real(8) :: x74
        real(8) :: x75
        real(8) :: x77
        real(8) :: x80
        real(8) :: x84
        real(8) :: x87
        real(8) :: x88
        real(8) :: x89
        real(8) :: x90
        real(8) :: x92
        real(8) :: x93
        real(8) :: x94
        real(8) :: x96
        real(8) :: x98
        real(8) :: x100
        real(8) :: x102
        real(8) :: x104
        real(8) :: x105
        real(8) :: x106
        real(8) :: x108
        real(8) :: x109
        real(8) :: x110
        real(8) :: x112
        real(8) :: x113
        real(8) :: x115
        real(8) :: x116
        real(8) :: x117
        real(8) :: x118
        real(8) :: x119
        real(8) :: x123
        real(8) :: x127
        real(8) :: x129
        real(8) :: x130
        real(8) :: x137
        real(8) :: x145
        real(8) :: x147
        real(8) :: x148
        real(8) :: x152
        real(8) :: x155
        real(8) :: x157
        real(8) :: x158
        real(8) :: x159
        real(8) :: x163
        real(8) :: x164
        real(8) :: x166
        real(8) :: x171
        real(8) :: x178
        real(8) :: x180
        real(8) :: x188
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0/x3
            x5 = x0*x4
            x6 = 5.0d0*x5 - 1.0d0
            x8 = x3**(-2)
            x9 = ai**3
            x10 = x3**1.5d0
            x11 = exp(-x10*x9)
            x12 = x11*x9
            x13 = x12*x8
            x14 = 36.0d0*x13
            x15 = 7.0d0*x5 - 3.0d0
            x17 = 1.0d0 - x11
            x18 = x17*x3**(-3.5d0)
            x19 = 20.0d0*x18
            x20 = x12/x3**3
            x21 = x0*x20
            x24 = 21.0d0*x8
            x25 = 15.0d0*x18
            x26 = x3**(-2.5d0)
            x27 = sqrt(x3)
            x28 = 1d0/x27
            x29 = 3.0d0*x9
            x30 = x29*x3
            x31 = x0*x28 - x0*x30 + x27
            x32 = x26*x31
            x34 = 18.0d0*x9
            x35 = 1d0/x10
            x36 = x0*x35
            x37 = 9.0d0*x9
            x38 = ai**6
            x39 = 9.0d0*x38
            x42 = -3.0d0*x28 + x3*x37
            x43 = -x0*x10*x39 + x0*x37 + x36 + x42
            x44 = x12*x35
            x45 = 4.0d0*x44
            x46 = x0*x26
            x48 = 12.0d0*x12
            x49 = x4*x9
            x51 = 9.0d0*ai**9*x3**2
            x52 = 18.0d0*x38
            x53 = x27*x52
            x54 = x28 - x30
            x56 = x11*x29*x35
            x57 = 3.0d0*xi
            x58 = x12/x3**4
            x59 = 90.0d0*x58
            x60 = x11*x38
            x63 = x17*x3**(-4.5d0)
            x64 = x0*x63
            x65 = x17*x3**(-5.5d0)
            x66 = 70.0d0*x65
            x67 = x11*x37
            x70 = 35.0d0*x63
            x72 = 6.0d0*x9
            x74 = x0*x72 - x28 + x36
            x75 = x11*x39
            x77 = x12*x3**(-3.5d0)
            x80 = x60*x8
            x84 = x11*x29
            x87 = 70.0d0*x63
            x88 = ai**9*x11
            x89 = 18.0d0*x88
            x90 = x1*x35
            x92 = x6*x70
            x93 = 36.0d0*x60
            x94 = 70.0d0*x20
            x96 = 140.0d0*x0*x65
            x98 = 204.0d0*x0*x58
            x100 = 126.0d0*x46*x60
            x102 = 54.0d0*x5*x88
            x104 = 27.0d0*x20*x6
            x105 = x6*x75
            x106 = x1*x26
            x108 = x1*x4
            x109 = 6.0d0*x60
            x110 = x109*x74
            x112 = 9.0d0*x88
            x113 = x112*x31
            x115 = 15.0d0*x31*x77
            x116 = 15.0d0*x31*x80
            x117 = 5.0d0*x18
            x118 = 3.0d0*x60
            x119 = x11*x52
            x123 = x11*x26
            x127 = x2*x35
            x129 = x2*x4
            x130 = x2*x28
            x137 = 35.0d0*x20
            x145 = 3.0d0*yi
            x147 = 3.0d0*zi
            x148 = ai**12*x11*x2*x27
            x152 = 7.0d0*x108 - 3.0d0
            x155 = yi**4
            x157 = x123*x72
            x158 = 5.0d0*x108 - 1.0d0
            x159 = x1*x28 - x1*x30 + x27
            x163 = -x1*x10*x39 + x1*x37 + x42 + x90
            x164 = x106*x163
            x166 = 9.0d0*xi
            x171 = x1*x72 - x28 + x90
            x178 = x1*x20
            x180 = x1*x63
            x188 = 54.0d0*x108*x2*x88
            T(i, 1) = -x57*(x11*x32*x34*x6 - 36.0d0*x12*x32 + x14*x6 - x15*x19 - 60.0d0*x15*x21 + x25*(x24*xi**4 - 14.0d0*x5 + 1.0d0) - x43*x45 + x43*x46*x48 + x56*(x0*x10*x52 - x0*x34 + x26*xi**4 - 2.0d0*x36 + x49*xi**4 - x51*xi**4 + x53*xi**4 + x54))
            T(i, 2) = -yi*(3.0d0*x0*x15*x70 + 135.0d0*x0*x31*x77 + 81.0d0*x0*x31*x80 - 45.0d0*x15*x21 - 108.0d0*x21*x6 + 216.0d0*x21 - 3.0d0*x25*x6 - 3.0d0*x31*x4*x75 - 3.0d0*x32*x67 - 81.0d0*x36*x6*x60 + 162.0d0*x36*x60 + 3.0d0*x36*x84*(x0*x27*x39 - x35 + x46 - x72) + 3.0d0*x43*x46*x84 + 9.0d0*x43*x5*x60 + 3.0d0*x46*x67*x74 - 3.0d0*x56*x74 - 3.0d0*x59*xi**4 + 3.0d0*x6*x67*x8 - 90.0d0*x64 + 3.0d0*x66*xi**4)
            T(i, 3) = -zi*(3.0d0*x0*x15*x70 + 135.0d0*x0*x31*x77 + 81.0d0*x0*x31*x80 - 45.0d0*x15*x21 - 108.0d0*x21*x6 + 216.0d0*x21 - 3.0d0*x25*x6 - 3.0d0*x31*x4*x75 - 3.0d0*x32*x67 - 81.0d0*x36*x6*x60 + 162.0d0*x36*x60 + 3.0d0*x36*x84*(x0*x27*x39 - x35 + x46 - x72) + 3.0d0*x43*x46*x84 + 9.0d0*x43*x5*x60 + 3.0d0*x46*x67*x74 - 3.0d0*x56*x74 - 3.0d0*x59*xi**4 + 3.0d0*x6*x67*x8 - 90.0d0*x64 + 3.0d0*x66*xi**4)
            T(i, 4) = -x57*(-x1*x100 - x1*x102 - x1*x104 + x1*x113*x28 + x1*x115 + x1*x116 - x1*x87 + x1*x89 + x1*x92 + x1*x94 + x1*x96 - x1*x98 - x105*x90 + x106*x11*x72*x74 + x108*x110 - x109*x28 - x117*x6 - x118*x31*x4 + x119*x36 - 10.0d0*x13 + 10.0d0*x18 + 24.0d0*x21 - x32*x84 - x44*(-3.0d0*x1*x46 + x74 + x90) + x6*x8*x84 - 10.0d0*x64 + x90*x93)
            T(i, 5) = -x57*yi*zi*(-x100 - x102 - x104 - x105*x35 + x110*x4 + x113*x28 + x115 + x116 + x123*x72*x74 + x20*(3.0d0*x5 - 1.0d0) + x35*x93 - x87 + x89 + x92 + x94 + x96 - x98)
            T(i, 6) = -x57*(-x100*x2 - x102*x2 - x104*x2 - x105*x127 - x109*x28 + x110*x129 + x113*x130 + x115*x2 + x116*x2 - x117*x6 - x118*x31*x4 + x119*x36 + x123*x2*x72*x74 + x127*x93 - 10.0d0*x13 + 10.0d0*x18 - x2*x87 + x2*x89 + x2*x92 + x2*x94 + x2*x96 - x2*x98 + 24.0d0*x21 - x32*x84 - x44*(x127 - 3.0d0*x2*x46 + x74) + x6*x8*x84 - 10.0d0*x64)
            T(i, 7) = x145*(27.0d0*ai**12*x0*x1*x11*x27 + 315.0d0*x0*x1*x58 - 315.0d0*x0*x1*x65 - 27.0d0*x0*x88 + x1*x102 - x1*x112 - x1*x137 + 159.0d0*x1*x46*x60 + x1*x70 - x119*x90 + 15.0d0*x13 - 105.0d0*x21 - x25 + x28*x75 - 54.0d0*x36*x60 + 105.0d0*x64)
            T(i, 8) = x147*(27.0d0*ai**12*x0*x1*x11*x27 + 315.0d0*x0*x1*x58 - 315.0d0*x0*x1*x65 - x0*x112 - x0*x137 + x0*x70 + x1*x102 - x1*x112 - x1*x137 + 159.0d0*x1*x46*x60 + x1*x70 - x117 - x119*x36 - x119*x90 + 5.0d0*x13 + 3.0d0*x28*x60)
            T(i, 9) = x145*(-x0*x112 - x0*x137 + 27.0d0*x0*x148 + 315.0d0*x0*x2*x58 - 315.0d0*x0*x2*x65 + x0*x70 + x102*x2 - x112*x2 - x117 - x119*x127 - x119*x36 + 5.0d0*x13 - x137*x2 + 159.0d0*x2*x46*x60 + x2*x70 + 3.0d0*x28*x60)
            T(i, 10) = x147*(27.0d0*x0*x148 + 315.0d0*x0*x2*x58 - 315.0d0*x0*x2*x65 - 27.0d0*x0*x88 + x102*x2 - x112*x2 - x119*x127 + 15.0d0*x13 - x137*x2 + 159.0d0*x2*x46*x60 + x2*x70 - 105.0d0*x21 - x25 + x28*x75 - 54.0d0*x36*x60 + 105.0d0*x64)
            T(i, 11) = -x166*(-20.0d0*x1*x152*x20 + x117*(-14.0d0*x108 + x155*x24 + 1.0d0) + 4.0d0*x12*x164 + x157*x158*x159 + x44*(x1*x10*x52 - x1*x34 + x155*x26 + x155*x49 - x155*x51 + x155*x53 + x54 - 2.0d0*x90))
            T(i, 12) = -x57*yi*zi*(-x1*x59 + x1*x66 + x118*x163*x4 - 15.0d0*x152*x20 + x152*x70 - 36.0d0*x158*x20 - 27.0d0*x158*x35*x60 + 45.0d0*x159*x77 + 27.0d0*x159*x80 + x163*x26*x84 + x171*x26*x67 + x56*(x1*x27*x39 + x106 - x35 - x72))
            T(i, 13) = x57*(204.0d0*x1*x2*x58 - 140.0d0*x1*x2*x65 + 126.0d0*x106*x2*x60 - x109*x129*x171 - x112*x130*x159 + x117*x158 + x118*x159*x4 - x119*x90 + x127*x158*x75 - x157*x171*x2 + 27.0d0*x158*x2*x20 - x158*x2*x70 - x158*x8*x84 - 15.0d0*x159*x2*x77 - 15.0d0*x159*x2*x80 + x159*x26*x84 - 24.0d0*x178 + 10.0d0*x180 + x188 + x44*(-3.0d0*x106*x2 + x127 + x171))
            T(i, 14) = x166*yi*zi*(-x11*x35*x52 - x112 + x129*x89 - x137 + 9.0d0*x148 + 53.0d0*x2*x26*x60 + 105.0d0*x2*x58 - 105.0d0*x2*x65 + x70)
            T(i, 15) = -x166*(x117*(-14.0d0*x129 + x24*zi**4 + 1.0d0) + 4.0d0*x12*x2*x26*(-x10*x2*x39 + x127 + x2*x37 + x42) + x157*(5.0d0*x129 - 1.0d0)*(x130 - x2*x30 + x27) - 20.0d0*x2*x20*(7.0d0*x129 - 3.0d0) + x44*(x10*x2*x52 - 2.0d0*x127 - x2*x34 + x26*zi**4 + x49*zi**4 - x51*zi**4 + x53*zi**4 + x54))
            T(i, 16) = -x145*(-60.0d0*x1*x152*x20 - 36.0d0*x12*x159*x26 + x123*x158*x159*x34 + x14*x158 - x152*x19 - x163*x45 + x164*x48 + x25*(-14.0d0*x108 + x155*x24 + 1.0d0) + x56*(x1*x10*x52 - x1*x34 + x155*x26 + x155*x49 - x155*x51 + x155*x53 + x54 - 2.0d0*x90))
            T(i, 17) = -x147*(-15.0d0*x1*x152*x20 + x1*x152*x70 - 36.0d0*x1*x158*x20 + 45.0d0*x1*x159*x77 + 27.0d0*x1*x159*x80 + x106*x171*x67 + x108*x118*x163 - x155*x59 + x155*x66 - x158*x25 - 27.0d0*x158*x60*x90 + x158*x67*x8 - x159*x26*x67 - x159*x4*x75 + x164*x84 - x171*x56 + 72.0d0*x178 - 30.0d0*x180 + 54.0d0*x60*x90 + x84*x90*(x1*x27*x39 + x106 - x35 - x72))
            T(i, 18) = -x145*(-204.0d0*x1*x2*x58 + 140.0d0*x1*x2*x65 - 126.0d0*x106*x2*x60 + x109*x129*x171 - x109*x28 + x112*x130*x159 - x117*x158 - x118*x159*x4 + x119*x90 - x127*x158*x75 + x127*x93 - 10.0d0*x13 + x157*x171*x2 - 27.0d0*x158*x2*x20 + x158*x2*x70 + x158*x8*x84 + 15.0d0*x159*x2*x77 + 15.0d0*x159*x2*x80 - x159*x26*x84 + 24.0d0*x178 + 10.0d0*x18 - 10.0d0*x180 - x188 - x2*x87 + x2*x89 + x2*x94 - x44*(-3.0d0*x106*x2 + x127 + x171))
            T(i, 19) = x147*(27.0d0*x1*x148 + 315.0d0*x1*x2*x58 - 315.0d0*x1*x2*x65 - 27.0d0*x1*x88 + 159.0d0*x106*x2*x60 - x112*x2 - x119*x127 + 15.0d0*x13 - x137*x2 - 105.0d0*x178 + 105.0d0*x180 + x188 + x2*x70 - x25 + x28*x75 - 54.0d0*x60*x90)
            T(i, 20) = -9.0d0*yi*(x117*(-14.0d0*x129 + x24*zi**4 + 1.0d0) + 4.0d0*x12*x2*x26*(-x10*x2*x39 + x127 + x2*x37 + x42) + x157*(5.0d0*x129 - 1.0d0)*(x130 - x2*x30 + x27) - 20.0d0*x2*x20*(7.0d0*x129 - 3.0d0) + x44*(x10*x2*x52 - 2.0d0*x127 - x2*x34 + x26*zi**4 + x49*zi**4 - x51*zi**4 + x53*zi**4 + x54))
            T(i, 21) = -x147*(-36.0d0*x12*x26*(x130 - x2*x30 + x27) + x123*x34*(5.0d0*x129 - 1.0d0)*(x130 - x2*x30 + x27) + x14*(5.0d0*x129 - 1.0d0) - x19*(7.0d0*x129 - 3.0d0) - 60.0d0*x2*x20*(7.0d0*x129 - 3.0d0) + x2*x26*x48*(-x10*x2*x39 + x127 + x2*x37 + x42) + x25*(-14.0d0*x129 + x24*zi**4 + 1.0d0) - x45*(-x10*x2*x39 + x127 + x2*x37 + x42) + x56*(x10*x2*x52 - 2.0d0*x127 - x2*x34 + x26*zi**4 + x49*zi**4 - x51*zi**4 + x53*zi**4 + x54))
        end do
    end subroutine T5_damp_amoeba
    pure subroutine T6_damp_amoeba(x, y, z, a, T)
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
        real(8) :: x9
        real(8) :: x10
        real(8) :: x11
        real(8) :: x12
        real(8) :: x13
        real(8) :: x15
        real(8) :: x17
        real(8) :: x18
        real(8) :: x20
        real(8) :: x21
        real(8) :: x23
        real(8) :: x24
        real(8) :: x25
        real(8) :: x26
        real(8) :: x27
        real(8) :: x29
        real(8) :: x30
        real(8) :: x31
        real(8) :: x32
        real(8) :: x33
        real(8) :: x35
        real(8) :: x36
        real(8) :: x37
        real(8) :: x39
        real(8) :: x40
        real(8) :: x41
        real(8) :: x42
        real(8) :: x43
        real(8) :: x44
        real(8) :: x48
        real(8) :: x49
        real(8) :: x50
        real(8) :: x51
        real(8) :: x54
        real(8) :: x55
        real(8) :: x56
        real(8) :: x58
        real(8) :: x59
        real(8) :: x60
        real(8) :: x61
        real(8) :: x62
        real(8) :: x65
        real(8) :: x67
        real(8) :: x68
        real(8) :: x71
        real(8) :: x73
        real(8) :: x74
        real(8) :: x75
        real(8) :: x76
        real(8) :: x78
        real(8) :: x79
        real(8) :: x81
        real(8) :: x82
        real(8) :: x83
        real(8) :: x85
        real(8) :: x86
        real(8) :: x88
        real(8) :: x89
        real(8) :: x90
        real(8) :: x93
        real(8) :: x94
        real(8) :: x97
        real(8) :: x98
        real(8) :: x99
        real(8) :: x100
        real(8) :: x101
        real(8) :: x103
        real(8) :: x104
        real(8) :: x106
        real(8) :: x107
        real(8) :: x109
        real(8) :: x110
        real(8) :: x112
        real(8) :: x113
        real(8) :: x116
        real(8) :: x117
        real(8) :: x119
        real(8) :: x121
        real(8) :: x122
        real(8) :: x124
        real(8) :: x128
        real(8) :: x130
        real(8) :: x132
        real(8) :: x134
        real(8) :: x135
        real(8) :: x137
        real(8) :: x140
        real(8) :: x141
        real(8) :: x142
        real(8) :: x145
        real(8) :: x147
        real(8) :: x148
        real(8) :: x149
        real(8) :: x150
        real(8) :: x152
        real(8) :: x153
        real(8) :: x154
        real(8) :: x155
        real(8) :: x156
        real(8) :: x157
        real(8) :: x158
        real(8) :: x159
        real(8) :: x162
        real(8) :: x167
        real(8) :: x168
        real(8) :: x169
        real(8) :: x170
        real(8) :: x171
        real(8) :: x172
        real(8) :: x173
        real(8) :: x185
        real(8) :: x187
        real(8) :: x188
        real(8) :: x204
        real(8) :: x210
        real(8) :: x211
        real(8) :: x213
        real(8) :: x215
        real(8) :: x218
        real(8) :: x219
        real(8) :: x230
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
        real(8) :: x248
        real(8) :: x249
        real(8) :: x250
        real(8) :: x251
        real(8) :: x253
        real(8) :: x254
        real(8) :: x255
        real(8) :: x258
        real(8) :: x259
        real(8) :: x260
        real(8) :: x262
        real(8) :: x263
        real(8) :: x265
        real(8) :: x273
        real(8) :: x276
        real(8) :: x279
        real(8) :: x284
        real(8) :: x285
        real(8) :: x289
        real(8) :: x294
        real(8) :: x296
        real(8) :: x300
        real(8) :: x302
        real(8) :: x314
        real(8) :: x319
        real(8) :: x320
        real(8) :: x321
        real(8) :: x322
        real(8) :: x324
        real(8) :: x326
        real(8) :: x330
        real(8) :: x338
        real(8) :: x340
        real(8) :: x341
        real(8) :: x344
        real(8) :: x351
        real(8) :: x352
        real(8) :: x353
        real(8) :: x355
        real(8) :: x358
        real(8) :: x361
        real(8) :: x364
        real(8) :: x374
        real(8) :: x379
        real(8) :: x382
        real(8) :: x384
        real(8) :: x385
        real(8) :: x387
        real(8) :: x388
        real(8) :: x389
        real(8) :: x391
        real(8) :: x395
        real(8) :: x398
        real(8) :: x399
        real(8) :: x400
        real(8) :: x402
        real(8) :: x404
        real(8) :: x408
        real(8) :: x409
        real(8) :: x415
        real(8) :: x421
        real(8) :: x424
        real(8) :: x425
        real(8) :: x428
        real(8) :: x429
        real(8) :: x435
        real(8) :: x436
        real(8) :: x441
        real(8) :: x450
        real(8) :: x456
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x1 = yi**2
            x2 = zi**2
            x3 = x0 + x1 + x2
            x4 = 1d0/x3
            x5 = x0*x4
            x6 = 7.0d0*x5 - 3.0d0
            x7 = ai**3
            x8 = x0*x7
            x9 = x3**(-3)
            x10 = x3**1.5d0
            x11 = exp(-x10*x7)
            x12 = x11*x9
            x13 = x12*x8
            x15 = xi**4
            x17 = x3**(-2)
            x18 = x15*x17
            x20 = x3**(-3.5d0)
            x21 = 1.0d0 - x11
            x23 = 25.0d0*x20*x21
            x24 = x3**(-4.5d0)
            x25 = x21*x24
            x26 = 35.0d0*x25
            x27 = x0*x26
            x29 = x3**(-2.5d0)
            x30 = 5.0d0*x5 - 1.0d0
            x31 = sqrt(x3)
            x32 = 1d0/x31
            x33 = x0*x32
            x35 = -3.0d0*x3*x8 + x31 + x33
            x36 = x11*x7
            x37 = x35*x36
            x39 = x11*x8
            x40 = 50.0d0*x20
            x41 = 1d0/x10
            x42 = x0*x41
            x43 = ai**6
            x44 = x10*x43
            x48 = 9.0d0*x3*x7 - 3.0d0*x32
            x49 = -9.0d0*x0*x44 + x42 + x48 + 9.0d0*x8
            x50 = x0*x29
            x51 = x36*x50
            x54 = 5.0d0*x41
            x55 = x4*x7
            x56 = ai**9
            x58 = 9.0d0*x3**2*x56
            x59 = x15*x43
            x60 = 18.0d0*x31
            x61 = 3.0d0*x3*x7
            x62 = x32 - x61
            x65 = 5.0d0*x20
            x67 = 90.0d0*x3**2*x56
            x68 = x3*x56
            x71 = 10.0d0*x7
            x73 = ai**12
            x74 = 27.0d0*x73
            x75 = x3**2.5d0*x74
            x76 = x0*x43
            x78 = 45.0d0*x7
            x79 = 45.0d0*x44 - x54 - x78
            x81 = xi*yi
            x82 = x3**(-4)
            x83 = x39*x82
            x85 = x3**(-5)
            x86 = x36*x85
            x88 = x21*x3**(-5.5d0)
            x89 = x0*x88
            x90 = x36*x9
            x93 = x11*x43
            x94 = x41*x93
            x97 = x50*x93
            x98 = x6*x97
            x99 = 3.0d0*x5
            x100 = x99 - 1.0d0
            x101 = 420.0d0*x89
            x103 = 105.0d0*x25
            x104 = 6.0d0*x7
            x106 = x0*x104 - x32 + x42
            x107 = x29*x36
            x109 = 180.0d0*x20
            x110 = x35*x93
            x112 = x24*x35
            x113 = x112*x39
            x116 = -x104 + 9.0d0*x31*x76 - x41 + x50
            x117 = x36*x41
            x119 = 90.0d0*x20
            x121 = x4*x93
            x122 = x39*x49
            x124 = x11*x76
            x128 = 3.0d0*x117
            x130 = 9.0d0*x121
            x132 = xi*zi
            x134 = x1*x41
            x135 = 3.0d0*x50
            x137 = x21*x3**(-6.5d0)
            x140 = x1*x7
            x141 = x11*x140
            x142 = 1500.0d0*x85
            x145 = 540.0d0*x20
            x147 = x1*x11
            x148 = x20*x78
            x149 = x148*x35
            x150 = x1*x17
            x152 = x1*x32
            x153 = x11*x56
            x154 = 27.0d0*x153
            x155 = x154*x35
            x156 = 18.0d0*x106
            x157 = x1*x29
            x158 = x157*x36
            x159 = x1*x4
            x162 = 3.0d0*x157
            x167 = x134*x93
            x168 = 27.0d0*x30
            x169 = x12*x140
            x170 = 81.0d0*x30
            x171 = x153*x159
            x172 = 162.0d0*x0
            x173 = 315.0d0*x89
            x185 = 9.0d0*x153
            x187 = x33*x49
            x188 = 15.0d0*x20
            x204 = 70.0d0*x88
            x210 = x36*x82
            x211 = x0*x11
            x213 = 36.0d0*x30
            x215 = x104*x11
            x218 = 5.0d0*x0
            x219 = x11*x82
            x230 = 54.0d0*x153
            x236 = 9.0d0*yi*zi
            x237 = x2*x41
            x239 = x2*x7
            x240 = x11*x239
            x241 = x11*x2
            x242 = x17*x2
            x243 = x2*x32
            x244 = x2*x29
            x245 = x244*x36
            x246 = x2*x4
            x248 = x2*x43
            x249 = x237*x93
            x250 = x12*x239
            x251 = x153*x246
            x253 = x185*x2
            x254 = 15.0d0*x242
            x255 = 18.0d0*x239
            x258 = x1*x88
            x259 = 1.0d0 - x99
            x260 = x0*x137
            x262 = x140*x219
            x263 = x157*x93
            x265 = x147*x73
            x273 = 18.0d0*x153
            x276 = x106*x93
            x279 = x265*x33
            x284 = x39*x85
            x285 = x1*x284
            x289 = 70.0d0*x25
            x294 = 9.0d0*x81
            x296 = 3.0d0*x159
            x300 = 54.0d0*x31
            x302 = 108.0d0*x171
            x314 = x2*x88
            x319 = x241*x73
            x320 = x300*x319
            x321 = 108.0d0*x251
            x322 = x244*x93
            x324 = x219*x239
            x326 = x2*x284
            x330 = x319*x33
            x338 = 105.0d0*x314
            x340 = 9.0d0*x132
            x341 = yi**4
            x344 = x11*x31
            x351 = 1155.0d0*x0*x137
            x352 = 580.0d0*x20
            x353 = x17*x341
            x355 = 195.0d0*x0*x153
            x358 = 27.0d0*ai**15*x3
            x361 = -x0*x185 - 35.0d0*x13 + 5.0d0*x17*x36 - x21*x65 + x27 + 3.0d0*x32*x93 - 18.0d0*x42*x93
            x364 = x352*x76
            x374 = x1*x2
            x379 = x150*x2
            x382 = x1*x26
            x384 = zi**4
            x385 = 1155.0d0*x137*x384
            x387 = x11*x384
            x388 = x17*x384
            x389 = 195.0d0*x153*x388
            x391 = 45.0d0*x11*x384*x73
            x395 = -14.0d0*x159 + 21.0d0*x353 + 1.0d0
            x398 = -x1*x61 + x152 + x31
            x399 = x36*x398
            x400 = 7.0d0*x159 - 3.0d0
            x402 = 5.0d0*x159 - 1.0d0
            x404 = -9.0d0*x1*x44 + x134 + 9.0d0*x140 + x48
            x408 = 18.0d0*x1*x44 - 2.0d0*x134 - 18.0d0*x140 + x29*x341 + x341*x43*x60 + x341*x55 - x341*x58 + x62
            x409 = x107*x408
            x415 = x1*x31*x43
            x421 = x11*x398
            x424 = x1*x104 + x134 - x32
            x425 = -x104 + x157 - x41 + 9.0d0*x415
            x428 = x398*x93
            x429 = x17*x428
            x435 = 3.0d0*x121
            x436 = x140*x241
            x441 = 27.0d0*x402
            x450 = x428*x9
            x456 = -x162*x2 + x237 + x424
            T(i, 1) = 900.0d0*x13*x6 - 675.0d0*x13*(21.0d0*x18 - 14.0d0*x5 + 1.0d0) - 9.0d0*x23*(21.0d0*x18 - 14.0d0*x5 + 1.0d0) + 9.0d0*x27*(33.0d0*x18 - 30.0d0*x5 + 5.0d0) - 270.0d0*x29*x30*x37 + 90.0d0*x30*x36*x49*x50 + 9.0d0*x35*x39*x40*x6 - 9.0d0*x36*x42*(-x0*x67 - x15*x65 - 90.0d0*x15*x68 + x15*x75 - 5.0d0*x18*x7 + 180.0d0*x31*x76 + 15.0d0*x32*x59 + x5*x71 + 10.0d0*x50 + x79) + 135.0d0*x36*x50*(18.0d0*x0*x44 + x15*x29 + x15*x55 - x15*x58 - 2.0d0*x42 + x59*x60 + x62 - 18.0d0*x8) - 9.0d0*x36*x54*(18.0d0*x0*x44 + x15*x29 + x15*x55 - x15*x58 - 2.0d0*x42 + x59*x60 + x62 - 18.0d0*x8) - 180.0d0*x49*x51
            T(i, 2) = x81*(3.0d0*x100*x101 + 3.0d0*x103*(21.0d0*x18 - 14.0d0*x5 + 1.0d0) + 54.0d0*x106*x107*x30 - 108.0d0*x106*x107 - 36.0d0*x107*x49 - 3.0d0*x109*x37 + 162.0d0*x110*x17*x30 - 324.0d0*x110*x17 + 540.0d0*x113 - 36.0d0*x116*x117 + 108.0d0*x116*x51 + 3.0d0*x119*x30*x37 - 3.0d0*x12*x78*(21.0d0*x18 - 14.0d0*x5 + 1.0d0) - 36.0d0*x121*x49 + 180.0d0*x122*x20 + 108.0d0*x124*x17*x49 + 3.0d0*x128*(x104 + x15*x65 + 36.0d0*x15*x68 + 2.0d0*x18*x7 - 54.0d0*x31*x76 - 18.0d0*x32*x59 + x41 - 6.0d0*x50) + 3.0d0*x130*(18.0d0*x0*x44 + x15*x29 + x15*x55 - x15*x58 - 2.0d0*x42 + x59*x60 + x62 - 18.0d0*x8) - 2520.0d0*x15*x86 - 420.0d0*x25*x6 + 27.0d0*x29*x36*(18.0d0*x0*x44 + x15*x29 + x15*x55 - x15*x58 - 2.0d0*x42 + x59*x60 + x62 - 18.0d0*x8) + 432.0d0*x30*x90 + 324.0d0*x30*x94 - 1080.0d0*x6*x83 + 180.0d0*x6*x90 + 1080.0d0*x83 - 840.0d0*x89 - 540.0d0*x98)
            T(i, 3) = x132*(3.0d0*x100*x101 + 3.0d0*x103*(21.0d0*x18 - 14.0d0*x5 + 1.0d0) + 54.0d0*x106*x107*x30 - 108.0d0*x106*x107 - 36.0d0*x107*x49 - 3.0d0*x109*x37 + 162.0d0*x110*x17*x30 - 324.0d0*x110*x17 + 540.0d0*x113 - 36.0d0*x116*x117 + 108.0d0*x116*x51 + 3.0d0*x119*x30*x37 - 3.0d0*x12*x78*(21.0d0*x18 - 14.0d0*x5 + 1.0d0) - 36.0d0*x121*x49 + 180.0d0*x122*x20 + 108.0d0*x124*x17*x49 + 3.0d0*x128*(x104 + x15*x65 + 36.0d0*x15*x68 + 2.0d0*x18*x7 - 54.0d0*x31*x76 - 18.0d0*x32*x59 + x41 - 6.0d0*x50) + 3.0d0*x130*(18.0d0*x0*x44 + x15*x29 + x15*x55 - x15*x58 - 2.0d0*x42 + x59*x60 + x62 - 18.0d0*x8) - 2520.0d0*x15*x86 - 420.0d0*x25*x6 + 27.0d0*x29*x36*(18.0d0*x0*x44 + x15*x29 + x15*x55 - x15*x58 - 2.0d0*x42 + x59*x60 + x62 - 18.0d0*x8) + 432.0d0*x30*x90 + 324.0d0*x30*x94 - 1080.0d0*x6*x83 + 180.0d0*x6*x90 + 1080.0d0*x83 - 840.0d0*x89 - 540.0d0*x98)
            T(i, 4) = -3.0d0*x0*x170*x171 + 90.0d0*x0*x25 - 3.0d0*x1*x101 - 3.0d0*x1*x103*x30 + 3.0d0*x1*x106*x119*x39 - 3.0d0*x1*x11*x145*x59 + 945.0d0*x1*x113 + 729.0d0*x1*x12*x35*x76 + 3.0d0*x1*x122*x188 + 3780.0d0*x1*x137*x15 + 243.0d0*x1*x153*x35*x42 + 3.0d0*x1*x173*x6 + 3.0d0*x1*x185*x187 - 648.0d0*x1*x30*x83 - 567.0d0*x1*x30*x97 - 585.0d0*x1*x6*x83 + 1836.0d0*x1*x83 + 1134.0d0*x1*x97 - 135.0d0*x1*x98 + 162.0d0*x106*x124*x150 + 54.0d0*x11*x116*x140*x50 - 135.0d0*x110*x150 + 27.0d0*x110*x4 + 54.0d0*x116*x124*x159 + 45.0d0*x124*x150*x49 - 81.0d0*x124*x17*x35 + 3.0d0*x128*(-x1*x135 + x106 + x134) + 3.0d0*x13*x213 + 45.0d0*x13*x6 - 216.0d0*x13 - 3.0d0*x135*x36*x49 - 3.0d0*x141*x142*x15 - 3.0d0*x147*x149 - 3.0d0*x149*x211 - 3.0d0*x15*x204 + 270.0d0*x15*x210 - 3.0d0*x152*x155 - 3.0d0*x156*x158 - 3.0d0*x156*x159*x93 + 3.0d0*x167*x168 + 3.0d0*x168*x42*x93 + 3.0d0*x169*x170 - 27.0d0*x17*x30*x36 + 3.0d0*x171*x172 + 45.0d0*x20*x21*x30 - 3.0d0*x27*x6 + 27.0d0*x29*x37 - 9.0d0*x36*x42*(-x0*x1*x65 + 9.0d0*x1*x33*x43 + x116 + x162) - 162.0d0*x42*x93 - 3.0d0*x49*x93*x99 - 27.0d0*x51*(-x1*x135 + x106 + x134)
            T(i, 5) = x236*(3.0d0*x100*x219*x8 - x100*x90 - 6.0d0*x106*x121 + 30.0d0*x106*x20*x39 - x106*x215*x29 - x109*x11*x59 - 15.0d0*x110*x17 + 105.0d0*x113 + x116*x215*x50 + 6.0d0*x116*x5*x93 + 81.0d0*x12*x35*x76 + x122*x65 + x124*x156*x17 + 420.0d0*x137*x15 - 500.0d0*x15*x86 + 3.0d0*x153*x187 - x154*x30*x5 + x155*x42 + x17*x218*x49*x93 - x17*x39*(3.0d0*x17 - x218*x9 + 9.0d0*x76) - x185*x32*x35 - x188*x37 + x230*x5 - x26*x30 - 72.0d0*x30*x83 + 27.0d0*x30*x90 + 9.0d0*x30*x94 - 63.0d0*x30*x97 - 65.0d0*x6*x83 + 105.0d0*x6*x89 + 204.0d0*x83 - 140.0d0*x89 + 126.0d0*x97 - 15.0d0*x98)
            T(i, 6) = -3.0d0*x0*x170*x251 + 90.0d0*x0*x25 - 3.0d0*x101*x2 - 3.0d0*x103*x2*x30 + 3.0d0*x106*x119*x2*x39 + 162.0d0*x106*x124*x242 + 3.0d0*x11*x116*x255*x50 - 3.0d0*x11*x145*x2*x59 - 135.0d0*x110*x242 + 27.0d0*x110*x4 + 945.0d0*x113*x2 + 54.0d0*x116*x124*x246 + 729.0d0*x12*x2*x35*x76 + 3.0d0*x122*x188*x2 - 81.0d0*x124*x17*x35 + 3.0d0*x124*x254*x49 + 3.0d0*x128*(x106 - x135*x2 + x237) + 3.0d0*x13*x213 + 45.0d0*x13*x6 - 216.0d0*x13 - 3.0d0*x135*x36*x49 + 3780.0d0*x137*x15*x2 - 3.0d0*x142*x15*x240 - 3.0d0*x149*x211 - 3.0d0*x149*x241 - 3.0d0*x15*x204 + 270.0d0*x15*x210 + 243.0d0*x153*x2*x35*x42 - 3.0d0*x155*x243 - 3.0d0*x156*x245 - 3.0d0*x156*x246*x93 + 3.0d0*x168*x249 + 3.0d0*x168*x42*x93 - 27.0d0*x17*x30*x36 + 3.0d0*x170*x250 + 3.0d0*x172*x251 + 3.0d0*x173*x2*x6 + 3.0d0*x187*x253 - 648.0d0*x2*x30*x83 - 567.0d0*x2*x30*x97 - 585.0d0*x2*x6*x83 + 1836.0d0*x2*x83 + 1134.0d0*x2*x97 - 135.0d0*x2*x98 + 45.0d0*x20*x21*x30 - 3.0d0*x27*x6 + 27.0d0*x29*x37 - 9.0d0*x36*x42*(-x0*x2*x65 + x116 + 3.0d0*x244 + 9.0d0*x248*x33) - 162.0d0*x42*x93 - 3.0d0*x49*x93*x99 - 27.0d0*x51*(x106 - x135*x2 + x237)
            T(i, 7) = -x294*(-35.0d0*x1*x110*x9 - 630.0d0*x1*x260 + 3.0d0*x106*x107 + 3.0d0*x106*x121 - x106*x141*x188 - x106*x152*x185 + 3.0d0*x107*(-x1*x135 + x106 + x134) + 15.0d0*x110*x17 - 35.0d0*x112*x141 + 3.0d0*x121*(-x1*x135 + x106 + x134) - x134*x273*x35 + 444.0d0*x147*x20*x76 + x150*x153*x172 - 15.0d0*x150*x276 + x159*x185*x30 - 36.0d0*x171 + x185*x32*x35 + x188*x37 + x213*x263 - x230*x5 - 105.0d0*x258*x30 + 210.0d0*x258 + x26*x30 + 89.0d0*x262*x30 - 210.0d0*x262 - 106.0d0*x263 - 9.0d0*x265*x35 - x265*x60 + x273 + 54.0d0*x279 + 774.0d0*x285 - x289 - 27.0d0*x30*x90 - 9.0d0*x30*x94 - 204.0d0*x83 + 140.0d0*x89 - x90*(x150*x218 - x159 + x259) + 70.0d0*x90 + 36.0d0*x94 - 126.0d0*x97)
            T(i, 8) = 3.0d0*x132*(-486.0d0*x0*x150*x153 + x1*x100*x104*x219 + 105.0d0*x1*x110*x9 + 1890.0d0*x1*x260 + 6.0d0*x100*x157*x93 - 3.0d0*x106*x107 - 3.0d0*x106*x121 + x106*x147*x148 + x106*x152*x154 - 3.0d0*x107*(-x1*x135 + x106 + x134) - 15.0d0*x110*x17 + 105.0d0*x112*x141 - 3.0d0*x121*(-x1*x135 + x106 + x134) + x134*x230*x35 - 1332.0d0*x147*x20*x76 + x147*x35*x74 + 45.0d0*x150*x276 - x168*x171 - x185*x32*x35 - x188*x37 + x230*x5 + 315.0d0*x258*x30 - 630.0d0*x258 - x26*x30 - 267.0d0*x262*x30 + 630.0d0*x262 - 108.0d0*x263*x30 + 318.0d0*x263 + x265*x300 - x273 - 162.0d0*x279 - 2322.0d0*x285 + x289 + 27.0d0*x30*x90 + 9.0d0*x30*x94 + x302 + 204.0d0*x83 - 140.0d0*x89 + x90*(15.0d0*x0*x150 + x259 - x296) - 70.0d0*x90 - 36.0d0*x94 + 126.0d0*x97)
            T(i, 9) = 3.0d0*x81*(-486.0d0*x0*x153*x242 + x100*x2*x215*x82 + 6.0d0*x100*x322 - 3.0d0*x106*x107 - 3.0d0*x106*x121 + x106*x148*x241 + x106*x154*x243 - 3.0d0*x107*(x106 - x135*x2 + x237) - 15.0d0*x110*x17 + 105.0d0*x110*x2*x9 + 105.0d0*x112*x240 - 3.0d0*x121*(x106 - x135*x2 + x237) - x168*x251 - x185*x32*x35 - x188*x37 + 1890.0d0*x2*x260 - 1332.0d0*x20*x241*x76 + x230*x237*x35 + x230*x5 + x241*x35*x74 + 45.0d0*x242*x276 - x26*x30 - x273 + x289 + 315.0d0*x30*x314 - 108.0d0*x30*x322 - 267.0d0*x30*x324 + 27.0d0*x30*x90 + 9.0d0*x30*x94 - 630.0d0*x314 + x320 + x321 + 318.0d0*x322 + 630.0d0*x324 - 2322.0d0*x326 - 162.0d0*x330 + 204.0d0*x83 - 140.0d0*x89 + x90*(15.0d0*x0*x242 - 3.0d0*x246 + x259) - 70.0d0*x90 - 36.0d0*x94 + 126.0d0*x97)
            T(i, 10) = -x340*(3.0d0*x106*x107 + 3.0d0*x106*x121 - x106*x185*x243 - x106*x188*x240 + 3.0d0*x107*(x106 - x135*x2 + x237) + 15.0d0*x110*x17 - 35.0d0*x110*x2*x9 - 35.0d0*x112*x240 + 3.0d0*x121*(x106 - x135*x2 + x237) + x153*x172*x242 + x185*x246*x30 + x185*x32*x35 + x188*x37 - 630.0d0*x2*x260 + 444.0d0*x20*x241*x76 + x213*x322 - x230*x5 - x237*x273*x35 - 36.0d0*x251 - x254*x276 + x26*x30 + x273 - x289 + 89.0d0*x30*x324 - x30*x338 - 27.0d0*x30*x90 - 9.0d0*x30*x94 + 210.0d0*x314 - 9.0d0*x319*x35 - x319*x60 - 106.0d0*x322 - 210.0d0*x324 + 774.0d0*x326 + 54.0d0*x330 - 204.0d0*x83 + 140.0d0*x89 - x90*(x218*x242 - x246 + x259) + 70.0d0*x90 + 36.0d0*x94 - 126.0d0*x97)
            T(i, 11) = 9.0d0*x0*x265*x300 + 9.0d0*x0*x302 - 9.0d0*x1*x273 + 9.0d0*x1*x289 + 5670.0d0*x1*x83 - 5670.0d0*x1*x89 + 2862.0d0*x1*x97 - 405.0d0*x11*x33*x341*x73 - 9.0d0*x124*x341*x352 - 324.0d0*x167 - 630.0d0*x169 + 945.0d0*x210*x341 - 9.0d0*x211*x341*x358 + 9.0d0*x273*x341*x4 - 10395.0d0*x284*x341 + 477.0d0*x29*x341*x93 + 81.0d0*x341*x344*x73 + 9.0d0*x341*x351 - 945.0d0*x341*x88 - 9.0d0*x353*x355 + 9.0d0*x361
            T(i, 12) = -x236*(x0*x147*x358 - x0*x344*x74 - x1*x351 + x147*x364 + x150*x355 - x159*x273 + x173 + x185 - x230*x5 + 105.0d0*x258 - x26 - 105.0d0*x262 - 53.0d0*x263 - 9.0d0*x265*x31 + 45.0d0*x279 + 1155.0d0*x285 - 315.0d0*x83 + 35.0d0*x90 + 18.0d0*x94 - 159.0d0*x97)
            T(i, 13) = -243.0d0*ai**15*x211*x3*x374 + 3.0d0*x0*x147*x31*x74 - 1755.0d0*x0*x153*x379 + 3.0d0*x0*x159*x230 + 3.0d0*x0*x230*x246 + 3.0d0*x0*x241*x31*x74 - 3.0d0*x1*x173 - 3.0d0*x1*x185 + 945.0d0*x1*x83 + 477.0d0*x1*x97 - 405.0d0*x11*x33*x374*x73 - 5220.0d0*x124*x20*x374 + 3.0d0*x159*x2*x230 - 54.0d0*x167 - 105.0d0*x169 - 3.0d0*x173*x2 - 945.0d0*x2*x258 + 3.0d0*x2*x26 + 945.0d0*x2*x262 + 477.0d0*x2*x263 + 945.0d0*x2*x83 + 477.0d0*x2*x97 - 54.0d0*x249 - 105.0d0*x250 - 3.0d0*x253 + 10395.0d0*x260*x374 - 10395.0d0*x284*x374 + 3.0d0*x344*x374*x74 + 3.0d0*x361 + 3.0d0*x382
            T(i, 14) = -x236*(x0*x241*x358 - x0*x344*x74 + x173 + x185 - x2*x351 - x230*x5 + x241*x364 + x242*x355 - x246*x273 - x26 - 9.0d0*x31*x319 - 53.0d0*x322 - 105.0d0*x324 + 1155.0d0*x326 + 45.0d0*x330 + x338 - 315.0d0*x83 + 35.0d0*x90 + 18.0d0*x94 - 159.0d0*x97)
            T(i, 15) = 9.0d0*x0*x320 + 9.0d0*x0*x321 - 9.0d0*x0*x358*x387 + 9.0d0*x0*x385 - 9.0d0*x0*x389 - 9.0d0*x2*x273 + 9.0d0*x2*x289 + 5670.0d0*x2*x83 - 5670.0d0*x2*x89 + 2862.0d0*x2*x97 + 945.0d0*x210*x384 - 324.0d0*x249 - 630.0d0*x250 + 9.0d0*x273*x384*x4 - 10395.0d0*x284*x384 + 477.0d0*x29*x384*x93 - 9.0d0*x33*x391 + 81.0d0*x344*x384*x73 + 9.0d0*x361 - 9.0d0*x364*x387 - 945.0d0*x384*x88
            T(i, 16) = x294*(x11*x29*x402*x404*x71 - x117*(-x1*x67 + 10.0d0*x157 + x159*x71 + 15.0d0*x32*x341*x43 - x341*x65 - 90.0d0*x341*x68 + x341*x75 - 5.0d0*x353*x7 + 180.0d0*x415 + x79) + x26*(-30.0d0*x159 + 33.0d0*x353 + 5.0d0) - 75.0d0*x395*x90 + x399*x40*x400 + 15.0d0*x409)
            T(i, 17) = x340*(x117*(x104 - 6.0d0*x157 - 18.0d0*x32*x341*x43 + x341*x65 + 36.0d0*x341*x68 + 2.0d0*x353*x7 + x41 - 54.0d0*x415) + 60.0d0*x140*x24*x421 + 20.0d0*x141*x20*x404 + 12.0d0*x150*x404*x93 + 12.0d0*x158*x425 + 30.0d0*x20*x399*x402 + x215*x29*x402*x424 + 140.0d0*x258*(x296 - 1.0d0) + x26*x395 - 120.0d0*x262*x400 - 60.0d0*x263*x400 - 280.0d0*x341*x86 - 15.0d0*x395*x90 + 18.0d0*x402*x429 + x408*x435 + 3.0d0*x409)
            T(i, 18) = 3.0d0*x81*(-x1*x204 - 3.0d0*x107*x404 - 9.0d0*x107*x456 + x119*x240*x424 - x128*(9.0d0*x152*x248 + 3.0d0*x244 - x374*x65 + x425) + 1260.0d0*x137*x374 - x142*x436 - x145*x374*x93 - x148*x421 + 81.0d0*x153*x237*x398 - 27.0d0*x17*x428 + x185*x243*x404 + x188*x240*x404 + 243.0d0*x2*x450 + 315.0d0*x239*x24*x421 + 54.0d0*x242*x424*x93 + 18.0d0*x245*x425 + 18.0d0*x246*x425*x93 - 81.0d0*x251*x402 + x254*x404*x93 - x26*x400 + 90.0d0*x262 + 315.0d0*x314*x400 - 45.0d0*x322*x400 - 189.0d0*x322*x402 - 195.0d0*x324*x400 - 216.0d0*x324*x402 + 15.0d0*x400*x90 + 36.0d0*x402*x90 - x404*x435 + x441*x94)
            T(i, 19) = -x340*(3.0d0*x107*x424 + 3.0d0*x107*x456 - 630.0d0*x137*x374 + 54.0d0*x152*x319 + 162.0d0*x153*x379 - x159*x230 - x185*x243*x424 + x185*x246*x402 + x185*x32*x398 - x188*x240*x424 + x188*x399 - 35.0d0*x2*x450 + 444.0d0*x20*x374*x93 - x237*x273*x398 - 35.0d0*x239*x24*x421 - x254*x424*x93 + 140.0d0*x258 + x26*x402 - 204.0d0*x262 - 126.0d0*x263 - 9.0d0*x319*x398 + 36.0d0*x322*x402 + 89.0d0*x324*x402 - x338*x402 - 9.0d0*x402*x94 + x424*x435 + 15.0d0*x429 + x435*x456 + 774.0d0*x436*x85 - x441*x90 - x90*(-x246 - x296 + 5.0d0*x379 + 1.0d0))
            T(i, 20) = -x294*(x185 - x26 + 630.0d0*x314 + x32*x391 - x320 - x321 - 318.0d0*x322 - 630.0d0*x324 + x352*x384*x93 + x358*x387 + 1155.0d0*x384*x86 - x385 + x389 + 35.0d0*x90 + 18.0d0*x94)
            T(i, 21) = x340*(15.0d0*x107*(18.0d0*x2*x44 - 2.0d0*x237 - x255 + x29*x384 + x384*x43*x60 + x384*x55 - x384*x58 + x62) + x11*x29*x71*(5.0d0*x246 - 1.0d0)*(-9.0d0*x2*x44 + x237 + 9.0d0*x239 + x48) - x117*(-x2*x67 + 10.0d0*x244 + x246*x71 + 180.0d0*x248*x31 + 15.0d0*x32*x384*x43 - x384*x65 - 90.0d0*x384*x68 + x384*x75 - 5.0d0*x388*x7 + x79) + x26*(-30.0d0*x246 + 33.0d0*x388 + 5.0d0) + x36*x40*(7.0d0*x246 - 3.0d0)*(-x2*x61 + x243 + x31) - 75.0d0*x90*(-14.0d0*x246 + 21.0d0*x388 + 1.0d0))
            T(i, 22) = -270.0d0*x107*x398*x402 - 9.0d0*x134*x36*(-x1*x67 + 10.0d0*x157 + x159*x71 + 15.0d0*x32*x341*x43 - x341*x65 - 90.0d0*x341*x68 + x341*x75 - 5.0d0*x353*x7 + 180.0d0*x415 + x79) + 9.0d0*x140*x40*x400*x421 + 90.0d0*x157*x36*x402*x404 - 180.0d0*x158*x404 + 135.0d0*x158*x408 - 675.0d0*x169*x395 + 900.0d0*x169*x400 - 9.0d0*x23*x395 - 9.0d0*x36*x408*x54 + 9.0d0*x382*(-30.0d0*x159 + 33.0d0*x353 + 5.0d0)
            T(i, 23) = 3.0d0*yi*zi*(x103*x395 + 18.0d0*x107*x402*x424 - 12.0d0*x107*x404 - 36.0d0*x107*x424 - x109*x399 - 12.0d0*x117*x425 + 3.0d0*x117*(x104 - 6.0d0*x157 - 18.0d0*x32*x341*x43 + x341*x65 + 36.0d0*x341*x68 + 2.0d0*x353*x7 + x41 - 54.0d0*x415) + x119*x399*x402 - x12*x395*x78 - 12.0d0*x121*x404 + x130*x408 + 180.0d0*x140*x24*x421 + 60.0d0*x141*x20*x404 + 36.0d0*x150*x404*x93 + 36.0d0*x158*x425 - 140.0d0*x25*x400 + 420.0d0*x258*(x296 - 1.0d0) - 280.0d0*x258 - 360.0d0*x262*x400 + 360.0d0*x262 - 180.0d0*x263*x400 - 840.0d0*x341*x86 + 60.0d0*x400*x90 + 54.0d0*x402*x429 + 144.0d0*x402*x90 + 108.0d0*x402*x94 + 9.0d0*x409 - 108.0d0*x429)
            T(i, 24) = 90.0d0*x1*x25 - 3.0d0*x103*x2*x402 + 27.0d0*x107*x398 + 3.0d0*x11*x157*x255*x425 + 3.0d0*x119*x424*x436 + 3.0d0*x128*x456 + 3.0d0*x130*x398 + 243.0d0*x134*x153*x2*x398 - 9.0d0*x134*x36*(9.0d0*x152*x248 + 3.0d0*x244 - x374*x65 + x425) + 3780.0d0*x137*x2*x341 - 3.0d0*x142*x240*x341 - 3.0d0*x145*x2*x341*x93 - 3.0d0*x147*x148*x398 - 3.0d0*x148*x241*x398 - 81.0d0*x150*x428 + 3.0d0*x152*x253*x404 - 3.0d0*x154*x243*x398 - 27.0d0*x158*x456 + 54.0d0*x159*x2*x425*x93 - 3.0d0*x162*x36*x404 + 3.0d0*x167*x441 - 162.0d0*x167 + 45.0d0*x169*x400 + 108.0d0*x169*x402 - 216.0d0*x169 - 27.0d0*x17*x36*x402 - 243.0d0*x171*x2*x402 + 486.0d0*x171*x2 + 3.0d0*x188*x404*x436 + 945.0d0*x2*x258*x400 - 1260.0d0*x2*x258 - 585.0d0*x2*x262*x400 - 648.0d0*x2*x262*x402 + 1836.0d0*x2*x262 - 135.0d0*x2*x263*x400 - 567.0d0*x2*x263*x402 + 1134.0d0*x2*x263 + 45.0d0*x20*x21*x402 - 3.0d0*x204*x341 + 270.0d0*x210*x341 + 945.0d0*x24*x398*x436 - 135.0d0*x242*x428 - 54.0d0*x245*x424 - 54.0d0*x246*x424*x93 + 3.0d0*x249*x441 + 243.0d0*x250*x402 - 3.0d0*x296*x404*x93 + 729.0d0*x374*x450 + 45.0d0*x379*x404*x93 + 162.0d0*x379*x424*x93 - 3.0d0*x382*x400
            T(i, 25) = -x236*(3.0d0*x107*x424 + 3.0d0*x107*x456 - 630.0d0*x137*x374 + 54.0d0*x152*x319 + 162.0d0*x153*x379 - x159*x230 - x185*x243*x424 + x185*x246*x402 + x185*x32*x398 - x188*x240*x424 + x188*x399 - 35.0d0*x2*x450 + 444.0d0*x20*x374*x93 - x237*x273*x398 - 35.0d0*x239*x24*x421 - 36.0d0*x251 - x254*x424*x93 + 140.0d0*x258 + x26*x402 - 204.0d0*x262 - 126.0d0*x263 + x273 - x289 + 210.0d0*x314 - 9.0d0*x319*x398 - x319*x60 + 36.0d0*x322*x402 - 106.0d0*x322 + 89.0d0*x324*x402 - 210.0d0*x324 - x338*x402 - 9.0d0*x402*x94 + x424*x435 + 15.0d0*x429 + x435*x456 + 774.0d0*x436*x85 - x441*x90 - x90*(-x246 - x296 + 5.0d0*x379 + 1.0d0) + 70.0d0*x90 + 36.0d0*x94)
            T(i, 26) = -9.0d0*x1*x185 - 9.0d0*x1*x352*x384*x93 + 9.0d0*x1*x385 - 9.0d0*x1*x389 + 9.0d0*x11*x300*x374*x73 - 10395.0d0*x140*x387*x85 - 9.0d0*x147*x358*x384 - 9.0d0*x152*x391 - 162.0d0*x167 - 315.0d0*x169 + 45.0d0*x17*x36 - 5670.0d0*x2*x258 + 5670.0d0*x2*x262 + 2862.0d0*x2*x263 - 9.0d0*x2*x273 + 9.0d0*x2*x289 + 9.0d0*x2*x302 - 9.0d0*x21*x65 + 945.0d0*x210*x384 - 324.0d0*x249 - 630.0d0*x250 + 9.0d0*x273*x384*x4 + 477.0d0*x29*x384*x93 + 27.0d0*x32*x93 + 81.0d0*x344*x384*x73 + 9.0d0*x382 - 945.0d0*x384*x88
            T(i, 27) = x236*(15.0d0*x107*(18.0d0*x2*x44 - 2.0d0*x237 - x255 + x29*x384 + x384*x43*x60 + x384*x55 - x384*x58 + x62) + x11*x29*x71*(5.0d0*x246 - 1.0d0)*(-9.0d0*x2*x44 + x237 + 9.0d0*x239 + x48) - x117*(-x2*x67 + 10.0d0*x244 + x246*x71 + 180.0d0*x248*x31 + 15.0d0*x32*x384*x43 - x384*x65 - 90.0d0*x384*x68 + x384*x75 - 5.0d0*x388*x7 + x79) + x26*(-30.0d0*x246 + 33.0d0*x388 + 5.0d0) + x36*x40*(7.0d0*x246 - 3.0d0)*(-x2*x61 + x243 + x31) - 75.0d0*x90*(-14.0d0*x246 + 21.0d0*x388 + 1.0d0))
            T(i, 28) = -270.0d0*x107*(5.0d0*x246 - 1.0d0)*(-x2*x61 + x243 + x31) + 9.0d0*x2*x26*(-30.0d0*x246 + 33.0d0*x388 + 5.0d0) - 9.0d0*x23*(-14.0d0*x246 + 21.0d0*x388 + 1.0d0) - 9.0d0*x237*x36*(-x2*x67 + 10.0d0*x244 + x246*x71 + 180.0d0*x248*x31 + 15.0d0*x32*x384*x43 - x384*x65 - 90.0d0*x384*x68 + x384*x75 - 5.0d0*x388*x7 + x79) + 9.0d0*x240*x40*(7.0d0*x246 - 3.0d0)*(-x2*x61 + x243 + x31) + 90.0d0*x244*x36*(5.0d0*x246 - 1.0d0)*(-9.0d0*x2*x44 + x237 + 9.0d0*x239 + x48) - 180.0d0*x245*(-9.0d0*x2*x44 + x237 + 9.0d0*x239 + x48) + 135.0d0*x245*(18.0d0*x2*x44 - 2.0d0*x237 - x255 + x29*x384 + x384*x43*x60 + x384*x55 - x384*x58 + x62) + 900.0d0*x250*(7.0d0*x246 - 3.0d0) - 675.0d0*x250*(-14.0d0*x246 + 21.0d0*x388 + 1.0d0) - 9.0d0*x36*x54*(18.0d0*x2*x44 - 2.0d0*x237 - x255 + x29*x384 + x384*x43*x60 + x384*x55 - x384*x58 + x62)
        end do
    end subroutine T6_damp_amoeba
end module T_tensor_damp_amoeba
