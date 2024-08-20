module T_tensor_damp_thole
    use precision
    implicit none
    private
    public Tn_damp_thole
    public T0_damp_thole
    public T1_damp_thole
    public T2_damp_thole
    public T3_damp_thole
    public T4_damp_thole
    public T5_damp_thole
    public T6_damp_thole
contains
    pure subroutine Tn_damp_thole(n, x, y, z, a, T)
        integer(ip), intent(in) :: n
        real(rp), intent(in) :: x(:), y(size(x)), z(size(x))
        real(rp), intent(in) :: a(:)
        real(rp), intent(inout) :: T(size(x), (n + 1) * (n + 2) * (n + 3) / 6)
        select case (n)
        case (0)
            call T0_damp_thole(x, y, z, a, T(:, 1:1))
        case (1)
            call T0_damp_thole(x, y, z, a, T(:, 1:1))
            call T1_damp_thole(x, y, z, a, T(:, 2:4))
        case (2)
            call T0_damp_thole(x, y, z, a, T(:, 1:1))
            call T1_damp_thole(x, y, z, a, T(:, 2:4))
            call T2_damp_thole(x, y, z, a, T(:, 5:10))
        case (3)
            call T0_damp_thole(x, y, z, a, T(:, 1:1))
            call T1_damp_thole(x, y, z, a, T(:, 2:4))
            call T2_damp_thole(x, y, z, a, T(:, 5:10))
            call T3_damp_thole(x, y, z, a, T(:, 11:20))
        case (4)
            call T0_damp_thole(x, y, z, a, T(:, 1:1))
            call T1_damp_thole(x, y, z, a, T(:, 2:4))
            call T2_damp_thole(x, y, z, a, T(:, 5:10))
            call T3_damp_thole(x, y, z, a, T(:, 11:20))
            call T4_damp_thole(x, y, z, a, T(:, 21:35))
        case (5)
            call T0_damp_thole(x, y, z, a, T(:, 1:1))
            call T1_damp_thole(x, y, z, a, T(:, 2:4))
            call T2_damp_thole(x, y, z, a, T(:, 5:10))
            call T3_damp_thole(x, y, z, a, T(:, 11:20))
            call T4_damp_thole(x, y, z, a, T(:, 21:35))
            call T5_damp_thole(x, y, z, a, T(:, 36:56))
        case (6)
            call T0_damp_thole(x, y, z, a, T(:, 1:1))
            call T1_damp_thole(x, y, z, a, T(:, 2:4))
            call T2_damp_thole(x, y, z, a, T(:, 5:10))
            call T3_damp_thole(x, y, z, a, T(:, 11:20))
            call T4_damp_thole(x, y, z, a, T(:, 21:35))
            call T5_damp_thole(x, y, z, a, T(:, 36:56))
            call T6_damp_thole(x, y, z, a, T(:, 57:84))
        case default
            call T0_damp_thole(x, y, z, a, T(:, 1:1))
            call T1_damp_thole(x, y, z, a, T(:, 2:4))
            call T2_damp_thole(x, y, z, a, T(:, 5:10))
            call T3_damp_thole(x, y, z, a, T(:, 11:20))
            call T4_damp_thole(x, y, z, a, T(:, 21:35))
            call T5_damp_thole(x, y, z, a, T(:, 36:56))
            call T6_damp_thole(x, y, z, a, T(:, 57:84))
        end select
    end subroutine Tn_damp_thole
    pure subroutine T0_damp_thole(x, y, z, a, T)
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
            T(i, 1) = -(0.5_rp * (ai * sqrt(xi**2 + yi**2 + zi**2) + 2.0_rp) * exp(-ai * sqrt(xi**2 + yi**2 + zi**2)) - 1.0_rp) *&
& (xi**2 + yi**2 + zi**2)**(-0.5_rp)
        end do
    end subroutine T0_damp_thole
    pure subroutine T1_damp_thole(x, y, z, a, T)
        real(rp), intent(in) :: x(:), y(size(x)), z(size(x))
        real(rp), intent(in) :: a(:)
        real(rp), intent(inout) :: T(size(x), 4)
        integer(ip) :: i
        real(rp) :: xi, yi, zi
        real(rp) :: ai
        real(rp) :: x4
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x4 = -0.5_rp * ai * (-ai * sqrt(xi**2 + yi**2 + zi**2) - 1.0_rp) * exp(-ai * sqrt(xi**2 + yi**2 + zi**2)) / (xi**2 +&
& yi**2 + zi**2) + 0.5_rp * ((ai * sqrt(xi**2 + yi**2 + zi**2) + 2.0_rp) * exp(-ai * sqrt(xi**2 + yi**2 + zi**2)) - 2.0_rp) * (xi*&
&*2 + yi**2 + zi**2)**(-1.5_rp)
            T(i, 1) = x4 * xi
            T(i, 2) = x4 * yi
            T(i, 3) = x4 * zi
        end do
    end subroutine T1_damp_thole
    pure subroutine T2_damp_thole(x, y, z, a, T)
        real(rp), intent(in) :: x(:), y(size(x)), z(size(x))
        real(rp), intent(in) :: a(:)
        real(rp), intent(inout) :: T(size(x), 10)
        integer(ip) :: i
        real(rp) :: xi, yi, zi
        real(rp) :: ai
        real(rp) :: x0
        real(rp) :: x3
        real(rp) :: x5
        real(rp) :: x8
        real(rp) :: x10
        real(rp) :: x11
        real(rp) :: x12
        real(rp) :: x14
        real(rp) :: x17
        real(rp) :: x18
        real(rp) :: x20
        real(rp) :: x21
        real(rp) :: x22
        do i = 1, size(x)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ai = a(i)
            x0 = xi**2
            x3 = x0 + yi**2 + zi**2
            x5 = ai * sqrt(x3)
            x8 = ai * (x5 + 1.0_rp) * exp(-x5) / x3**2
            x10 = 3.0_rp / x3
            x11 = x3**(-1.5_rp)
            x12 = x5 + 2.0_rp
            x14 = 0.5_rp * x11 * (x12 * exp(-x5) - 2.0_rp)
            x17 = 2.0_rp * ai / x3
            x18 = ai * x12 / x3
            x20 = -x12 * x3**(-0.5_rp) + x3**(-0.5_rp)
            x21 = 0.5_rp * ai * x3**(-0.5_rp) * exp(-x5)
            x22 = x11 * x12
            T(i, 1) = -x0 * x8 - x14 * (x0 * x10 - 1.0_rp) - x21 * (x0 * x11 * x12 - x0 * x11 - x0 * x17 + x0 * x18 + x20)
            T(i, 2) = -xi * yi * (x21 * (-x11 - x17 + x18 + x22) + 1.5_rp * x3**(-2.5_rp) * (x12 * exp(-x5) - 2.0_rp) + x8)
            T(i, 3) = -xi * zi * (x21 * (-x11 - x17 + x18 + x22) + 1.5_rp * x3**(-2.5_rp) * (x12 * exp(-x5) - 2.0_rp) + x8)
            T(i, 4) = -x14 * (x10 * yi**2 - 1.0_rp) - x21 * (-x11 * yi**2 - x17 * yi**2 + x18 * yi**2 + x20 + x22 * yi**2) - x8&
                     & * yi**2
            T(i, 5) = -yi * zi * (x21 * (-x11 - x17 + x18 + x22) + 1.5_rp * x3**(-2.5_rp) * (x12 * exp(-x5) - 2.0_rp) + x8)
            T(i, 6) = -x14 * (x10 * zi**2 - 1.0_rp) - x21 * (-x11 * zi**2 - x17 * zi**2 + x18 * zi**2 + x20 + x22 * zi**2) - x8&
                     & * zi**2
        end do
    end subroutine T2_damp_thole
    pure subroutine T3_damp_thole(x, y, z, a, T)
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
        real(rp) :: x8
        real(rp) :: x9
        real(rp) :: x11
        real(rp) :: x13
        real(rp) :: x14
        real(rp) :: x15
        real(rp) :: x17
        real(rp) :: x18
        real(rp) :: x19
        real(rp) :: x20
        real(rp) :: x21
        real(rp) :: x23
        real(rp) :: x24
        real(rp) :: x25
        real(rp) :: x27
        real(rp) :: x29
        real(rp) :: x31
        real(rp) :: x32
        real(rp) :: x33
        real(rp) :: x34
        real(rp) :: x37
        real(rp) :: x38
        real(rp) :: x39
        real(rp) :: x41
        real(rp) :: x43
        real(rp) :: x46
        real(rp) :: x47
        real(rp) :: x48
        real(rp) :: x50
        real(rp) :: x56
        real(rp) :: x60
        real(rp) :: x66
        real(rp) :: x67
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
            x8 = ai * sqrt(x3)
            x9 = exp(-x8)
            x11 = ai / x3**2
            x13 = 3.0_rp * x11 * x9 * (x8 + 1.0_rp)
            x14 = x8 + 2.0_rp
            x15 = x14 * x9 - 2.0_rp
            x17 = 3.0_rp * x3**(-2.5_rp)
            x18 = x15 * x17
            x19 = x3**(-1.5_rp)
            x20 = 3.0_rp * x19
            x21 = x0 * x19
            x23 = ai * x4
            x24 = 2.0_rp * x23
            x25 = x14 * x23
            x27 = -x14 * x3**(-0.5_rp) + x3**(-0.5_rp)
            x29 = ai * x9
            x31 = x0 * x17
            x32 = 6.0_rp * x11
            x33 = ai**2
            x34 = x20 * x33
            x37 = 3.0_rp * x11 * x14
            x38 = -x14 * x20 + x20 + 6.0_rp * x23 - 3.0_rp * x25
            x39 = x29 * x3**(-0.5_rp)
            x41 = x29 * (x8 + 1.0_rp) / x3**3
            x43 = x15 * x3**(-3.5_rp)
            x46 = 0.5_rp * x11 * x9 * (x8 + 1.0_rp)
            x47 = 1.5_rp * x15 * x3**(-2.5_rp)
            x48 = x14 * x19
            x50 = x29 * (-x19 - x24 + x25 + x48)
            x56 = 0.5_rp * x39
            x60 = x1 * x17
            x66 = x1 * x19
            x67 = -x1 * x24 + x1 * x25 + x1 * x48 + x27 - x66
            T(i, 1) = 0.5_rp * xi * (x13 * (3.0_rp * x0 * x4 - 1.0_rp) + x18 * (5.0_rp * x0 * x4 - 3.0_rp) + x20 * x29 * (-x0 * x24&
& + x0 * x25 + x14 * x21 - x21 + x27) + x39 * (-x0 * x32 - x0 * x34 + x0 * x37 + x14 * x21 * x33 + x14 * x31 - x31 + x38))
            T(i, 2) = yi * (3.0_rp * x0 * x41 + 3.0_rp * x0 * x43 + 0.5_rp * x19 * x29 * (-x0 * x24 + x0 * x25 + x14 * x21 - x21 &
&+ x27) + x21 * x50 + 0.5_rp * x33 * x4 * x9 * (-x0 * x24 + x0 * x25 + x14 * x21 - x21 + x27) + x46 * (3.0_rp * x0 * x4 - 1.0_rp) &
&+ x47 * (3.0_rp * x0 * x4 - 1.0_rp) - x56 * (-2.0_rp * x0 * x11 * x14 + 5.0_rp * x0 * x11 - x14 * x31 - x19 + x21 * x33 - x23 + x3&
&1 + x48))
            T(i, 3) = zi * (3.0_rp * x0 * x41 + 3.0_rp * x0 * x43 + 0.5_rp * x19 * x29 * (-x0 * x24 + x0 * x25 + x14 * x21 - x21 &
&+ x27) + x21 * x50 + 0.5_rp * x33 * x4 * x9 * (-x0 * x24 + x0 * x25 + x14 * x21 - x21 + x27) + x46 * (3.0_rp * x0 * x4 - 1.0_rp) &
&+ x47 * (3.0_rp * x0 * x4 - 1.0_rp) - x56 * (-2.0_rp * x0 * x11 * x14 + 5.0_rp * x0 * x11 - x14 * x31 - x19 + x21 * x33 - x23 + x3&
&1 + x48))
            T(i, 4) = xi * (4.5_rp * x1 * x41 + 7.5_rp * x1 * x43 + 0.5_rp * x19 * x29 * x67 - x46 - x47 + x50 * x66 + x56 * (-x1&
                                                                                                                          & * x32 + x1 * x33 * x48 - x1 * x34 + x1 * x37 + x14 * x60 + x19 + x24 - x25 - x48 - x60))
            T(i, 5) = 0.5_rp * xi * yi * zi * (x20 * x50 + x39 * (x14 * x17 - x17 - x32 + x33 * x48 - x34 + x37) + 9.0_rp * x41 &
&+ 15.0_rp * x43)
            T(i, 6) = xi * (x19 * x2 * x50 + 0.5_rp * x19 * x29 * (-x19 * x2 - x2 * x24 + x2 * x25 + x2 * x48 + x27) + 4.5_rp * &
&x2 * x41 + 7.5_rp * x2 * x43 - x46 - x47 + x56 * (x14 * x17 * x2 - x17 * x2 + x19 - x2 * x32 + x2 * x33 * x48 - x2 * x34 + x2 *&
& x37 + x24 - x25 - x48))
            T(i, 7) = 0.5_rp * yi * (x13 * (3.0_rp * x1 * x4 - 1.0_rp) + x18 * (5.0_rp * x1 * x4 - 3.0_rp) + x20 * x29 * x67 + x39 &
                                   &* (-x1 * x32 + x1 * x33 * x48 - x1 * x34 + x1 * x37 + x14 * x60 + x38 - x60))
            T(i, 8) = zi * (3.0_rp * x1 * x41 + 3.0_rp * x1 * x43 + 0.5_rp * x19 * x29 * x67 + 0.5_rp * x33 * x4 * x67 * x9 + x46 &
&* (3.0_rp * x1 * x4 - 1.0_rp) + x47 * (3.0_rp * x1 * x4 - 1.0_rp) + x50 * x66 - x56 * (-2.0_rp * x1 * x11 * x14 + 5.0_rp * x1 * x11 &
&- x14 * x60 - x19 - x23 + x33 * x66 + x48 + x60))
            T(i, 9) = yi * (x19 * x2 * x50 + 0.5_rp * x19 * x29 * (-x19 * x2 - x2 * x24 + x2 * x25 + x2 * x48 + x27) + 4.5_rp * &
&x2 * x41 + 7.5_rp * x2 * x43 - x46 - x47 + x56 * (x14 * x17 * x2 - x17 * x2 + x19 - x2 * x32 + x2 * x33 * x48 - x2 * x34 + x2 *&
& x37 + x24 - x25 - x48))
            T(i, 10) = 0.5_rp * zi * (x13 * (3.0_rp * x2 * x4 - 1.0_rp) + x18 * (5.0_rp * x2 * x4 - 3.0_rp) + x20 * x29 * (-x19 * x&
&2 - x2 * x24 + x2 * x25 + x2 * x48 + x27) + x39 * (x14 * x17 * x2 - x17 * x2 - x2 * x32 + x2 * x33 * x48 - x2 * x34 + x2 * x37&
& + x38))
        end do
    end subroutine T3_damp_thole
    pure subroutine T4_damp_thole(x, y, z, a, T)
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
        real(rp) :: x8
        real(rp) :: x10
        real(rp) :: x11
        real(rp) :: x12
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
        real(rp) :: x28
        real(rp) :: x29
        real(rp) :: x30
        real(rp) :: x32
        real(rp) :: x33
        real(rp) :: x35
        real(rp) :: x36
        real(rp) :: x37
        real(rp) :: x38
        real(rp) :: x39
        real(rp) :: x41
        real(rp) :: x42
        real(rp) :: x43
        real(rp) :: x44
        real(rp) :: x45
        real(rp) :: x47
        real(rp) :: x50
        real(rp) :: x52
        real(rp) :: x53
        real(rp) :: x55
        real(rp) :: x57
        real(rp) :: x59
        real(rp) :: x60
        real(rp) :: x62
        real(rp) :: x64
        real(rp) :: x65
        real(rp) :: x67
        real(rp) :: x69
        real(rp) :: x70
        real(rp) :: x71
        real(rp) :: x72
        real(rp) :: x73
        real(rp) :: x74
        real(rp) :: x75
        real(rp) :: x76
        real(rp) :: x79
        real(rp) :: x82
        real(rp) :: x85
        real(rp) :: x86
        real(rp) :: x87
        real(rp) :: x89
        real(rp) :: x91
        real(rp) :: x93
        real(rp) :: x94
        real(rp) :: x96
        real(rp) :: x100
        real(rp) :: x102
        real(rp) :: x104
        real(rp) :: x107
        real(rp) :: x108
        real(rp) :: x111
        real(rp) :: x117
        real(rp) :: x118
        real(rp) :: x120
        real(rp) :: x121
        real(rp) :: x122
        real(rp) :: x123
        real(rp) :: x127
        real(rp) :: x128
        real(rp) :: x129
        real(rp) :: x131
        real(rp) :: x138
        real(rp) :: x139
        real(rp) :: x142
        real(rp) :: x144
        real(rp) :: x147
        real(rp) :: x148
        real(rp) :: x153
        real(rp) :: x156
        real(rp) :: x157
        real(rp) :: x160
        real(rp) :: x163
        real(rp) :: x164
        real(rp) :: x167
        real(rp) :: x172
        real(rp) :: x175
        real(rp) :: x180
        real(rp) :: x186
        real(rp) :: x187
        real(rp) :: x203
        real(rp) :: x204
        real(rp) :: x209
        real(rp) :: x210
        real(rp) :: x214
        real(rp) :: x226
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
            x7 = ai / x3**3
            x8 = x0 * x7
            x10 = ai * sqrt(x3)
            x11 = exp(-x10)
            x12 = x10 + 1.0_rp
            x15 = xi**4
            x16 = x3**(-2)
            x17 = 35.0_rp * x16
            x18 = x10 + 2.0_rp
            x19 = x11 * x18 - 2.0_rp
            x20 = x3**(-2.5_rp)
            x21 = 1.5_rp * x20
            x22 = x19 * x21
            x23 = 3.0_rp * x5 - 1.0_rp
            x24 = x3**(-1.5_rp)
            x25 = 3.0_rp * x24
            x26 = x0 * x24
            x28 = ai * x4
            x29 = 2.0_rp * x28
            x30 = x18 * x28
            x32 = -x18 * x3**(-0.5_rp) + x3**(-0.5_rp)
            x33 = -x0 * x29 + x0 * x30 + x18 * x26 - x26 + x32
            x35 = ai * x11 * x33
            x36 = ai * x16
            x37 = 6.0_rp * x36
            x38 = ai**2
            x39 = x25 * x38
            x41 = 3.0_rp * x36
            x42 = x18 * x41
            x43 = 3.0_rp * x20
            x44 = x0 * x43
            x45 = x18 * x44
            x47 = x18 * x25
            x50 = x25 + 6.0_rp * x28 - 3.0_rp * x30 - x47
            x52 = ai * x11
            x53 = x26 * x52
            x55 = 15.0_rp * x3**(-3.5_rp)
            x57 = 18.0_rp * x20
            x59 = 30.0_rp * x7
            x60 = x38 * x57
            x62 = ai**3 * x16
            x64 = x26 * x38
            x65 = x0 * x36
            x67 = 18.0_rp * x18 * x36
            x69 = x20 * x38
            x70 = 6.0_rp * x69
            x71 = 15.0_rp * x7
            x72 = -x25 - 6.0_rp * x28 + 3.0_rp * x30 + x47
            x73 = 0.5_rp * x11
            x74 = ai * x3**(-0.5_rp) * x73
            x75 = 9.0_rp * x12
            x76 = x3**(-4)
            x79 = x19 * x3**(-4.5_rp)
            x82 = x11 * x12 * x7
            x85 = x19 * x3**(-3.5_rp)
            x86 = 7.5_rp * x85
            x87 = x18 * x24
            x89 = -x24 - x29 + x30 + x87
            x91 = ai * x11 * x24
            x93 = 1.5_rp * x89 * x91
            x94 = 4.5_rp * x20
            x96 = 1.5_rp * x11 * x16 * x33 * x38
            x100 = x11 * (-2.0_rp * x18 * x65 - x24 - x28 + x44 - x45 + x64 + 5.0_rp * x65 + x87)
            x102 = 0.5_rp * x91
            x104 = x38 * x4 * x73
            x107 = x0 * x55
            x108 = x107 * x18
            x111 = 9.0_rp * x20
            x117 = 5.0_rp * x0 * x69
            x118 = 23.0_rp * x8
            x120 = 8.0_rp * x18 * x8
            x121 = x24 - x87
            x122 = x1 * x43
            x123 = x122 * x18
            x127 = 30.0_rp * x0 * x79
            x128 = x23 * x86
            x129 = x1 * x24
            x131 = x1 * x4
            x138 = x21 * x35
            x139 = -x1 * x29 + x1 * x30 + x1 * x87 - x129 + x32
            x142 = ai**3 * x33 * x73
            x144 = 21.0_rp * x0 * x12 * x52 * x76
            x147 = 6.0_rp * x0 * x20 * x52 * x89
            x148 = x1 * x7
            x153 = 3.0_rp * x82
            x156 = x38 * x87
            x157 = x18 * x43
            x160 = yi * zi
            x163 = x2 * x24
            x164 = x2 * x4
            x167 = x121 + x156 * x2 + x157 * x2 - x2 * x37 - x2 * x39 + x2 * x42 - x2 * x43 + x29 - x30
            x172 = x1 * x55
            x175 = x172 * x18
            x180 = x1 * x52
            x186 = x1 * x18 * x41 - x1 * x37 + x1 * x38 * x87 - x1 * x39 - x122 + x123 + x50
            x187 = x139 * x52
            x203 = x52 * (-x163 - x2 * x29 + x2 * x30 + x2 * x87 + x32)
            x204 = x2 * x52
            x209 = yi**4
            x210 = 3.0_rp * x131 - 1.0_rp
            x214 = x129 * x38
            x226 = zi**4
            T(i, 1) = -6.0_rp * x11 * x12 * x8 * (5.0_rp * x5 - 3.0_rp) - x22 * (x15 * x17 - 30.0_rp * x5 + 3.0_rp) - x23 * x25 * x&
&35 - 2.0_rp * x53 * (-x0 * x37 - x0 * x39 + x0 * x42 + x18 * x26 * x38 - x44 + x45 + x50) - x74 * (-x0 * x18 * x57 + x0 * x57 -&
& x0 * x67 + x15 * x18 * x55 + x15 * x18 * x62 + x15 * x18 * x70 + x15 * x18 * x71 - x15 * x55 - x15 * x59 - x15 * x60 - 4.0_rp &
&* x15 * x62 - 6.0_rp * x18 * x26 * x38 + 18.0_rp * x64 + 36.0_rp * x65 + x72)
            T(i, 2) = -xi * yi * (-1.5_rp * ai * x100 * x24 + x0 * x52 * x75 * x76 + 15.0_rp * x0 * x79 + x102 * (-x0 * x37 - x0&
& * x39 + x0 * x42 + x18 * x26 * x38 - x44 + x45 + x50) + x104 * (-x0 * x37 - x0 * x39 + x0 * x42 + x18 * x26 * x38 - x44 + x45&
& + x50) + 4.5_rp * x23 * x82 + x23 * x93 + x35 * x94 - x74 * (-12.0_rp * x0 * x18 * x7 + x0 * x62 + 12.0_rp * x0 * x69 + x107 - x&
&108 + x111 * x18 - x111 + x18 * x37 - 15.0_rp * x36 - x38 * x45 - x39 + 27.0_rp * x8) + 1.5_rp * x82 * (5.0_rp * x5 - 3.0_rp) + x86&
& * (5.0_rp * x5 - 3.0_rp) + x96)
            T(i, 3) = -xi * zi * (-1.5_rp * ai * x100 * x24 + x0 * x52 * x75 * x76 + 15.0_rp * x0 * x79 + x102 * (-x0 * x37 - x0&
& * x39 + x0 * x42 + x18 * x26 * x38 - x44 + x45 + x50) + x104 * (-x0 * x37 - x0 * x39 + x0 * x42 + x18 * x26 * x38 - x44 + x45&
& + x50) + 4.5_rp * x23 * x82 + x23 * x93 + x35 * x94 - x74 * (-12.0_rp * x0 * x18 * x7 + x0 * x62 + 12.0_rp * x0 * x69 + x107 - x&
&108 + x111 * x18 - x111 + x18 * x37 - 15.0_rp * x36 - x38 * x45 - x39 + 27.0_rp * x8) + 1.5_rp * x82 * (5.0_rp * x5 - 3.0_rp) + x86&
& * (5.0_rp * x5 - 3.0_rp) + x96)
            T(i, 4) = ai * x100 * x129 + 3.0_rp * x0 * x85 - x1 * x127 - x1 * x128 - x1 * x138 - x1 * x144 - x1 * x147 - x1 * x&
&96 + x100 * x131 * x38 - x102 * x139 * x23 + x102 * x33 + x104 * x33 - 3.0_rp * x11 * x12 * x148 * x23 + 3.0_rp * x11 * x12 * x8&
& - x129 * x142 + x22 * x23 - x53 * (x1 * x18 * x41 - x1 * x37 + x1 * x38 * x87 - x1 * x39 + x121 - x122 + x123 + x29 - x30) + &
&x74 * (x1 * x107 - x1 * x108 + x1 * x117 + x1 * x118 - x1 * x120 - x1 * x41 + x121 - x122 + x123 + 2.0_rp * x18 * x65 + x28 - x&
&44 + x45 - x64 - 5.0_rp * x65)
            T(i, 5) = -x160 * (-ai * x100 * x24 - x100 * x38 * x4 + x102 * x23 * x89 + x127 + x128 + x138 + x142 * x24 + x144 &
&+ x147 + x153 * x23 + x53 * (x156 + x157 - x37 - x39 + x42 - x43) - x74 * (x107 - x108 + x117 + x118 - x120 + x157 - x41 - x43&
&) + x96)
            T(i, 6) = ai * x100 * x163 + 3.0_rp * x0 * x85 + x100 * x164 * x38 - x102 * x23 * (-x163 - x2 * x29 + x2 * x30 + x2&
& * x87 + x32) + x102 * x33 + x104 * x33 + 3.0_rp * x11 * x12 * x8 - x127 * x2 - x128 * x2 - x138 * x2 - x142 * x163 - x144 * x2&
& - x147 * x2 - x153 * x2 * x23 - x167 * x53 - x2 * x96 + x22 * x23 + x74 * (x107 * x2 - x108 * x2 + x117 * x2 + x118 * x2 - x1&
&20 * x2 + x121 + x157 * x2 + 2.0_rp * x18 * x65 - x2 * x41 - x2 * x43 + x28 - x44 + x45 - x64 - 5.0_rp * x65)
            T(i, 7) = -xi * yi * (52.5_rp * x1 * x79 + x102 * x186 - x11 * x7 * x75 + 30.0_rp * x12 * x180 * x76 + x180 * x89 * &
&x94 + x187 * x94 + x74 * (x1 * x18 * x62 + 6.0_rp * x1 * x18 * x69 + 15.0_rp * x1 * x18 * x7 - x1 * x38 * x57 - x1 * x59 - 4.0_rp&
& * x1 * x62 - x111 * x18 + x111 - x172 + x175 - 9.0_rp * x18 * x36 + 9.0_rp * x24 * x38 + 18.0_rp * x36 - x38 * x47) - 22.5_rp * x&
&85 + 1.5_rp * x91 * (x1 * x18 * x41 - x1 * x37 + x1 * x38 * x87 - x1 * x39 + x121 - x122 + x123 + x29 - x30) - x93)
            T(i, 8) = -xi * zi * (52.5_rp * x1 * x79 - x102 * x89 + 30.0_rp * x12 * x180 * x76 + x129 * x52 * (x156 + x157 - x37&
& - x39 + x42 - x43) - x153 + 7.5_rp * x180 * x20 * x89 + x187 * x21 + x74 * (x1 * x18 * x62 + 6.0_rp * x1 * x18 * x69 + 15.0_rp *&
& x1 * x18 * x7 - x1 * x38 * x57 - x1 * x59 - 4.0_rp * x1 * x62 - x156 - x157 - x172 + x175 + x37 + x39 - x42 + x43) - x86 + x91&
& * (x1 * x18 * x41 - x1 * x37 + x1 * x38 * x87 - x1 * x39 + x121 - x122 + x123 + x29 - x30))
            T(i, 9) = -xi * yi * (-x102 * x89 + 30.0_rp * x12 * x204 * x76 - x153 + x163 * x52 * (x156 + x157 - x37 - x39 + x42&
& - x43) + x167 * x91 + 52.5_rp * x2 * x79 + 7.5_rp * x20 * x204 * x89 + x203 * x21 + x74 * (-x156 - x157 + x18 * x2 * x55 + x18 &
&* x2 * x62 + x18 * x2 * x70 + x18 * x2 * x71 - x2 * x38 * x57 - x2 * x55 - x2 * x59 - 4.0_rp * x2 * x62 + x37 + x39 - x42 + x43&
&) - x86)
            T(i, 10) = -xi * zi * (x102 * (x156 * x2 + x157 * x2 - x2 * x37 - x2 * x39 + x2 * x42 - x2 * x43 + x50) - x11 * x7&
& * x75 + 30.0_rp * x12 * x204 * x76 + 1.5_rp * x167 * x91 + 52.5_rp * x2 * x79 + x203 * x94 + x204 * x89 * x94 + x74 * (-x111 * x&
&18 + x111 + x18 * x2 * x55 + x18 * x2 * x62 + x18 * x2 * x70 + x18 * x2 * x71 - 9.0_rp * x18 * x36 - x2 * x38 * x57 - x2 * x55 &
&- x2 * x59 - 4.0_rp * x2 * x62 + 9.0_rp * x24 * x38 + 18.0_rp * x36 - x38 * x47) - 22.5_rp * x85 - x93)
            T(i, 11) = -6.0_rp * x11 * x12 * x148 * (5.0_rp * x131 - 3.0_rp) - 2.0_rp * x129 * x186 * x52 - x187 * x210 * x25 - x2&
&2 * (-30.0_rp * x131 + x17 * x209 + 3.0_rp) - x74 * (-x1 * x18 * x57 + 36.0_rp * x1 * x36 - 6.0_rp * x1 * x38 * x87 + x1 * x57 - x&
&1 * x67 + x18 * x209 * x55 + x18 * x209 * x62 + x18 * x209 * x70 + x18 * x209 * x71 - x209 * x55 - x209 * x59 - x209 * x60 - 4&
&.0_rp * x209 * x62 + 18.0_rp * x214 + x72)
            T(i, 12) = -x160 * (15.0_rp * x1 * x79 + x102 * x186 + x104 * x186 + 1.5_rp * x11 * x139 * x16 * x38 + x180 * x75 * &
&x76 + x187 * x94 + 4.5_rp * x210 * x82 + x210 * x93 - x74 * (x1 * x62 + 12.0_rp * x1 * x69 + x111 * x18 - x111 - x123 * x38 - 12&
&.0_rp * x148 * x18 + 27.0_rp * x148 + x172 - x175 + x18 * x37 - 15.0_rp * x36 - x39) + 1.5_rp * x82 * (5.0_rp * x131 - 3.0_rp) + x86&
& * (5.0_rp * x131 - 3.0_rp) - 1.5_rp * x91 * (-2.0_rp * x1 * x18 * x36 + 5.0_rp * x1 * x36 + x122 - x123 + x214 - x24 - x28 + x87))
            T(i, 13) = -ai**3 * x139 * x163 * x73 - 21.0_rp * x1 * x12 * x2 * x52 * x76 - 6.0_rp * x1 * x2 * x20 * x52 * x89 - 3&
&0.0_rp * x1 * x2 * x79 + 3.0_rp * x1 * x85 + x102 * x139 - x102 * x210 * (-x163 - x2 * x29 + x2 * x30 + x2 * x87 + x32) + x104 *&
& x139 + 3.0_rp * x11 * x12 * x148 - 1.5_rp * x11 * x139 * x16 * x2 * x38 + x11 * x164 * x38 * (-2.0_rp * x1 * x18 * x36 + 5.0_rp *&
& x1 * x36 + x122 - x123 + x214 - x24 - x28 + x87) - x129 * x167 * x52 - x153 * x2 * x210 + x163 * x52 * (-2.0_rp * x1 * x18 * x&
&36 + 5.0_rp * x1 * x36 + x122 - x123 + x214 - x24 - x28 + x87) - x187 * x2 * x21 - x2 * x210 * x86 + x210 * x22 + x74 * (2.0_rp &
&* x1 * x18 * x36 + 5.0_rp * x1 * x2 * x69 - 5.0_rp * x1 * x36 + x121 - x122 + x123 - 8.0_rp * x148 * x18 * x2 + 23.0_rp * x148 * x&
&2 + x157 * x2 + x172 * x2 - x175 * x2 - x2 * x41 - x2 * x43 - x214 + x28)
            T(i, 14) = -x160 * (x102 * (x156 * x2 + x157 * x2 - x2 * x37 - x2 * x39 + x2 * x42 - x2 * x43 + x50) - x11 * x7 * &
&x75 + 30.0_rp * x12 * x204 * x76 + 1.5_rp * x167 * x91 + 52.5_rp * x2 * x79 + x203 * x94 + x204 * x89 * x94 + x74 * (-x111 * x18 &
&+ x111 + x18 * x2 * x55 + x18 * x2 * x62 + x18 * x2 * x70 + x18 * x2 * x71 - 9.0_rp * x18 * x36 - x2 * x38 * x57 - x2 * x55 - x&
&2 * x59 - 4.0_rp * x2 * x62 + 9.0_rp * x24 * x38 + 18.0_rp * x36 - x38 * x47) - 22.5_rp * x85 - x93)
            T(i, 15) = -2.0_rp * x163 * x52 * (x156 * x2 + x157 * x2 - x2 * x37 - x2 * x39 + x2 * x42 - x2 * x43 + x50) - 6.0_rp&
& * x2 * x82 * (5.0_rp * x164 - 3.0_rp) - x203 * x25 * (3.0_rp * x164 - 1.0_rp) - x22 * (-30.0_rp * x164 + x17 * x226 + 3.0_rp) - x74&
& * (-6.0_rp * x156 * x2 + 18.0_rp * x163 * x38 - x18 * x2 * x57 + x18 * x226 * x55 + x18 * x226 * x62 + x18 * x226 * x70 + x18 *&
& x226 * x71 + 36.0_rp * x2 * x36 + x2 * x57 - x2 * x67 - x226 * x55 - x226 * x59 - x226 * x60 - 4.0_rp * x226 * x62 + x72)
        end do
    end subroutine T4_damp_thole
    pure subroutine T5_damp_thole(x, y, z, a, T)
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
        real(rp) :: x11
        real(rp) :: x13
        real(rp) :: x14
        real(rp) :: x15
        real(rp) :: x16
        real(rp) :: x18
        real(rp) :: x19
        real(rp) :: x20
        real(rp) :: x21
        real(rp) :: x23
        real(rp) :: x24
        real(rp) :: x25
        real(rp) :: x27
        real(rp) :: x28
        real(rp) :: x29
        real(rp) :: x31
        real(rp) :: x32
        real(rp) :: x33
        real(rp) :: x34
        real(rp) :: x35
        real(rp) :: x36
        real(rp) :: x38
        real(rp) :: x39
        real(rp) :: x40
        real(rp) :: x41
        real(rp) :: x42
        real(rp) :: x43
        real(rp) :: x45
        real(rp) :: x46
        real(rp) :: x47
        real(rp) :: x48
        real(rp) :: x49
        real(rp) :: x51
        real(rp) :: x54
        real(rp) :: x55
        real(rp) :: x57
        real(rp) :: x58
        real(rp) :: x59
        real(rp) :: x61
        real(rp) :: x63
        real(rp) :: x64
        real(rp) :: x66
        real(rp) :: x68
        real(rp) :: x70
        real(rp) :: x71
        real(rp) :: x72
        real(rp) :: x73
        real(rp) :: x74
        real(rp) :: x75
        real(rp) :: x76
        real(rp) :: x77
        real(rp) :: x80
        real(rp) :: x81
        real(rp) :: x82
        real(rp) :: x83
        real(rp) :: x85
        real(rp) :: x86
        real(rp) :: x88
        real(rp) :: x91
        real(rp) :: x93
        real(rp) :: x94
        real(rp) :: x96
        real(rp) :: x97
        real(rp) :: x98
        real(rp) :: x99
        real(rp) :: x100
        real(rp) :: x101
        real(rp) :: x102
        real(rp) :: x106
        real(rp) :: x107
        real(rp) :: x108
        real(rp) :: x109
        real(rp) :: x110
        real(rp) :: x112
        real(rp) :: x113
        real(rp) :: x117
        real(rp) :: x118
        real(rp) :: x119
        real(rp) :: x120
        real(rp) :: x122
        real(rp) :: x123
        real(rp) :: x125
        real(rp) :: x126
        real(rp) :: x129
        real(rp) :: x131
        real(rp) :: x132
        real(rp) :: x134
        real(rp) :: x135
        real(rp) :: x140
        real(rp) :: x141
        real(rp) :: x143
        real(rp) :: x150
        real(rp) :: x151
        real(rp) :: x159
        real(rp) :: x160
        real(rp) :: x161
        real(rp) :: x162
        real(rp) :: x166
        real(rp) :: x168
        real(rp) :: x170
        real(rp) :: x172
        real(rp) :: x173
        real(rp) :: x174
        real(rp) :: x175
        real(rp) :: x178
        real(rp) :: x179
        real(rp) :: x180
        real(rp) :: x181
        real(rp) :: x182
        real(rp) :: x183
        real(rp) :: x184
        real(rp) :: x186
        real(rp) :: x191
        real(rp) :: x192
        real(rp) :: x193
        real(rp) :: x196
        real(rp) :: x198
        real(rp) :: x199
        real(rp) :: x201
        real(rp) :: x203
        real(rp) :: x204
        real(rp) :: x205
        real(rp) :: x206
        real(rp) :: x210
        real(rp) :: x213
        real(rp) :: x214
        real(rp) :: x217
        real(rp) :: x218
        real(rp) :: x220
        real(rp) :: x221
        real(rp) :: x223
        real(rp) :: x225
        real(rp) :: x226
        real(rp) :: x228
        real(rp) :: x229
        real(rp) :: x230
        real(rp) :: x231
        real(rp) :: x232
        real(rp) :: x233
        real(rp) :: x234
        real(rp) :: x235
        real(rp) :: x236
        real(rp) :: x237
        real(rp) :: x238
        real(rp) :: x240
        real(rp) :: x241
        real(rp) :: x244
        real(rp) :: x245
        real(rp) :: x246
        real(rp) :: x249
        real(rp) :: x251
        real(rp) :: x254
        real(rp) :: x257
        real(rp) :: x259
        real(rp) :: x261
        real(rp) :: x262
        real(rp) :: x266
        real(rp) :: x268
        real(rp) :: x271
        real(rp) :: x272
        real(rp) :: x273
        real(rp) :: x274
        real(rp) :: x275
        real(rp) :: x278
        real(rp) :: x280
        real(rp) :: x284
        real(rp) :: x288
        real(rp) :: x292
        real(rp) :: x298
        real(rp) :: x299
        real(rp) :: x301
        real(rp) :: x302
        real(rp) :: x303
        real(rp) :: x304
        real(rp) :: x312
        real(rp) :: x313
        real(rp) :: x314
        real(rp) :: x316
        real(rp) :: x317
        real(rp) :: x321
        real(rp) :: x322
        real(rp) :: x323
        real(rp) :: x327
        real(rp) :: x331
        real(rp) :: x332
        real(rp) :: x343
        real(rp) :: x345
        real(rp) :: x347
        real(rp) :: x357
        real(rp) :: x361
        real(rp) :: x363
        real(rp) :: x366
        real(rp) :: x368
        real(rp) :: x370
        real(rp) :: x371
        real(rp) :: x375
        real(rp) :: x380
        real(rp) :: x382
        real(rp) :: x388
        real(rp) :: x389
        real(rp) :: x391
        real(rp) :: x392
        real(rp) :: x395
        real(rp) :: x396
        real(rp) :: x398
        real(rp) :: x401
        real(rp) :: x405
        real(rp) :: x406
        real(rp) :: x407
        real(rp) :: x409
        real(rp) :: x417
        real(rp) :: x419
        real(rp) :: x424
        real(rp) :: x428
        real(rp) :: x430
        real(rp) :: x434
        real(rp) :: x436
        real(rp) :: x443
        real(rp) :: x458
        real(rp) :: x459
        real(rp) :: x464
        real(rp) :: x465
        real(rp) :: x467
        real(rp) :: x475
        real(rp) :: x485
        real(rp) :: x487
        real(rp) :: x491
        real(rp) :: x492
        real(rp) :: x496
        real(rp) :: x503
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
            x6 = xi**4
            x7 = x3**(-2)
            x8 = x6 * x7
            x10 = x3**(-3)
            x11 = ai * x10
            x13 = ai * sqrt(x3)
            x14 = exp(-x13)
            x15 = x13 + 1.0_rp
            x16 = x14 * x15
            x18 = 7.5_rp * x11 * x16
            x19 = x3**(-3.5_rp)
            x20 = x13 + 2.0_rp
            x21 = x14 * x20 - 2.0_rp
            x23 = 7.5_rp * x19 * x21
            x24 = x3**(-1.5_rp)
            x25 = x0 * x24
            x27 = ai * x4
            x28 = 2.0_rp * x27
            x29 = x20 * x27
            x31 = -x20 * x3**(-0.5_rp) + x3**(-0.5_rp)
            x32 = -x0 * x28 + x0 * x29 + x20 * x25 - x25 + x31
            x33 = x14 * x32
            x34 = x3**(-2.5_rp)
            x35 = ai * x34
            x36 = 5.0_rp * x5 - 3.0_rp
            x38 = 3.0_rp * x5 - 1.0_rp
            x39 = ai * x7
            x40 = 6.0_rp * x39
            x41 = ai**2
            x42 = 3.0_rp * x24
            x43 = x41 * x42
            x45 = 3.0_rp * x39
            x46 = x20 * x45
            x47 = 3.0_rp * x34
            x48 = x0 * x47
            x49 = x20 * x48
            x51 = x20 * x42
            x54 = 6.0_rp * x27 - 3.0_rp * x29 + x42 - x51
            x55 = -x0 * x40 - x0 * x43 + x0 * x46 + x20 * x25 * x41 - x48 + x49 + x54
            x57 = ai * x14 * x24
            x58 = 5.0_rp * x57
            x59 = 15.0_rp * x19
            x61 = 18.0_rp * x34
            x63 = 30.0_rp * x11
            x64 = x41 * x6
            x66 = ai**3
            x68 = x25 * x41
            x70 = 18.0_rp * x39
            x71 = x20 * x70
            x72 = x34 * x41
            x73 = 6.0_rp * x20
            x74 = x72 * x73
            x75 = x11 * x20
            x76 = 15.0_rp * x75
            x77 = -6.0_rp * x27 + 3.0_rp * x29 - x42 + x51
            x80 = x3**(-4.5_rp)
            x81 = 105.0_rp * x80
            x82 = x6 * x81
            x83 = x0 * x19
            x85 = ai / x3**4
            x86 = x6 * x85
            x88 = x10 * x66
            x91 = ai**4 * x34
            x93 = x66 * x7
            x94 = x0 * x93
            x96 = x0 * x72
            x97 = x0 * x11
            x98 = x20 * x97
            x99 = 60.0_rp * x20
            x100 = 10.0_rp * x20
            x101 = 45.0_rp * x19
            x102 = x101 * x20
            x106 = x20 * x24
            x107 = x106 * x41
            x108 = x20 * x39
            x109 = 15.0_rp * x107 + 45.0_rp * x108 + 45.0_rp * x20 * x34 - 45.0_rp * x24 * x41 - 45.0_rp * x34 - 90.0_rp * x39
            x110 = 0.5_rp * x14
            x112 = ai * x110 * x3**(-0.5_rp)
            x113 = ai * x14
            x117 = x0 * x85
            x118 = x117 * x16
            x119 = x21 * x80
            x120 = x0 * x119
            x122 = 1.5_rp * x14
            x123 = x11 * x122 * x15
            x125 = x106 - x24 - x28 + x29
            x126 = x14 * x35
            x129 = ai * x33
            x131 = 9.0_rp * x34
            x132 = ai * x131
            x134 = x41 * x7
            x135 = 3.0_rp * x134
            x140 = -2.0_rp * x0 * x108 + 5.0_rp * x0 * x39 + x106 - x24 - x27 + x48 - x49 + x68
            x141 = x14 * x140
            x143 = 6.0_rp * x126
            x150 = 15.0_rp * x83
            x151 = x150 * x20
            x159 = x14 * (x131 * x20 - x131 + x150 - x151 + x20 * x40 - 15.0_rp * x39 - x41 * x49 - x43 + x94 + 12.0_rp * x96 + &
&27.0_rp * x97 - 12.0_rp * x98)
            x160 = ai * x159
            x161 = x4 * x41
            x162 = x110 * x161
            x166 = -x131 * x20 + x131
            x168 = x166 - x20 * x40 + 15.0_rp * x39 + x43
            x170 = x21 * x3**(-5.5_rp)
            x172 = 210.0_rp * x0 * x170
            x173 = x1 * x101
            x174 = 69.0_rp * x11
            x175 = 15.0_rp * x72
            x178 = 7.0_rp * x0 * x88
            x179 = 24.0_rp * x75
            x180 = x41 * x83
            x181 = 72.0_rp * x180
            x182 = 177.0_rp * x117
            x183 = x117 * x20
            x184 = 72.0_rp * x183
            x186 = x0 * x81
            x191 = 52.5_rp * x119
            x192 = x191 * x36
            x193 = x1 * x24
            x196 = x1 * x4
            x198 = x132 * x141
            x199 = x135 * x141
            x201 = x122 * x35 * x55
            x203 = x1 * x106 - x1 * x28 + x1 * x29 - x193 + x31
            x204 = x122 * x35
            x205 = x203 * x204
            x206 = -x106 + x24
            x210 = x1 * x106 * x41
            x213 = x1 * x47
            x214 = x20 * x213
            x217 = x1 * x20 * x45 - x1 * x40 - x1 * x43 + x206 + x210 - x213 + x214 + x28 - x29
            x218 = 1.5_rp * x57
            x220 = x122 * x134 * x55
            x221 = 1.5_rp * x33
            x223 = x221 * x34 * x66
            x225 = x10 * x33 * x41
            x226 = 10.5_rp * x225
            x228 = 22.5_rp * x129 * x19
            x229 = x132 * x14
            x230 = x125 * x229
            x231 = x230 * x38
            x232 = x16 * x85
            x233 = 15.0_rp * x232 * x36
            x234 = x113 * x125
            x235 = 18.0_rp * x234 * x83
            x236 = x0 * x1
            x237 = x113 * x15 / x3**5
            x238 = 120.0_rp * x237
            x240 = 22.5_rp * x232
            x241 = x240 * x38
            x244 = 5.0_rp * x96
            x245 = 23.0_rp * x97
            x246 = 8.0_rp * x98
            x249 = 2.0_rp * x0 * x108 - 5.0_rp * x0 * x39 + x1 * x150 - x1 * x151 + x1 * x244 + x1 * x245 - x1 * x246 - x1 * x45&
& + x206 - x213 + x214 + x27 - x48 + x49 - x68
            x251 = 4.5_rp * x11 * x16
            x254 = 0.5_rp * x57
            x257 = x20 * x47
            x259 = x107 + x257 - x40 - x43 + x46 - x47
            x261 = x150 - x151 + x244 + x245 - x246 + x257 - x45 - x47
            x262 = x24 * x66
            x266 = xi * yi * zi
            x268 = x101 * x2
            x271 = x2 * x24
            x272 = x271 * x66
            x273 = x2 * x4
            x274 = x106 * x2 - x2 * x28 + x2 * x29 - x271 + x31
            x275 = x204 * x274
            x278 = x107 * x2
            x280 = x2 * x47
            x284 = x2 * x257 - x2 * x40 - x2 * x43 + x2 * x46 + x206 + x278 + x28 - x280 - x29
            x288 = 2.0_rp * x0 * x108 - 5.0_rp * x0 * x39 + x150 * x2 - x151 * x2 + x2 * x244 + x2 * x245 - x2 * x246 + x2 * x25&
&7 - x2 * x45 + x206 + x27 - x280 - x48 + x49 - x68
            x292 = x1 * x11
            x298 = -x257 + x47
            x299 = -x150 + x151 - x244 - x245 + x246 + x45
            x301 = ai * x122 * x3**(-0.5_rp)
            x302 = x122 * x161
            x303 = x1 * x59
            x304 = x20 * x303
            x312 = x1 * x61
            x313 = x312 * x41
            x314 = x1 * x93
            x316 = x20 * x314
            x317 = x1 * x72
            x321 = -x1 * x63 + x1 * x76 - 9.0_rp * x108 + x166 + 9.0_rp * x24 * x41 - x303 + x304 - x313 - 4.0_rp * x314 + x316 +&
& x317 * x73 - x41 * x51 + x70
            x322 = x113 * x25
            x323 = x1 * x20 * x45 - x1 * x40 - x1 * x43 + x210 - x213 + x214 + x54
            x327 = x217 * x229
            x331 = 22.5_rp * x19 * x21
            x332 = 1.5_rp * x141
            x343 = 195.0_rp * x237
            x345 = x298 - x303 + x304
            x347 = -x1 * x63 + x1 * x76 - x107 - x313 - 4.0_rp * x314 + x316 + x317 * x73 + x345 + x40 + x43 - x46
            x357 = x11 * x2
            x361 = x2 * x59
            x363 = x20 * x361
            x366 = x2 * x61
            x368 = x2 * x93
            x370 = x20 * x368
            x371 = x2 * x72
            x375 = -x107 - x2 * x63 + x2 * x76 + x298 - x361 + x363 - x366 * x41 - 4.0_rp * x368 + x370 + x371 * x73 + x40 + x4&
&3 - x46
            x380 = x191 * x2
            x382 = 7.5_rp * x2
            x388 = -9.0_rp * x108 + x166 - x2 * x63 + x2 * x76 + 9.0_rp * x24 * x41 - x361 + x363 - x366 * x41 - 4.0_rp * x368 + &
&x370 + x371 * x73 - x41 * x51 + x70
            x389 = x2 * x257 - x2 * x40 - x2 * x43 + x2 * x46 + x278 - x280 + x54
            x391 = x229 * x284
            x392 = 4.5_rp * x126
            x395 = yi**4
            x396 = 472.5_rp * x170
            x398 = 90.0_rp * x1 * x19
            x401 = 90.0_rp * x11
            x405 = x395 * x81
            x406 = 210.0_rp * x85
            x407 = x395 * x41
            x409 = 40.0_rp * x88
            x417 = x113 * x42
            x419 = x113 * x193
            x424 = 262.5_rp * x237
            x428 = x193 * x41
            x430 = 36.0_rp * x1 * x39 - x1 * x71 - x20 * x312 + x20 * x395 * x59 + x20 * x395 * x93 - 6.0_rp * x210 + x312 - x39&
&5 * x59 - x395 * x63 + x395 * x74 + x395 * x76 - 4.0_rp * x395 * x93 - x407 * x61 + 18.0_rp * x428 + x77
            x434 = x1 * x81
            x436 = 135.0_rp * x19 * x41
            x443 = x100 * x88
            x458 = x1 * x2
            x459 = x19 * x458
            x464 = x2 * x392
            x465 = x113 * x271
            x467 = x2 * x304
            x475 = zi**4
            x485 = 5.0_rp * x196 - 3.0_rp
            x487 = 3.0_rp * x196 - 1.0_rp
            x491 = x1 * x75
            x492 = x1 * x119
            x496 = -2.0_rp * x1 * x108 + 5.0_rp * x1 * x39 + x106 + x213 - x214 - x24 - x27 + x428
            x503 = x1 * x232
            T(i, 1) = xi * (x112 * (x100 * x6 * x88 - x100 * x94 + x102 * x64 + x109 - 135.0_rp * x19 * x64 + x20 * x6 * x91 + &
&x20 * x82 - 150.0_rp * x20 * x83 + 105.0_rp * x20 * x86 - 40.0_rp * x6 * x88 - 5.0_rp * x6 * x91 - x82 + 150.0_rp * x83 - 210.0_rp *&
& x86 + 40.0_rp * x94 - x96 * x99 + 180.0_rp * x96 + 300.0_rp * x97 - 150.0_rp * x98) + x18 * (-30.0_rp * x5 + 35.0_rp * x8 + 3.0_rp) &
&+ x23 * (-70.0_rp * x5 + 63.0_rp * x8 + 15.0_rp) + 15.0_rp * x33 * x35 * x36 + x38 * x55 * x58 + 2.5_rp * x57 * (-x0 * x20 * x61 + &
&36.0_rp * x0 * x39 + x0 * x61 - x0 * x71 - 6.0_rp * x20 * x25 * x41 + x20 * x59 * x6 + x20 * x66 * x8 - x59 * x6 - x6 * x63 + x6&
& * x74 + x6 * x76 - x61 * x64 - 4.0_rp * x66 * x8 + 18.0_rp * x68 + x77))
            T(i, 2) = yi * (-ai * x141 * x38 * x42 + 6.0_rp * x0 * x125 * x126 * x36 + 2.0_rp * x0 * x134 * x14 * x55 + x0 * x14&
&3 * x55 - x112 * (x0 * x20 * x41 * x61 + x168 - 30.0_rp * x19 * x20 * x64 + 105.0_rp * x19 * x64 - 4.0_rp * x20 * x6 * x88 - x20 &
&* x82 + 90.0_rp * x20 * x83 - 90.0_rp * x20 * x86 + 22.0_rp * x6 * x88 + x6 * x91 + x82 - 90.0_rp * x83 + 195.0_rp * x86 - 6.0_rp * &
&x94 - 72.0_rp * x96 - 162.0_rp * x97 + 72.0_rp * x98) + 60.0_rp * x113 * x15 * x6 / x3**5 + 30.0_rp * x118 * x36 + 30.0_rp * x120 * &
&(7.0_rp * x5 - 3.0_rp) + x123 * (-30.0_rp * x5 + 35.0_rp * x8 + 3.0_rp) + 18.0_rp * x129 * x83 + x132 * x33 * x38 + x135 * x33 * x38&
& - 2.0_rp * x160 * x25 + x162 * (-x0 * x20 * x61 + 36.0_rp * x0 * x39 + x0 * x61 - x0 * x71 - 6.0_rp * x20 * x25 * x41 + x20 * x5&
&9 * x6 + x20 * x66 * x8 - x59 * x6 - x6 * x63 + x6 * x74 + x6 * x76 - x61 * x64 - 4.0_rp * x66 * x8 + 18.0_rp * x68 + x77) + x23&
& * (-30.0_rp * x5 + 35.0_rp * x8 + 3.0_rp) + 0.5_rp * x57 * (-x0 * x20 * x61 + 36.0_rp * x0 * x39 + x0 * x61 - x0 * x71 - 6.0_rp * x&
&20 * x25 * x41 + x20 * x59 * x6 + x20 * x66 * x8 - x59 * x6 - x6 * x63 + x6 * x74 + x6 * x76 - x61 * x64 - 4.0_rp * x66 * x8 + &
&18.0_rp * x68 + x77))
            T(i, 3) = zi * (-ai * x141 * x38 * x42 + 6.0_rp * x0 * x125 * x126 * x36 + 2.0_rp * x0 * x134 * x14 * x55 + x0 * x14&
&3 * x55 - x112 * (x0 * x20 * x41 * x61 + x168 - 30.0_rp * x19 * x20 * x64 + 105.0_rp * x19 * x64 - 4.0_rp * x20 * x6 * x88 - x20 &
&* x82 + 90.0_rp * x20 * x83 - 90.0_rp * x20 * x86 + 22.0_rp * x6 * x88 + x6 * x91 + x82 - 90.0_rp * x83 + 195.0_rp * x86 - 6.0_rp * &
&x94 - 72.0_rp * x96 - 162.0_rp * x97 + 72.0_rp * x98) + 60.0_rp * x113 * x15 * x6 / x3**5 + 30.0_rp * x118 * x36 + 30.0_rp * x120 * &
&(7.0_rp * x5 - 3.0_rp) + x123 * (-30.0_rp * x5 + 35.0_rp * x8 + 3.0_rp) + 18.0_rp * x129 * x83 + x132 * x33 * x38 + x135 * x33 * x38&
& - 2.0_rp * x160 * x25 + x162 * (-x0 * x20 * x61 + 36.0_rp * x0 * x39 + x0 * x61 - x0 * x71 - 6.0_rp * x20 * x25 * x41 + x20 * x5&
&9 * x6 + x20 * x66 * x8 - x59 * x6 - x6 * x63 + x6 * x74 + x6 * x76 - x61 * x64 - 4.0_rp * x66 * x8 + 18.0_rp * x68 + x77) + x23&
& * (-30.0_rp * x5 + 35.0_rp * x8 + 3.0_rp) + 0.5_rp * x57 * (-x0 * x20 * x61 + 36.0_rp * x0 * x39 + x0 * x61 - x0 * x71 - 6.0_rp * x&
&20 * x25 * x41 + x20 * x59 * x6 + x20 * x66 * x8 - x59 * x6 - x6 * x63 + x6 * x74 + x6 * x76 - x61 * x64 - 4.0_rp * x66 * x8 + &
&18.0_rp * x68 + x77))
            T(i, 4) = xi * (x1 * x172 + x1 * x192 - x1 * x198 - x1 * x199 + x1 * x201 + x1 * x220 + x1 * x223 + x1 * x226 + x1&
& * x228 + x1 * x231 + x1 * x233 + x1 * x235 + x1 * x241 + x110 * x193 * x55 * x66 - x112 * (-x1 * x151 * x41 - x1 * x174 - x1 &
&* x175 + x1 * x178 + x1 * x179 + x1 * x181 + x1 * x182 - x1 * x184 - x1 * x186 * x20 + x1 * x186 - x150 + x151 + x168 + x173 *&
& x20 - x173 + x41 * x49 - x94 - 12.0_rp * x96 - 27.0_rp * x97 + 12.0_rp * x98) - 9.0_rp * x118 - 15.0_rp * x120 - x134 * x221 - x15&
&9 * x196 * x41 - x160 * x193 - x162 * x55 + x205 * x36 + x217 * x218 * x38 - x218 * x249 - x23 * x36 + x236 * x238 - x251 * x3&
&8 - x254 * x55 - 4.5_rp * x33 * x35)
            T(i, 5) = x266 * (-ai * x159 * x24 + x0 * x238 + x110 * x262 * x55 - x112 * (-x101 + x102 - x151 * x41 - x174 - x1&
&75 + x178 + x179 + x181 + x182 - x184 - x186 * x20 + x186) + x125 * x204 * x36 - x159 * x161 + x172 + x192 - x198 - x199 + x20&
&1 + x218 * x259 * x38 - x218 * x261 + x220 + x223 + x226 + x228 + x231 + x233 + x235 + x241)
            T(i, 6) = xi * (x0 * x2 * x238 + x110 * x272 * x55 - x112 * (x102 * x2 - x150 - x151 * x2 * x41 + x151 + x168 - x1&
&74 * x2 - x175 * x2 + x178 * x2 + x179 * x2 + x181 * x2 + x182 * x2 - x184 * x2 - x186 * x2 * x20 + x186 * x2 - x268 + x41 * x&
&49 - x94 - 12.0_rp * x96 - 27.0_rp * x97 + 12.0_rp * x98) - 9.0_rp * x118 - 15.0_rp * x120 - x134 * x221 - x159 * x273 * x41 - x160&
& * x271 - x162 * x55 + x172 * x2 + x192 * x2 - x198 * x2 - x199 * x2 + x2 * x201 + x2 * x220 + x2 * x223 + x2 * x226 + x2 * x2&
&28 + x2 * x231 + x2 * x233 + x2 * x235 + x2 * x241 + x218 * x284 * x38 - x218 * x288 - x23 * x36 - x251 * x38 - x254 * x55 + x&
&275 * x36 - 4.5_rp * x33 * x35)
            T(i, 7) = yi * (ai**4 * x1 * x110 * x32 * x7 + 315.0_rp * x0 * x1 * x170 + x0 * x173 * x234 - x0 * x230 + x0 * x327&
& + 7.5_rp * x1 * x129 * x19 - 4.5_rp * x1 * x134 * x141 - 4.5_rp * x1 * x141 * x35 + x1 * x191 * x38 + 7.5_rp * x1 * x225 + x1 * x&
&241 + 9.0_rp * x113 * x203 * x83 - 54.0_rp * x118 - 90.0_rp * x120 + 4.5_rp * x126 * x203 * x38 - 4.5_rp * x134 * x33 + x140 * x218&
& + x161 * x332 - x193 * x332 * x66 + x213 * x33 * x66 - x218 * x249 - x221 * x262 + x236 * x343 - x249 * x302 - x251 * x38 + x&
&254 * x323 * x38 - x301 * (51.0_rp * x1 * x117 + 11.0_rp * x1 * x180 - 16.0_rp * x1 * x183 + 5.0_rp * x1 * x19 * x20 - 5.0_rp * x1 &
&* x19 - 35.0_rp * x20 * x236 * x80 + 35.0_rp * x236 * x80 - 5.0_rp * x292 + x298 + x299) + x321 * x322 - 4.5_rp * x33 * x35 - x331&
& * x38)
            T(i, 8) = zi * (ai**4 * x1 * x110 * x32 * x7 + 315.0_rp * x0 * x1 * x170 + 7.5_rp * x1 * x129 * x19 - 4.5_rp * x1 * x&
&134 * x141 - 4.5_rp * x1 * x141 * x35 + x1 * x191 * x38 + 7.5_rp * x1 * x225 + 51.0_rp * x1 * x234 * x83 + x1 * x241 - x110 * x26&
&2 * x32 - x112 * (153.0_rp * x1 * x117 + 33.0_rp * x1 * x180 - 48.0_rp * x1 * x183 - x1 * x186 * x20 + x1 * x186 - 15.0_rp * x292 &
&+ x299 + x345) - x113 * x193 * x261 + 3.0_rp * x113 * x203 * x83 + x113 * x217 * x48 - 18.0_rp * x118 - 30.0_rp * x120 - x123 * x&
&38 - x134 * x221 - x14 * x196 * x261 * x41 + x140 * x162 + x140 * x254 + x143 * x236 * x259 - x162 * x249 - x193 * x332 * x66 &
&+ x205 * x38 + x213 * x234 * x38 + x213 * x33 * x66 + x217 * x254 * x38 - x221 * x35 - x23 * x38 - x234 * x48 + x236 * x343 - &
&x249 * x254 + x322 * x347)
            T(i, 9) = yi * (ai**4 * x110 * x2 * x32 * x7 + x0 * x143 * x2 * x259 + 315.0_rp * x0 * x170 * x2 + x0 * x2 * x343 -&
& x110 * x262 * x32 - x112 * (153.0_rp * x117 * x2 + 33.0_rp * x180 * x2 - 48.0_rp * x183 * x2 - x186 * x2 * x20 + x186 * x2 + x29&
&8 + x299 - 15.0_rp * x357 - x361 + x363) - x113 * x261 * x271 + 3.0_rp * x113 * x274 * x83 + x113 * x284 * x48 - 18.0_rp * x118 -&
& 30.0_rp * x120 - x123 * x38 + x129 * x19 * x382 - 4.5_rp * x134 * x141 * x2 - x134 * x221 - x14 * x261 * x273 * x41 + x140 * x1&
&62 + x140 * x254 - 4.5_rp * x141 * x2 * x35 - x162 * x288 + 51.0_rp * x2 * x234 * x83 + x2 * x241 - x221 * x35 + x225 * x382 - x&
&23 * x38 + x234 * x280 * x38 - x234 * x48 + x254 * x284 * x38 - x254 * x288 - x272 * x332 + x275 * x38 + x280 * x33 * x66 + x3&
&22 * x375 + x38 * x380)
            T(i, 10) = zi * (ai**4 * x110 * x2 * x32 * x7 + 315.0_rp * x0 * x170 * x2 + x0 * x2 * x343 - x0 * x230 + x0 * x234 &
&* x268 + x0 * x391 + 9.0_rp * x113 * x274 * x83 - 54.0_rp * x118 - 90.0_rp * x120 + x129 * x19 * x382 - 4.5_rp * x134 * x141 * x2 &
&- 4.5_rp * x134 * x33 + x140 * x218 - 4.5_rp * x141 * x2 * x35 + x161 * x332 + x2 * x241 - x218 * x288 - x221 * x262 + x225 * x3&
&82 - x251 * x38 + x254 * x38 * x389 - x272 * x332 + x274 * x38 * x392 + x280 * x33 * x66 - x288 * x302 - x301 * (-35.0_rp * x0 &
&* x2 * x20 * x80 + 35.0_rp * x0 * x2 * x80 + 51.0_rp * x117 * x2 + 11.0_rp * x180 * x2 - 16.0_rp * x183 * x2 + 5.0_rp * x19 * x2 * &
&x20 - 5.0_rp * x19 * x2 + x298 + x299 - 5.0_rp * x357) + x322 * x388 - 4.5_rp * x33 * x35 - x331 * x38 + x38 * x380)
            T(i, 11) = xi * (-315.0_rp * x1 * x119 + x1 * x143 * x323 - 135.0_rp * x1 * x232 + x1 * x327 + x112 * (-x1 * x20 * x&
&401 + 108.0_rp * x1 * x72 + x100 * x395 * x88 + x102 * x407 + 9.0_rp * x108 + x131 * x20 - x131 - 135.0_rp * x19 * x407 - 36.0_rp &
&* x20 * x317 + 105.0_rp * x20 * x395 * x85 + x20 * x395 * x91 - x20 * x398 + x20 * x405 - 9.0_rp * x24 * x41 + 180.0_rp * x292 + &
&24.0_rp * x314 - 6.0_rp * x316 - x395 * x406 - x395 * x409 - 5.0_rp * x395 * x91 + x398 - x405 + x41 * x51 - x70) + x113 * x173 *&
& x203 + 30.0_rp * x19 * x234 * x395 - x203 * x229 - x217 * x417 - x234 * x312 + x251 + x254 * x430 + 2.0_rp * x321 * x419 + x331&
& + x395 * x396 + x395 * x424)
            T(i, 12) = x266 * (52.5_rp * x1 * x19 * x234 + x1 * x259 * x392 + x1 * x396 + x1 * x424 + x112 * (105.0_rp * x1 * x2&
&0 * x85 + x1 * x20 * x91 - x1 * x406 - x1 * x409 - x1 * x436 + x1 * x443 - 5.0_rp * x1 * x91 + x101 - x102 + x173 * x20 * x41 -&
& x20 * x41 * x61 + x20 * x434 - 3.0_rp * x20 * x93 + x401 - x434 + 54.0_rp * x72 - 45.0_rp * x75 + 12.0_rp * x93) + 22.5_rp * x113 &
&* x19 * x203 - 157.5_rp * x119 - 13.5_rp * x125 * x126 + x204 * x323 - x218 * x259 + x218 * x347 - 67.5_rp * x232 + x321 * x57 + &
&x327)
            T(i, 13) = xi * (7.5_rp * x1 * x113 * x19 * x274 - x1 * x191 + x1 * x2 * x396 - x1 * x240 + x1 * x284 * x392 + x112&
& * (105.0_rp * x1 * x2 * x20 * x85 + x1 * x2 * x20 * x91 - x1 * x2 * x406 - x1 * x2 * x409 - x1 * x2 * x436 + x1 * x2 * x443 - &
&5.0_rp * x1 * x2 * x91 + x1 * x63 - x1 * x76 + x173 * x2 * x20 * x41 + x2 * x20 * x434 - x2 * x434 + x2 * x63 - x2 * x76 + x259&
& + x303 - x304 + x313 + 4.0_rp * x314 - x316 - x317 * x73 + x361 - x363 + x366 * x41 + 4.0_rp * x368 - x370 - x371 * x73) + 60.0&
&_rp * x113 * x125 * x459 + x113 * x19 * x203 * x382 + x123 + x143 * x259 * x458 - x2 * x240 - x205 - x213 * x234 - x217 * x254 &
&+ x217 * x464 + x23 - x234 * x280 - x254 * x284 + x254 * (-x1 * x2 * x63 + x1 * x2 * x76 - x1 * x20 * x45 + x1 * x40 + x1 * x4&
&3 + x125 - x2 * x257 - x2 * x303 - x2 * x313 - 4.0_rp * x2 * x314 + x2 * x316 + x2 * x317 * x73 + x2 * x40 + x2 * x43 - x2 * x4&
&6 - x210 + x213 - x214 - x278 + x280 + x467) - x275 + x347 * x465 + x375 * x419 - x380 + x424 * x458)
            T(i, 14) = x266 * (x112 * (x101 + x102 * x2 * x41 - x102 + x2 * x20 * x81 + 105.0_rp * x2 * x20 * x85 + x2 * x20 * &
&x91 - x2 * x406 - x2 * x409 - x2 * x436 + x2 * x443 - x2 * x81 - 5.0_rp * x2 * x91 - x20 * x41 * x61 - 3.0_rp * x20 * x93 + x401&
& + 54.0_rp * x72 - 45.0_rp * x75 + 12.0_rp * x93) + 22.5_rp * x113 * x19 * x274 - 157.5_rp * x119 - 13.5_rp * x125 * x126 + 52.5_rp *&
& x19 * x2 * x234 + x2 * x396 + x2 * x424 + x204 * x389 - x218 * x259 + x218 * x375 - 67.5_rp * x232 + x259 * x464 + x388 * x57 &
&+ x391)
            T(i, 15) = xi * (x112 * (x102 * x41 * x475 + 9.0_rp * x108 + x131 * x20 - x131 - 90.0_rp * x19 * x2 * x20 + 90.0_rp *&
& x19 * x2 - x2 * x20 * x401 + 108.0_rp * x2 * x72 - 36.0_rp * x20 * x371 + x20 * x475 * x81 + 105.0_rp * x20 * x475 * x85 + x20 *&
& x475 * x91 - 9.0_rp * x24 * x41 + 180.0_rp * x357 + 24.0_rp * x368 - 6.0_rp * x370 - x406 * x475 - x409 * x475 + x41 * x51 - x436&
& * x475 + x443 * x475 - x475 * x81 - 5.0_rp * x475 * x91 - x70) + x113 * x268 * x274 - 315.0_rp * x119 * x2 + x143 * x2 * x389 +&
& 30.0_rp * x19 * x234 * x475 - 135.0_rp * x2 * x232 + x2 * x391 - x229 * x274 - x234 * x366 + x251 + x254 * (36.0_rp * x2 * x39 -&
& x2 * x71 - x20 * x366 + x20 * x475 * x59 + x20 * x475 * x93 + 18.0_rp * x271 * x41 - 6.0_rp * x278 + x366 - x41 * x475 * x61 - &
&x475 * x59 - x475 * x63 + x475 * x74 + x475 * x76 - 4.0_rp * x475 * x93 + x77) - x284 * x417 + x331 + 2.0_rp * x388 * x465 + x39&
&6 * x475 + x424 * x475)
            T(i, 16) = yi * (x112 * (-150.0_rp * x1 * x19 * x20 + 150.0_rp * x1 * x19 + x100 * x395 * x88 + x102 * x407 + x109 -&
& 135.0_rp * x19 * x407 + 105.0_rp * x20 * x395 * x85 + x20 * x395 * x91 + x20 * x405 + 300.0_rp * x292 + 40.0_rp * x314 - 10.0_rp *&
& x316 - x317 * x99 + 180.0_rp * x317 - x395 * x406 - x395 * x409 - 5.0_rp * x395 * x91 - x405 - 150.0_rp * x491) + 15.0_rp * x126 &
&* x203 * x485 + x18 * (-30.0_rp * x196 + 35.0_rp * x395 * x7 + 3.0_rp) + x23 * (-70.0_rp * x196 + 63.0_rp * x395 * x7 + 15.0_rp) + x&
&323 * x487 * x58 + 2.5_rp * x430 * x57)
            T(i, 17) = zi * (18.0_rp * x1 * x113 * x19 * x203 + 6.0_rp * x1 * x125 * x126 * x485 + 2.0_rp * x1 * x134 * x14 * x32&
&3 + x1 * x143 * x323 - x112 * (x168 - 30.0_rp * x19 * x20 * x407 + 105.0_rp * x19 * x407 + x20 * x313 - 90.0_rp * x20 * x395 * x8&
&5 - 4.0_rp * x20 * x395 * x88 + x20 * x398 - x20 * x405 - 162.0_rp * x292 - 6.0_rp * x314 - 72.0_rp * x317 + 195.0_rp * x395 * x85 &
&+ 22.0_rp * x395 * x88 + x395 * x91 - x398 + x405 + 72.0_rp * x491) + 60.0_rp * x113 * x15 * x395 / x3**5 + x123 * (-30.0_rp * x19&
&6 + 35.0_rp * x395 * x7 + 3.0_rp) + x135 * x14 * x203 * x487 + x162 * x430 + x203 * x229 * x487 + x23 * (-30.0_rp * x196 + 35.0_rp&
& * x395 * x7 + 3.0_rp) + x254 * x430 - x417 * x487 * x496 - 2.0_rp * x419 * (x131 * x20 - x131 + x20 * x40 - x214 * x41 + 27.0_rp&
& * x292 + x303 - x304 + x314 + 12.0_rp * x317 - 15.0_rp * x39 - x43 - 12.0_rp * x491) + 30.0_rp * x485 * x503 + 30.0_rp * x492 * (7&
&.0_rp * x196 - 3.0_rp))
            T(i, 18) = yi * (10.5_rp * x10 * x14 * x2 * x203 * x41 + x110 * x272 * x323 - x112 * (-72.0_rp * x1 * x2 * x20 * x85&
& + x102 * x2 + x166 - x174 * x2 - x175 * x2 + x179 * x2 - x2 * x20 * x434 + x2 * x434 - x20 * x40 + x214 * x41 - x268 - 27.0_rp&
& * x292 - x303 + x304 - x314 - 12.0_rp * x317 + 15.0_rp * x39 + 72.0_rp * x41 * x459 - x41 * x467 + x43 + 177.0_rp * x458 * x85 + &
&7.0_rp * x458 * x88 + 12.0_rp * x491) + 22.5_rp * x113 * x19 * x2 * x203 + x122 * x134 * x2 * x323 - x122 * x134 * x203 + x122 * &
&x2 * x203 * x34 * x66 - 4.5_rp * x126 * x203 - x135 * x14 * x2 * x496 - x14 * x273 * x41 * (x131 * x20 - x131 + x20 * x40 - x21&
&4 * x41 + 27.0_rp * x292 + x303 - x304 + x314 + 12.0_rp * x317 - 15.0_rp * x39 - x43 - 12.0_rp * x491) - x162 * x323 + 210.0_rp * x&
&170 * x458 + x2 * x204 * x323 - x2 * x229 * x496 + x2 * x230 * x487 + 15.0_rp * x2 * x232 * x485 + x2 * x240 * x487 + x218 * x2&
&84 * x487 - x218 * (2.0_rp * x1 * x108 - 5.0_rp * x1 * x39 + x2 * x257 + 23.0_rp * x2 * x292 + x2 * x303 + 5.0_rp * x2 * x317 - x2&
& * x45 - 8.0_rp * x2 * x491 + x206 - x213 + x214 + x27 - x280 - x428 - x467) - x23 * x485 + 18.0_rp * x234 * x459 + x238 * x458 &
&- x251 * x487 - x254 * x323 + x275 * x485 + x380 * x485 - x465 * (x131 * x20 - x131 + x20 * x40 - x214 * x41 + 27.0_rp * x292 +&
& x303 - x304 + x314 + 12.0_rp * x317 - 15.0_rp * x39 - x43 - 12.0_rp * x491) - 15.0_rp * x492 - 9.0_rp * x503)
            T(i, 19) = zi * (ai**4 * x110 * x2 * x203 * x7 + 9.0_rp * x1 * x113 * x19 * x274 - x1 * x230 + x1 * x391 + x10 * x1&
&4 * x203 * x382 * x41 + x113 * x19 * x203 * x382 - x122 * x203 * x262 - x122 * x272 * x496 - 4.5_rp * x126 * x203 - 4.5_rp * x13&
&4 * x14 * x2 * x496 - 4.5_rp * x134 * x14 * x203 + x14 * x203 * x280 * x66 + 315.0_rp * x170 * x458 + x173 * x2 * x234 + x2 * x2&
&40 * x487 + x218 * x496 - x218 * (2.0_rp * x1 * x108 - 5.0_rp * x1 * x39 + x2 * x257 + 23.0_rp * x2 * x292 + x2 * x303 + 5.0_rp * &
&x2 * x317 - x2 * x45 - 8.0_rp * x2 * x491 + x206 - x213 + x214 + x27 - x280 - x428 - x467) - x251 * x487 + x254 * x389 * x487 +&
& x274 * x392 * x487 - x301 * (-35.0_rp * x1 * x2 * x20 * x80 - 16.0_rp * x1 * x2 * x20 * x85 + 35.0_rp * x1 * x2 * x80 + 5.0_rp * &
&x19 * x2 * x20 - 5.0_rp * x19 * x2 - 23.0_rp * x292 - 5.0_rp * x317 + x345 - 5.0_rp * x357 + 11.0_rp * x41 * x459 + x45 + 51.0_rp * &
&x458 * x85 + 8.0_rp * x491) + x302 * x496 - x302 * (2.0_rp * x1 * x108 - 5.0_rp * x1 * x39 + x2 * x257 + 23.0_rp * x2 * x292 + x2 &
&* x303 + 5.0_rp * x2 * x317 - x2 * x45 - 8.0_rp * x2 * x491 + x206 - x213 + x214 + x27 - x280 - x428 - x467) - x331 * x487 + x34&
&3 * x458 + x380 * x487 + x388 * x419 - x464 * x496 - 90.0_rp * x492 - 54.0_rp * x503)
            T(i, 20) = yi * (x112 * (x102 * x41 * x475 + 9.0_rp * x108 + x131 * x20 - x131 - 90.0_rp * x19 * x2 * x20 + 90.0_rp *&
& x19 * x2 - x2 * x20 * x401 + 108.0_rp * x2 * x72 - 36.0_rp * x20 * x371 + x20 * x475 * x81 + 105.0_rp * x20 * x475 * x85 + x20 *&
& x475 * x91 - 9.0_rp * x24 * x41 + 180.0_rp * x357 + 24.0_rp * x368 - 6.0_rp * x370 - x406 * x475 - x409 * x475 + x41 * x51 - x436&
& * x475 + x443 * x475 - x475 * x81 - 5.0_rp * x475 * x91 - x70) + x113 * x268 * x274 - 315.0_rp * x119 * x2 + x143 * x2 * x389 +&
& 30.0_rp * x19 * x234 * x475 - 135.0_rp * x2 * x232 + x2 * x391 - x229 * x274 - x234 * x366 + x251 + x254 * (36.0_rp * x2 * x39 -&
& x2 * x71 - x20 * x366 + x20 * x475 * x59 + x20 * x475 * x93 + 18.0_rp * x271 * x41 - 6.0_rp * x278 + x366 - x41 * x475 * x61 - &
&x475 * x59 - x475 * x63 + x475 * x74 + x475 * x76 - 4.0_rp * x475 * x93 + x77) - x284 * x417 + x331 + 2.0_rp * x388 * x465 + x39&
&6 * x475 + x424 * x475)
            T(i, 21) = zi * (x112 * (x102 * x41 * x475 + x109 - 150.0_rp * x19 * x2 * x20 + 150.0_rp * x19 * x2 - 150.0_rp * x2 *&
& x75 + x20 * x475 * x81 + 105.0_rp * x20 * x475 * x85 + x20 * x475 * x91 + 300.0_rp * x357 + 40.0_rp * x368 - 10.0_rp * x370 - x37&
&1 * x99 + 180.0_rp * x371 - x406 * x475 - x409 * x475 - x436 * x475 + x443 * x475 - x475 * x81 - 5.0_rp * x475 * x91) + 15.0_rp *&
& x126 * x274 * (5.0_rp * x273 - 3.0_rp) + x18 * (-30.0_rp * x273 + 35.0_rp * x475 * x7 + 3.0_rp) + x23 * (-70.0_rp * x273 + 63.0_rp *&
& x475 * x7 + 15.0_rp) + x389 * x58 * (3.0_rp * x273 - 1.0_rp) + 2.5_rp * x57 * (36.0_rp * x2 * x39 - x2 * x71 - x20 * x366 + x20 * &
&x475 * x59 + x20 * x475 * x93 + 18.0_rp * x271 * x41 - 6.0_rp * x278 + x366 - x41 * x475 * x61 - x475 * x59 - x475 * x63 + x475 &
&* x74 + x475 * x76 - 4.0_rp * x475 * x93 + x77))
        end do
    end subroutine T5_damp_thole
    pure subroutine T6_damp_thole(x, y, z, a, T)
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
        real(rp) :: x14
        real(rp) :: x15
        real(rp) :: x17
        real(rp) :: x18
        real(rp) :: x20
        real(rp) :: x21
        real(rp) :: x22
        real(rp) :: x23
        real(rp) :: x24
        real(rp) :: x25
        real(rp) :: x27
        real(rp) :: x28
        real(rp) :: x29
        real(rp) :: x30
        real(rp) :: x31
        real(rp) :: x32
        real(rp) :: x34
        real(rp) :: x35
        real(rp) :: x36
        real(rp) :: x38
        real(rp) :: x39
        real(rp) :: x40
        real(rp) :: x42
        real(rp) :: x43
        real(rp) :: x44
        real(rp) :: x45
        real(rp) :: x46
        real(rp) :: x47
        real(rp) :: x49
        real(rp) :: x50
        real(rp) :: x51
        real(rp) :: x52
        real(rp) :: x53
        real(rp) :: x55
        real(rp) :: x58
        real(rp) :: x59
        real(rp) :: x61
        real(rp) :: x63
        real(rp) :: x64
        real(rp) :: x66
        real(rp) :: x67
        real(rp) :: x68
        real(rp) :: x69
        real(rp) :: x70
        real(rp) :: x72
        real(rp) :: x74
        real(rp) :: x76
        real(rp) :: x77
        real(rp) :: x78
        real(rp) :: x79
        real(rp) :: x80
        real(rp) :: x81
        real(rp) :: x83
        real(rp) :: x84
        real(rp) :: x86
        real(rp) :: x87
        real(rp) :: x88
        real(rp) :: x90
        real(rp) :: x91
        real(rp) :: x92
        real(rp) :: x93
        real(rp) :: x94
        real(rp) :: x95
        real(rp) :: x96
        real(rp) :: x97
        real(rp) :: x98
        real(rp) :: x99
        real(rp) :: x100
        real(rp) :: x102
        real(rp) :: x103
        real(rp) :: x104
        real(rp) :: x107
        real(rp) :: x108
        real(rp) :: x109
        real(rp) :: x110
        real(rp) :: x111
        real(rp) :: x112
        real(rp) :: x115
        real(rp) :: x119
        real(rp) :: x120
        real(rp) :: x121
        real(rp) :: x123
        real(rp) :: x125
        real(rp) :: x127
        real(rp) :: x128
        real(rp) :: x129
        real(rp) :: x130
        real(rp) :: x134
        real(rp) :: x135
        real(rp) :: x136
        real(rp) :: x137
        real(rp) :: x138
        real(rp) :: x140
        real(rp) :: x141
        real(rp) :: x143
        real(rp) :: x145
        real(rp) :: x146
        real(rp) :: x147
        real(rp) :: x149
        real(rp) :: x150
        real(rp) :: x151
        real(rp) :: x153
        real(rp) :: x154
        real(rp) :: x155
        real(rp) :: x157
        real(rp) :: x158
        real(rp) :: x159
        real(rp) :: x160
        real(rp) :: x161
        real(rp) :: x162
        real(rp) :: x163
        real(rp) :: x165
        real(rp) :: x167
        real(rp) :: x168
        real(rp) :: x169
        real(rp) :: x170
        real(rp) :: x172
        real(rp) :: x174
        real(rp) :: x175
        real(rp) :: x177
        real(rp) :: x178
        real(rp) :: x179
        real(rp) :: x180
        real(rp) :: x185
        real(rp) :: x186
        real(rp) :: x187
        real(rp) :: x188
        real(rp) :: x189
        real(rp) :: x191
        real(rp) :: x192
        real(rp) :: x196
        real(rp) :: x199
        real(rp) :: x200
        real(rp) :: x205
        real(rp) :: x208
        real(rp) :: x209
        real(rp) :: x210
        real(rp) :: x214
        real(rp) :: x219
        real(rp) :: x223
        real(rp) :: x228
        real(rp) :: x230
        real(rp) :: x231
        real(rp) :: x233
        real(rp) :: x234
        real(rp) :: x235
        real(rp) :: x237
        real(rp) :: x239
        real(rp) :: x240
        real(rp) :: x241
        real(rp) :: x243
        real(rp) :: x244
        real(rp) :: x245
        real(rp) :: x246
        real(rp) :: x247
        real(rp) :: x248
        real(rp) :: x250
        real(rp) :: x252
        real(rp) :: x253
        real(rp) :: x256
        real(rp) :: x261
        real(rp) :: x263
        real(rp) :: x268
        real(rp) :: x270
        real(rp) :: x282
        real(rp) :: x283
        real(rp) :: x287
        real(rp) :: x288
        real(rp) :: x291
        real(rp) :: x292
        real(rp) :: x293
        real(rp) :: x295
        real(rp) :: x297
        real(rp) :: x299
        real(rp) :: x302
        real(rp) :: x303
        real(rp) :: x304
        real(rp) :: x305
        real(rp) :: x307
        real(rp) :: x308
        real(rp) :: x313
        real(rp) :: x314
        real(rp) :: x317
        real(rp) :: x318
        real(rp) :: x319
        real(rp) :: x320
        real(rp) :: x321
        real(rp) :: x322
        real(rp) :: x326
        real(rp) :: x327
        real(rp) :: x329
        real(rp) :: x331
        real(rp) :: x332
        real(rp) :: x333
        real(rp) :: x335
        real(rp) :: x336
        real(rp) :: x338
        real(rp) :: x340
        real(rp) :: x341
        real(rp) :: x342
        real(rp) :: x344
        real(rp) :: x345
        real(rp) :: x347
        real(rp) :: x348
        real(rp) :: x349
        real(rp) :: x350
        real(rp) :: x351
        real(rp) :: x353
        real(rp) :: x359
        real(rp) :: x362
        real(rp) :: x363
        real(rp) :: x364
        real(rp) :: x366
        real(rp) :: x367
        real(rp) :: x368
        real(rp) :: x369
        real(rp) :: x371
        real(rp) :: x372
        real(rp) :: x373
        real(rp) :: x374
        real(rp) :: x375
        real(rp) :: x376
        real(rp) :: x377
        real(rp) :: x378
        real(rp) :: x379
        real(rp) :: x380
        real(rp) :: x381
        real(rp) :: x382
        real(rp) :: x383
        real(rp) :: x384
        real(rp) :: x385
        real(rp) :: x389
        real(rp) :: x391
        real(rp) :: x393
        real(rp) :: x394
        real(rp) :: x395
        real(rp) :: x397
        real(rp) :: x400
        real(rp) :: x401
        real(rp) :: x404
        real(rp) :: x405
        real(rp) :: x406
        real(rp) :: x407
        real(rp) :: x408
        real(rp) :: x410
        real(rp) :: x411
        real(rp) :: x413
        real(rp) :: x414
        real(rp) :: x416
        real(rp) :: x421
        real(rp) :: x422
        real(rp) :: x426
        real(rp) :: x427
        real(rp) :: x428
        real(rp) :: x429
        real(rp) :: x430
        real(rp) :: x432
        real(rp) :: x435
        real(rp) :: x436
        real(rp) :: x437
        real(rp) :: x440
        real(rp) :: x443
        real(rp) :: x444
        real(rp) :: x445
        real(rp) :: x446
        real(rp) :: x450
        real(rp) :: x453
        real(rp) :: x455
        real(rp) :: x456
        real(rp) :: x457
        real(rp) :: x459
        real(rp) :: x460
        real(rp) :: x462
        real(rp) :: x470
        real(rp) :: x471
        real(rp) :: x473
        real(rp) :: x474
        real(rp) :: x475
        real(rp) :: x477
        real(rp) :: x478
        real(rp) :: x480
        real(rp) :: x481
        real(rp) :: x482
        real(rp) :: x483
        real(rp) :: x485
        real(rp) :: x488
        real(rp) :: x489
        real(rp) :: x490
        real(rp) :: x494
        real(rp) :: x495
        real(rp) :: x496
        real(rp) :: x499
        real(rp) :: x500
        real(rp) :: x502
        real(rp) :: x505
        real(rp) :: x506
        real(rp) :: x508
        real(rp) :: x509
        real(rp) :: x510
        real(rp) :: x511
        real(rp) :: x514
        real(rp) :: x515
        real(rp) :: x516
        real(rp) :: x518
        real(rp) :: x522
        real(rp) :: x526
        real(rp) :: x527
        real(rp) :: x528
        real(rp) :: x529
        real(rp) :: x536
        real(rp) :: x537
        real(rp) :: x538
        real(rp) :: x539
        real(rp) :: x540
        real(rp) :: x542
        real(rp) :: x543
        real(rp) :: x545
        real(rp) :: x547
        real(rp) :: x548
        real(rp) :: x551
        real(rp) :: x552
        real(rp) :: x557
        real(rp) :: x558
        real(rp) :: x560
        real(rp) :: x565
        real(rp) :: x566
        real(rp) :: x567
        real(rp) :: x569
        real(rp) :: x571
        real(rp) :: x575
        real(rp) :: x582
        real(rp) :: x583
        real(rp) :: x585
        real(rp) :: x589
        real(rp) :: x591
        real(rp) :: x593
        real(rp) :: x595
        real(rp) :: x596
        real(rp) :: x598
        real(rp) :: x603
        real(rp) :: x604
        real(rp) :: x605
        real(rp) :: x611
        real(rp) :: x622
        real(rp) :: x623
        real(rp) :: x624
        real(rp) :: x628
        real(rp) :: x629
        real(rp) :: x630
        real(rp) :: x632
        real(rp) :: x633
        real(rp) :: x635
        real(rp) :: x639
        real(rp) :: x641
        real(rp) :: x643
        real(rp) :: x649
        real(rp) :: x650
        real(rp) :: x651
        real(rp) :: x652
        real(rp) :: x654
        real(rp) :: x657
        real(rp) :: x663
        real(rp) :: x664
        real(rp) :: x665
        real(rp) :: x667
        real(rp) :: x671
        real(rp) :: x674
        real(rp) :: x680
        real(rp) :: x685
        real(rp) :: x687
        real(rp) :: x688
        real(rp) :: x697
        real(rp) :: x701
        real(rp) :: x706
        real(rp) :: x707
        real(rp) :: x708
        real(rp) :: x709
        real(rp) :: x710
        real(rp) :: x712
        real(rp) :: x717
        real(rp) :: x719
        real(rp) :: x720
        real(rp) :: x725
        real(rp) :: x727
        real(rp) :: x733
        real(rp) :: x736
        real(rp) :: x738
        real(rp) :: x739
        real(rp) :: x740
        real(rp) :: x742
        real(rp) :: x744
        real(rp) :: x745
        real(rp) :: x747
        real(rp) :: x748
        real(rp) :: x749
        real(rp) :: x753
        real(rp) :: x754
        real(rp) :: x757
        real(rp) :: x758
        real(rp) :: x764
        real(rp) :: x766
        real(rp) :: x769
        real(rp) :: x771
        real(rp) :: x772
        real(rp) :: x776
        real(rp) :: x778
        real(rp) :: x782
        real(rp) :: x786
        real(rp) :: x789
        real(rp) :: x790
        real(rp) :: x793
        real(rp) :: x794
        real(rp) :: x797
        real(rp) :: x802
        real(rp) :: x804
        real(rp) :: x807
        real(rp) :: x809
        real(rp) :: x810
        real(rp) :: x811
        real(rp) :: x812
        real(rp) :: x813
        real(rp) :: x814
        real(rp) :: x815
        real(rp) :: x816
        real(rp) :: x817
        real(rp) :: x818
        real(rp) :: x819
        real(rp) :: x820
        real(rp) :: x821
        real(rp) :: x827
        real(rp) :: x830
        real(rp) :: x831
        real(rp) :: x835
        real(rp) :: x838
        real(rp) :: x839
        real(rp) :: x843
        real(rp) :: x849
        real(rp) :: x850
        real(rp) :: x856
        real(rp) :: x859
        real(rp) :: x863
        real(rp) :: x865
        real(rp) :: x867
        real(rp) :: x871
        real(rp) :: x872
        real(rp) :: x875
        real(rp) :: x876
        real(rp) :: x893
        real(rp) :: x895
        real(rp) :: x897
        real(rp) :: x898
        real(rp) :: x902
        real(rp) :: x904
        real(rp) :: x907
        real(rp) :: x909
        real(rp) :: x911
        real(rp) :: x913
        real(rp) :: x918
        real(rp) :: x922
        real(rp) :: x923
        real(rp) :: x929
        real(rp) :: x935
        real(rp) :: x941
        real(rp) :: x942
        real(rp) :: x944
        real(rp) :: x945
        real(rp) :: x946
        real(rp) :: x947
        real(rp) :: x958
        real(rp) :: x961
        real(rp) :: x963
        real(rp) :: x965
        real(rp) :: x966
        real(rp) :: x969
        real(rp) :: x975
        real(rp) :: x983
        real(rp) :: x984
        real(rp) :: x1010
        real(rp) :: x1011
        real(rp) :: x1013
        real(rp) :: x1015
        real(rp) :: x1019
        real(rp) :: x1020
        real(rp) :: x1022
        real(rp) :: x1026
        real(rp) :: x1040
        real(rp) :: x1045
        real(rp) :: x1060
        real(rp) :: x1067
        real(rp) :: x1069
        real(rp) :: x1081
        real(rp) :: x1087
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
            x6 = xi**4
            x7 = x3**(-2)
            x8 = x6 * x7
            x10 = x3**(-4)
            x11 = ai * x10
            x12 = x0 * x11
            x14 = ai * sqrt(x3)
            x15 = exp(-x14)
            x17 = x15 * (x14 + 1.0_rp)
            x18 = x12 * x17
            x20 = xi**6
            x21 = x3**(-3)
            x22 = 231.0_rp * x21
            x23 = x3**(-3.5_rp)
            x24 = x14 + 2.0_rp
            x25 = x15 * x24 - 2.0_rp
            x27 = 22.5_rp * x23 * x25
            x28 = -30.0_rp * x5 + 35.0_rp * x8 + 3.0_rp
            x29 = x3**(-2.5_rp)
            x30 = ai * x29
            x31 = x3**(-1.5_rp)
            x32 = x0 * x31
            x34 = ai * x4
            x35 = 2.0_rp * x34
            x36 = x24 * x34
            x38 = -x24 * x3**(-0.5_rp) + x3**(-0.5_rp)
            x39 = -x0 * x35 + x0 * x36 + x24 * x32 - x32 + x38
            x40 = x15 * x39
            x42 = x15 * x30
            x43 = ai * x7
            x44 = 6.0_rp * x43
            x45 = ai**2
            x46 = 3.0_rp * x31
            x47 = x45 * x46
            x49 = 3.0_rp * x43
            x50 = x24 * x49
            x51 = 3.0_rp * x29
            x52 = x0 * x51
            x53 = x24 * x52
            x55 = x24 * x46
            x58 = 6.0_rp * x34 - 3.0_rp * x36 + x46 - x55
            x59 = -x0 * x44 - x0 * x47 + x0 * x50 + x24 * x32 * x45 - x52 + x53 + x58
            x61 = 5.0_rp * x5 - 3.0_rp
            x63 = 3.0_rp * x5 - 1.0_rp
            x64 = 15.0_rp * x23
            x66 = 18.0_rp * x29
            x67 = x0 * x66
            x68 = ai * x21
            x69 = 30.0_rp * x68
            x70 = x45 * x6
            x72 = ai**3
            x74 = x32 * x45
            x76 = 18.0_rp * x43
            x77 = x0 * x24
            x78 = x29 * x45
            x79 = x24 * x78
            x80 = 6.0_rp * x79
            x81 = 15.0_rp * x68
            x83 = -6.0_rp * x34 + 3.0_rp * x36 - x46 + x55
            x84 = 36.0_rp * x0 * x43 - 6.0_rp * x24 * x32 * x45 + x24 * x6 * x64 + x24 * x6 * x81 - x24 * x67 + x24 * x72 * x8 -&
& x6 * x64 - x6 * x69 + x6 * x80 - x66 * x70 + x67 - 4.0_rp * x72 * x8 + 18.0_rp * x74 - x76 * x77 + x83
            x86 = ai * x15 * x31
            x87 = 7.5_rp * x86
            x88 = x0 * x23
            x90 = x11 * x6
            x91 = x23 * x70
            x92 = x21 * x72
            x93 = x6 * x92
            x94 = ai**4
            x95 = x29 * x94
            x96 = x6 * x95
            x97 = x7 * x72
            x98 = x0 * x97
            x99 = 180.0_rp * x78
            x100 = x0 * x68
            x102 = 150.0_rp * x24
            x103 = x0 * x78
            x104 = 60.0_rp * x24
            x107 = 45.0_rp * x23
            x108 = x107 * x24
            x109 = x24 * x90
            x110 = x3**(-4.5_rp)
            x111 = x110 * x6
            x112 = 105.0_rp * x111
            x115 = 45.0_rp * x29
            x119 = x115 * x24
            x120 = x24 * x31
            x121 = x120 * x45
            x123 = x24 * x43
            x125 = -x115 + x119 + 15.0_rp * x121 + 45.0_rp * x123 - 45.0_rp * x31 * x45 - 90.0_rp * x43
            x127 = ai * x46
            x128 = x127 * x15
            x129 = x3**(-5.5_rp)
            x130 = 945.0_rp * x129
            x134 = ai / x3**5
            x135 = 1890.0_rp * x134
            x136 = 1260.0_rp * x45
            x137 = x110 * x136
            x138 = x10 * x72
            x140 = x23 * x94
            x141 = 75.0_rp * x140
            x143 = ai**5 * x21
            x145 = 675.0_rp * x24
            x146 = x24 * x64
            x147 = x146 * x94
            x149 = 105.0_rp * x138
            x150 = 270.0_rp * x78
            x151 = x24 * x45
            x153 = 945.0_rp * x134
            x154 = x115 - x119 - 15.0_rp * x121 - 45.0_rp * x123 + 45.0_rp * x31 * x45 + 90.0_rp * x43
            x155 = 0.5_rp * x15
            x157 = ai * x155 * x3**(-0.5_rp)
            x158 = 7.0_rp * x5
            x159 = x158 - 3.0_rp
            x160 = x134 * x17
            x161 = x0 * x160
            x162 = x0 * x25
            x163 = x129 * x162
            x165 = x11 * x17
            x167 = 7.5_rp * x165
            x168 = x110 * x25
            x169 = 52.5_rp * x168
            x170 = 7.5_rp * x42
            x172 = x120 - x31 - x35 + x36
            x174 = x0 * x110
            x175 = ai * x40
            x177 = ai * x23
            x178 = x177 * x40
            x179 = x21 * x45
            x180 = x179 * x40
            x185 = -2.0_rp * x0 * x123 + 5.0_rp * x0 * x43 + x120 - x31 - x34 + x52 - x53 + x74
            x186 = 15.0_rp * x42
            x187 = ai * x15
            x188 = x187 * x88
            x189 = 30.0_rp * x188 * x59
            x191 = x45 * x7
            x192 = x15 * x191
            x196 = x100 * x24
            x199 = x0 * x64
            x200 = x199 * x24
            x205 = 9.0_rp * x29
            x208 = x205 * x24 - x205 + x24 * x44 - 15.0_rp * x43 - x47
            x209 = 27.0_rp * x100 + 12.0_rp * x103 - 12.0_rp * x196 + x199 - x200 + x208 - x45 * x53 + x98
            x210 = 5.0_rp * x86
            x214 = 90.0_rp * x88
            x219 = x214 * x24
            x223 = 30.0_rp * x24
            x228 = -x205 * x24 + x205
            x230 = x228 - x24 * x44 + 15.0_rp * x43 + x47
            x231 = x15 * (-162.0_rp * x100 - 72.0_rp * x103 - 90.0_rp * x109 - x112 * x24 + x112 + 72.0_rp * x196 - x214 + x219 - &
&x223 * x91 + x230 + x24 * x45 * x67 - 4.0_rp * x24 * x93 + 195.0_rp * x90 + 105.0_rp * x91 + 22.0_rp * x93 + x96 - 6.0_rp * x98)
            x233 = 0.5_rp * x86
            x234 = x4 * x45
            x235 = x155 * x234
            x237 = x130 * x6
            x239 = x237 * x24
            x240 = x0 * x92
            x241 = 10.0_rp * x95
            x243 = x138 * x6
            x244 = x134 * x6
            x245 = x24 * x244
            x246 = x111 * x45
            x247 = x24 * x246
            x248 = 5.0_rp * x23
            x250 = x45 * x88
            x252 = x12 * x24
            x253 = 225.0_rp * x23
            x256 = x24 * x68
            x261 = 15.0_rp * x78
            x263 = 69.0_rp * x68
            x268 = 24.0_rp * x256
            x270 = x1 * x23
            x282 = x1 * x107
            x283 = x24 * x282
            x287 = x0 * x1
            x288 = 42.0_rp * x7
            x291 = 30.0_rp * x174 * x25
            x292 = x1 * x110
            x293 = x25 * x292
            x295 = x1 * x31
            x297 = x1 * x4
            x299 = x129 * x25
            x302 = 7.0_rp * x240
            x303 = 72.0_rp * x250
            x304 = 177.0_rp * x12
            x305 = 72.0_rp * x252
            x307 = 105.0_rp * x174
            x308 = x1 * x307
            x313 = -x1 * x200 * x45 - x1 * x261 - x1 * x263 + x1 * x268 + x1 * x302 + x1 * x303 + x1 * x304 - x1 * x305 - 27.0&
&_rp * x100 - 12.0_rp * x103 + 12.0_rp * x196 - x199 + x200 + x230 - x24 * x308 - x282 + x283 + x308 + x45 * x53 - x98
            x314 = x187 * x32
            x317 = 5.0_rp * x78
            x318 = 23.0_rp * x68
            x319 = 8.0_rp * x196
            x320 = -x120 + x31
            x321 = x1 * x51
            x322 = x24 * x321
            x326 = 2.0_rp * x0 * x123 - 5.0_rp * x0 * x43 + x1 * x199 - x1 * x200 - x1 * x319 - x1 * x49 + x287 * x317 + x287 * &
&x318 + x320 - x321 + x322 + x34 - x52 + x53 - x74
            x327 = x15 * x326
            x329 = 1.5_rp * x15 * x30 * x84
            x331 = x1 * x120 - x1 * x35 + x1 * x36 - x295 + x38
            x332 = 1.5_rp * x15
            x333 = x30 * x332
            x335 = 1.5_rp * x15 * x191 * x84
            x336 = x295 * x72
            x338 = x187 * (x14 + 1.0_rp) / x3**6
            x340 = 840.0_rp * x338 * x6
            x341 = 180.0_rp * x174 * x175
            x342 = x172 * x187
            x344 = 120.0_rp * x111 * x342
            x345 = x175 * x63
            x347 = x10 * x40 * x45
            x348 = 36.0_rp * x347
            x349 = 21.0_rp * x180 * x63
            x350 = x1 * x11
            x351 = x17 * x350
            x353 = 14.0_rp * x179
            x359 = x1 * x120 * x45
            x362 = x1 * x24 * x49 - x1 * x44 - x1 * x47 + x320 - x321 + x322 + x35 + x359 - x36
            x363 = 6.0_rp * x42
            x364 = x0 * x363 * x61
            x366 = x40 * x63 * x72
            x367 = x29 * x72
            x368 = 2.0_rp * x367
            x369 = 4.0_rp * x192
            x371 = x185 * x192
            x372 = 6.0_rp * x371 * x63
            x373 = 12.0_rp * x42
            x374 = x1 * x66
            x375 = x185 * x187
            x376 = x375 * x63
            x377 = 36.0_rp * x185 * x188
            x378 = 210.0_rp * x61
            x379 = x134 * x287
            x380 = x17 * x379
            x381 = 60.0_rp * x17
            x382 = x159 * x381
            x383 = x172 * x188
            x384 = 60.0_rp * x383 * x61
            x385 = 7.5_rp * x23 * x25
            x389 = 18.0_rp * x88
            x391 = x24 * x51
            x393 = x121 + x391 - x44 - x47 + x50 - x51
            x394 = x0 * x317
            x395 = x0 * x318
            x397 = x199 - x200 - x319 + x391 + x394 + x395 - x49 - x51
            x400 = x24 * x307
            x401 = x200 * x45
            x404 = ai * x15 * (-x107 + x108 - x261 - x263 + x268 + x302 + x303 + x304 - x305 + x307 - x400 - x401)
            x405 = x31 * x72
            x406 = x155 * x405
            x407 = 210.0_rp * x174
            x408 = 315.0_rp * x129
            x410 = 8.0_rp * x256
            x411 = x24 * x407
            x413 = ai * x3**(-0.5_rp) * x332
            x414 = yi * zi
            x416 = x2 * x23
            x421 = x107 * x2
            x422 = x108 * x2
            x426 = x0 * x2
            x427 = x110 * x2
            x428 = x25 * x427
            x429 = x2 * x31
            x430 = x2 * x4
            x432 = x2 * x400
            x435 = -27.0_rp * x100 - 12.0_rp * x103 + 12.0_rp * x196 - x199 - x2 * x261 - x2 * x263 + x2 * x268 + x2 * x302 + x2 &
                  &* x303 + x2 * x304 - x2 * x305 + x2 * x307 - x2 * x401 + x200 + x230 - x421 + x422 - x432 + x45 * x53 - x98
            x436 = x2 * x51
            x437 = x2 * x391
            x440 = 2.0_rp * x0 * x123 - 5.0_rp * x0 * x43 + x199 * x2 - x2 * x200 - x2 * x319 + x2 * x394 + x2 * x395 - x2 * x49&
& + x320 + x34 - x436 + x437 - x52 + x53 - x74
            x443 = x120 * x2 - x2 * x35 + x2 * x36 + x38 - x429
            x444 = x429 * x72
            x445 = x11 * x2
            x446 = x17 * x445
            x450 = x121 * x2
            x453 = -x2 * x44 - x2 * x47 + x2 * x50 + x320 + x35 - x36 - x436 + x437 + x450
            x455 = x2 * x66
            x456 = x134 * x426
            x457 = x17 * x456
            x459 = x25 * x3**(-6.5_rp)
            x460 = 2835.0_rp * x459
            x462 = 13.5_rp * x42
            x470 = -x391 + x51
            x471 = -x199 + x200 + x319 - x394 - x395 + x49
            x473 = 51.0_rp * x1 * x12 - 35.0_rp * x1 * x174 * x24 + 35.0_rp * x1 * x174 + x1 * x24 * x248 - x1 * x248 + 11.0_rp * &
&x1 * x250 - 16.0_rp * x1 * x252 - 5.0_rp * x1 * x68 + x470 + x471
            x474 = 4.5_rp * x86
            x475 = 4.5_rp * x192
            x477 = 1.5_rp * x86
            x478 = x270 * x45
            x480 = x24 * x350
            x481 = x174 * x45
            x482 = x1 * x481
            x483 = x24 * x379
            x485 = -x307 + x400
            x488 = 105.0_rp * x292
            x489 = x24 * x488
            x490 = x107 - x108
            x494 = x234 * x332
            x495 = x1 * x299
            x496 = 472.5_rp * x495
            x499 = x7 * x94
            x500 = x155 * x499
            x502 = x15 * x72
            x505 = x40 * x72
            x506 = 18.0_rp * x270
            x508 = x187 * x362
            x509 = x508 * x88
            x510 = x187 * x331
            x511 = x174 * x510
            x514 = 31.5_rp * x179
            x515 = x15 * x185
            x516 = x1 * x515
            x518 = 4.5_rp * x42
            x522 = 4.5_rp * x367
            x526 = x1 * x24 * x49 - x1 * x44 - x1 * x47 - x321 + x322 + x359 + x58
            x527 = x333 * x526
            x528 = x1 * x64
            x529 = x1 * x146
            x536 = x1 * x69
            x537 = x374 * x45
            x538 = x1 * x97
            x539 = 4.0_rp * x538
            x540 = x24 * x538
            x542 = x1 * x81
            x543 = x24 * x542
            x545 = x1 * x80 - 9.0_rp * x123 + x228 + 9.0_rp * x31 * x45 - x45 * x55 - x528 + x529 - x536 - x537 - x539 + x540 + &
&x543 + x76
            x547 = x332 * x39
            x548 = x21 * x94
            x551 = 7.5_rp * x59
            x552 = x187 * x270
            x557 = x462 * x63
            x558 = x15 * x331
            x560 = 22.5_rp * x177 * x558
            x565 = x174 * x342
            x566 = x1 * x565
            x567 = 1575.0_rp * x338
            x569 = x270 * x342
            x571 = x1 * x160
            x575 = 157.5_rp * x168
            x582 = xi * yi
            x583 = 315.0_rp * x292
            x585 = x130 * x287
            x589 = 33.0_rp * x250
            x591 = 153.0_rp * x12
            x593 = 48.0_rp * x252
            x595 = x470 - x528 + x529
            x596 = x1 * x589 + x1 * x591 - x1 * x593 - x24 * x308 + x308 + x471 - x542 + x595
            x598 = x191 * x332
            x603 = 3.0_rp * x192
            x604 = x187 * x205
            x605 = x177 * x558
            x611 = x1 * x80 - x121 + x44 + x47 - x50 - x536 - x537 - x539 + x540 + x543 + x595
            x622 = xi * zi
            x623 = 315.0_rp * x427
            x624 = x416 * x45
            x628 = x24 * x445
            x629 = x2 * x481
            x630 = x24 * x456
            x632 = x2 * x299
            x633 = 472.5_rp * x632
            x635 = x2 * x81
            x639 = x2 * x64
            x641 = x146 * x2
            x643 = x2 * x307 + x2 * x589 + x2 * x591 - x2 * x593 - x432 + x470 + x471 - x635 - x639 + x641
            x649 = x187 * x443
            x650 = x174 * x649
            x651 = x187 * x453
            x652 = x651 * x88
            x654 = x2 * x397
            x657 = x187 * x416
            x663 = x2 * x69
            x664 = x45 * x455
            x665 = x2 * x97
            x667 = x24 * x665
            x671 = -x121 + x2 * x80 + x24 * x635 + x44 + x47 + x470 - x50 - x639 + x641 - x663 - x664 - 4.0_rp * x665 + x667
            x674 = x2 * x500
            x680 = x2 * x515
            x685 = x2 * x565
            x687 = x160 * x2
            x688 = 157.5_rp * x687
            x697 = 51.0_rp * x12 * x2 - 35.0_rp * x174 * x2 * x24 + 35.0_rp * x174 * x2 + x2 * x24 * x248 - x2 * x248 + 11.0_rp * &
&x2 * x250 - 16.0_rp * x2 * x252 - 5.0_rp * x2 * x68 + x470 + x471
            x701 = 105.0_rp * x427
            x706 = -x2 * x44 - x2 * x47 + x2 * x50 - x436 + x437 + x450 + x58
            x707 = x333 * x706
            x708 = -9.0_rp * x123 + x2 * x80 + x228 + x24 * x635 + 9.0_rp * x31 * x45 - x45 * x55 - x639 + x641 - x663 - x664 - &
&4.0_rp * x665 + x667 + x76
            x709 = 22.5_rp * x15 * x177 * x443
            x710 = yi**4
            x712 = 3.0_rp * x234
            x717 = x110 * x710
            x719 = 30.0_rp * x270
            x720 = x11 * x710
            x725 = x45 * x710
            x727 = x134 * x710
            x733 = 90.0_rp * x68
            x736 = 105.0_rp * x717
            x738 = 90.0_rp * x270
            x739 = x24 * x738
            x740 = x24 * x736
            x742 = 9.0_rp * x123 + x205 * x24 - x205 - 9.0_rp * x31 * x45 + x45 * x55 - x76
            x744 = 40.0_rp * x92
            x745 = x710 * x95
            x747 = x710 * x92
            x748 = 10.0_rp * x24
            x749 = 105.0_rp * x24
            x753 = 18.0_rp * x45 * x7
            x754 = x185 * x502
            x757 = 5.0_rp * x40 * x548
            x758 = x7 * x710
            x764 = x187 * x295
            x766 = x15 * x45
            x769 = 9.0_rp * x192
            x771 = 12.0_rp * x367
            x772 = 30.0_rp * x177
            x776 = 52.5_rp * x175
            x778 = 22.5_rp * x40 * x72
            x782 = x295 * x45
            x786 = -x1 * x24 * x76 + 36.0_rp * x1 * x43 + x146 * x710 - x24 * x374 + x24 * x710 * x81 + x24 * x710 * x97 - 6.0_rp * x359 + x374 - x64 * x710 - x66 * x725 - x69 * x710 + x710 * x80 - 4.0_rp * x710 * x97 + 18.0_rp * x782 + x83
            x789 = ai**5 * x155 * x29 * x39
            x790 = 2205.0_rp * x338
            x793 = 210.0_rp * x17
            x794 = x63 * x793
            x797 = x1 * x526
            x802 = 3780.0_rp * x459
            x804 = 35.0_rp * x292
            x807 = x332 * x405
            x809 = 210.0_rp * x350
            x810 = 135.0_rp * x478
            x811 = x1 * x744
            x812 = x1 * x95
            x813 = 5.0_rp * x812
            x814 = x24 * x812
            x815 = x1 * x92
            x816 = x748 * x815
            x817 = x283 * x45
            x818 = 105.0_rp * x480
            x819 = 12.0_rp * x97
            x820 = 54.0_rp * x78
            x821 = 45.0_rp * x256
            x827 = -x24 * x45 * x66 - 3.0_rp * x24 * x97 - x488 + x489 + x490 + x733 - x809 - x810 - x811 - x813 + x814 + x816 &
&+ x817 + x818 + x819 + x820 - x821
            x830 = 30.0_rp * x179
            x831 = 2.0_rp * x499
            x835 = x0 * x604
            x838 = 52.5_rp * x347
            x839 = x187 * x393
            x843 = 22.5_rp * x172
            x849 = x1 * x2
            x850 = x134 * x849
            x856 = x2 * x292
            x859 = 3.0_rp * x2
            x863 = x187 * x429
            x865 = x430 * x766
            x867 = x2 * x529
            x871 = x2 * x350
            x872 = x2 * x24
            x875 = x2 * x489
            x876 = x528 - x529
            x893 = 135.0_rp * x624
            x895 = x2 * x95
            x897 = x24 * x895
            x898 = x2 * x92
            x902 = -x2 * x744 - x24 * x45 * x66 + x24 * x701 - 3.0_rp * x24 * x97 + x422 * x45 - 210.0_rp * x445 + x490 + 105.0_rp * x628 - x701 + x733 + x748 * x898 + x819 + x820 - x821 - x893 - 5.0_rp * x895 + x897
            x904 = 30.0_rp * x416
            x907 = zi**4
            x909 = x45 * x907
            x911 = x110 * x907
            x913 = x11 * x907
            x918 = x907 * x95
            x922 = x108 * x909 - x2 * x24 * x733 + 180.0_rp * x2 * x68 + 108.0_rp * x2 * x78 - 36.0_rp * x2 * x79 - 135.0_rp * x23&
                  & * x909 - 90.0_rp * x24 * x416 + 105.0_rp * x24 * x911 + x24 * x918 + 90.0_rp * x416 + 24.0_rp * x665 - 6.0_rp * x667 + x742 - x744&
                  & * x907 + x748 * x907 * x92 + x749 * x913 - 105.0_rp * x911 - 210.0_rp * x913 - 5.0_rp * x918
            x923 = x7 * x907
            x929 = x146 * x907 - x2 * x24 * x76 + 36.0_rp * x2 * x43 - x24 * x455 + x24 * x81 * x907 + x24 * x907 * x97 + 18.0_rp * x429 * x45 - x45 * x66 * x907 - 6.0_rp * x450 + x455 - x64 * x907 - x69 * x907 + x80 * x907 + x83 - 4.0_rp * x907 * x97
            x935 = 1050.0_rp * x292
            x941 = x130 * x710
            x942 = 420.0_rp * x138
            x944 = x24 * x941
            x945 = x138 * x710
            x946 = 420.0_rp * x151
            x947 = x24 * x727
            x958 = -150.0_rp * x1 * x256 + 300.0_rp * x1 * x68 - 60.0_rp * x1 * x79 + x1 * x99 + x108 * x725 + x125 - 135.0_rp * x&
&23 * x725 - 150.0_rp * x24 * x270 + x24 * x745 + 150.0_rp * x270 + 40.0_rp * x538 - 10.0_rp * x540 - x710 * x744 + x720 * x749 - 2&
&10.0_rp * x720 - x736 + x740 - 5.0_rp * x745 + x747 * x748
            x961 = 5197.5_rp * x459
            x963 = 2835.0_rp * x338
            x965 = x187 * x526 * x719
            x966 = 630.0_rp * x292
            x969 = -x107 + x108 + x24 * x45 * x66 + 3.0_rp * x24 * x97 - x733 - x819 - x820 + x821
            x975 = 45.0_rp * x165 + x172 * x518 + x575
            x983 = x2 * x518
            x984 = x292 * x649
            x1010 = yi**6
            x1011 = -30.0_rp * x297 + 35.0_rp * x758 + 3.0_rp
            x1013 = 5.0_rp * x297 - 3.0_rp
            x1015 = 3.0_rp * x297 - 1.0_rp
            x1019 = x1 * x68
            x1020 = x1 * x78
            x1022 = x23 * x725
            x1026 = x1 * x256
            x1040 = -2.0_rp * x1 * x123 + 5.0_rp * x1 * x43 + x120 - x31 + x321 - x322 - x34 + x782
            x1045 = 27.0_rp * x1019 + 12.0_rp * x1020 - 12.0_rp * x1026 + x208 - x322 * x45 + x538 + x876
            x1060 = 2.0_rp * x1 * x123 - 5.0_rp * x1 * x43 - x2 * x49 + x2 * x528 + x317 * x849 + x318 * x849 + x320 - x321 + x3&
&22 + x34 - x410 * x849 - x436 + x437 - x782 - x867
            x1067 = x2 * x478
            x1069 = x2 * x480
            x1081 = x1 * x907
            x1087 = zi**6
            T(i, 1) = -x0 * x128 * (x0 * x99 - x100 * x102 + 300.0_rp * x100 - x103 * x104 + x108 * x70 + 105.0_rp * x109 + x112&
& * x24 - x112 + x125 - 150.0_rp * x24 * x88 + 10.0_rp * x24 * x93 + x24 * x96 - 10.0_rp * x24 * x98 + 150.0_rp * x88 - 210.0_rp * x&
&90 - 135.0_rp * x91 - 40.0_rp * x93 - 5.0_rp * x96 + 40.0_rp * x98) - 30.0_rp * x0 * x42 * x59 * x61 - x157 * (x100 * x145 - 1350.0&
&_rp * x100 - x102 * x93 - 810.0_rp * x103 - 1575.0_rp * x109 + 420.0_rp * x110 * x151 * x20 - 1575.0_rp * x111 * x24 + 1575.0_rp * x&
&111 + x130 * x20 * x24 - x130 * x20 - x135 * x20 - x137 * x20 - 420.0_rp * x138 * x20 - x141 * x20 + x143 * x20 * x24 - 6.0_rp *&
& x143 * x20 - x145 * x91 + x147 * x20 + x149 * x20 * x24 + x150 * x77 + x153 * x20 * x24 + x154 + 675.0_rp * x24 * x88 - 15.0_rp&
& * x24 * x96 + 45.0_rp * x24 * x98 - 675.0_rp * x88 + 3150.0_rp * x90 + 2025.0_rp * x91 + 600.0_rp * x93 + 75.0_rp * x96 - 180.0_rp *&
& x98) - 45.0_rp * x18 * (-70.0_rp * x5 + 63.0_rp * x8 + 15.0_rp) - x27 * (x20 * x22 + 105.0_rp * x5 - 315.0_rp * x8 - 5.0_rp) - 22.5_rp * x28 * x30 * x40 - x63 * x84 * x87
            T(i, 2) = -xi * yi * (-2.5_rp * ai * x231 * x31 + 2.5_rp * x15 * x191 * x84 + 7.5_rp * x15 * x30 * x84 - x157 * (-x0 &
&* x241 - x104 * x243 + 1050.0_rp * x111 * x45 - x119 * x45 - 1950.0_rp * x12 + 35.0_rp * x140 * x6 + x143 * x6 + 1050.0_rp * x174 &
&* x24 - 1050.0_rp * x174 + x237 - x239 + 40.0_rp * x24 * x240 - x24 * x248 * x6 * x94 + 300.0_rp * x24 * x250 - x24 * x253 - 220.&
&0_rp * x240 + 285.0_rp * x243 + 1785.0_rp * x244 - 840.0_rp * x245 - 315.0_rp * x247 + 900.0_rp * x252 + x253 - 180.0_rp * x256 - 105&
&0.0_rp * x45 * x88 + 405.0_rp * x68 + 15.0_rp * x97 + x99) + 150.0_rp * x159 * x161 + 210.0_rp * x163 * (9.0_rp * x5 - 5.0_rp) + 37.5&
&_rp * x165 * x28 + x167 * (-70.0_rp * x5 + 63.0_rp * x8 + 15.0_rp) + x169 * (-70.0_rp * x5 + 63.0_rp * x8 + 15.0_rp) + x170 * x172 * &
&x28 + 150.0_rp * x174 * x175 + 75.0_rp * x178 * x61 + 15.0_rp * x180 * x61 - x185 * x186 * x61 + x186 * x59 * x63 + x189 + 5.0_rp &
&* x192 * x59 * x63 - x209 * x210 * x63 + x233 * (x0 * x99 - x100 * x102 + 300.0_rp * x100 - x103 * x104 + x108 * x70 + 105.0_rp &
&* x109 + x112 * x24 - x112 + x125 - 150.0_rp * x24 * x88 + 10.0_rp * x24 * x93 + x24 * x96 - 10.0_rp * x24 * x98 + 150.0_rp * x88 &
&- 210.0_rp * x90 - 135.0_rp * x91 - 40.0_rp * x93 - 5.0_rp * x96 + 40.0_rp * x98) + x235 * (x0 * x99 - x100 * x102 + 300.0_rp * x100&
& - x103 * x104 + x108 * x70 + 105.0_rp * x109 + x112 * x24 - x112 + x125 - 150.0_rp * x24 * x88 + 10.0_rp * x24 * x93 + x24 * x96&
& - 10.0_rp * x24 * x98 + 150.0_rp * x88 - 210.0_rp * x90 - 135.0_rp * x91 - 40.0_rp * x93 - 5.0_rp * x96 + 40.0_rp * x98))
            T(i, 3) = -xi * zi * (-2.5_rp * ai * x231 * x31 + 2.5_rp * x15 * x191 * x84 + 7.5_rp * x15 * x30 * x84 - x157 * (-x0 &
&* x241 - x104 * x243 + 1050.0_rp * x111 * x45 - x119 * x45 - 1950.0_rp * x12 + 35.0_rp * x140 * x6 + x143 * x6 + 1050.0_rp * x174 &
&* x24 - 1050.0_rp * x174 + x237 - x239 + 40.0_rp * x24 * x240 - x24 * x248 * x6 * x94 + 300.0_rp * x24 * x250 - x24 * x253 - 220.&
&0_rp * x240 + 285.0_rp * x243 + 1785.0_rp * x244 - 840.0_rp * x245 - 315.0_rp * x247 + 900.0_rp * x252 + x253 - 180.0_rp * x256 - 105&
&0.0_rp * x45 * x88 + 405.0_rp * x68 + 15.0_rp * x97 + x99) + 150.0_rp * x159 * x161 + 210.0_rp * x163 * (9.0_rp * x5 - 5.0_rp) + 37.5&
&_rp * x165 * x28 + x167 * (-70.0_rp * x5 + 63.0_rp * x8 + 15.0_rp) + x169 * (-70.0_rp * x5 + 63.0_rp * x8 + 15.0_rp) + x170 * x172 * &
&x28 + 150.0_rp * x174 * x175 + 75.0_rp * x178 * x61 + 15.0_rp * x180 * x61 - x185 * x186 * x61 + x186 * x59 * x63 + x189 + 5.0_rp &
&* x192 * x59 * x63 - x209 * x210 * x63 + x233 * (x0 * x99 - x100 * x102 + 300.0_rp * x100 - x103 * x104 + x108 * x70 + 105.0_rp &
&* x109 + x112 * x24 - x112 + x125 - 150.0_rp * x24 * x88 + 10.0_rp * x24 * x93 + x24 * x96 - 10.0_rp * x24 * x98 + 150.0_rp * x88 &
&- 210.0_rp * x90 - 135.0_rp * x91 - 40.0_rp * x93 - 5.0_rp * x96 + 40.0_rp * x98) + x235 * (x0 * x99 - x100 * x102 + 300.0_rp * x100&
& - x103 * x104 + x108 * x70 + 105.0_rp * x109 + x112 * x24 - x112 + x125 - 150.0_rp * x24 * x88 + 10.0_rp * x24 * x93 + x24 * x96&
& - 10.0_rp * x24 * x98 + 150.0_rp * x88 - 210.0_rp * x90 - 135.0_rp * x91 - 40.0_rp * x93 - 5.0_rp * x96 + 40.0_rp * x98))
            T(i, 4) = ai * x205 * x40 * x63 + ai * x231 * x295 + 2.0_rp * x0 * x192 * x59 + x0 * x363 * x59 - x1 * x189 - x1 * &
&x329 - x1 * x335 - x1 * x340 - x1 * x341 - x1 * x344 - x1 * x349 + x1 * x372 + x1 * x377 - x1 * x384 + x127 * x327 * x63 - x15&
& * x287 * x353 * x59 - x15 * x287 * x368 * x59 - x155 * x336 * x84 + x157 * (-1062.0_rp * x1 * x12 + 630.0_rp * x1 * x174 * x24 &
&- 630.0_rp * x1 * x174 + x1 * x219 * x45 + x1 * x237 - x1 * x239 - 24.0_rp * x1 * x24 * x243 - 42.0_rp * x1 * x240 + 162.0_rp * x1&
& * x243 + 1665.0_rp * x1 * x244 - 720.0_rp * x1 * x245 + 825.0_rp * x1 * x246 - 210.0_rp * x1 * x247 - 432.0_rp * x1 * x250 + 432.0&
&_rp * x1 * x252 + x1 * x261 + x1 * x263 - x1 * x268 + 162.0_rp * x100 + 72.0_rp * x103 + 90.0_rp * x109 + x112 * x24 - x112 - 72.0&
&_rp * x196 + x208 + x214 - x219 + x223 * x91 - x24 * x45 * x67 + 4.0_rp * x24 * x93 + 9.0_rp * x270 * x6 * x94 + x282 - x283 - 19&
&5.0_rp * x90 - 105.0_rp * x91 - 22.0_rp * x93 - x96 + 6.0_rp * x98) - 300.0_rp * x159 * x287 * x299 + x175 * x389 + 30.0_rp * x18 * &
&x61 + 3.0_rp * x191 * x40 * x63 + x209 * x287 * x369 + x209 * x287 * x373 + x231 * x297 * x45 + x233 * x84 + x235 * x84 + x244 &
&* x381 - 52.5_rp * x28 * x293 - x28 * x331 * x333 - 15.0_rp * x28 * x351 + x28 * x385 - x282 * x345 - x287 * x348 - x291 * (-12.&
&0_rp * x1 * x4 - x158 + x287 * x288 + 3.0_rp) + 2.0_rp * x313 * x314 - x321 * x366 - x362 * x364 + x374 * x376 - x378 * x380 - x3&
&79 * x382
            T(i, 5) = -x414 * (-ai * x231 * x31 + x0 * x134 * x382 + x0 * x15 * x353 * x59 + x0 * x15 * x368 * x59 - x0 * x209&
& * x369 - x0 * x209 * x373 + x0 * x348 + x107 * x345 - x128 * x397 * x63 + 300.0_rp * x159 * x163 + x161 * x378 + 180.0_rp * x16&
&3 * (x158 - 2.0_rp) + 15.0_rp * x165 * x28 + x169 * x28 + x172 * x28 * x333 + x189 - x231 * x234 - 2.0_rp * x32 * x404 + x329 + x&
&335 + x340 + x341 + x344 + x349 + x364 * x393 + x366 * x51 - x372 - x376 * x66 - x377 + x384 + x406 * x84 - x413 * (-354.0_rp *&
& x12 + 3.0_rp * x140 * x6 - x146 + x223 * x250 - 8.0_rp * x24 * x243 - x24 * x408 * x6 - 14.0_rp * x240 + 54.0_rp * x243 + 555.0_rp&
& * x244 - 240.0_rp * x245 + 275.0_rp * x246 - 70.0_rp * x247 - 144.0_rp * x250 + 144.0_rp * x252 + x317 + x318 - x407 + x408 * x6 -&
& x410 + x411 + x64))
            T(i, 6) = ai * x205 * x40 * x63 + ai * x231 * x429 + 2.0_rp * x0 * x192 * x59 + x0 * x363 * x59 + x127 * x15 * x440&
& * x63 - x15 * x353 * x426 * x59 - x15 * x368 * x426 * x59 - x155 * x444 * x84 + x157 * (162.0_rp * x100 + 72.0_rp * x103 + 90.0&
&_rp * x109 + x112 * x24 - x112 - 1062.0_rp * x12 * x2 + 630.0_rp * x174 * x2 * x24 - 630.0_rp * x174 * x2 - 72.0_rp * x196 + x2 * x&
&219 * x45 + x2 * x237 - x2 * x239 - 24.0_rp * x2 * x24 * x243 - 42.0_rp * x2 * x240 + 162.0_rp * x2 * x243 + 1665.0_rp * x2 * x244&
& - 720.0_rp * x2 * x245 + 825.0_rp * x2 * x246 - 210.0_rp * x2 * x247 - 432.0_rp * x2 * x250 + 432.0_rp * x2 * x252 + x2 * x261 + x&
&2 * x263 - x2 * x268 + x208 + x214 - x219 + x223 * x91 - x24 * x45 * x67 + 4.0_rp * x24 * x93 + 9.0_rp * x416 * x6 * x94 + x421 &
&- x422 - 195.0_rp * x90 - 105.0_rp * x91 - 22.0_rp * x93 - x96 + 6.0_rp * x98) - 300.0_rp * x159 * x299 * x426 + x175 * x389 + 30.0&
&_rp * x18 * x61 - x189 * x2 + 3.0_rp * x191 * x40 * x63 - x2 * x329 - x2 * x335 - x2 * x340 - x2 * x341 - x2 * x344 - x2 * x349 &
&+ x2 * x372 + x2 * x377 - x2 * x384 + x209 * x369 * x426 + x209 * x373 * x426 + x231 * x430 * x45 + x233 * x84 + x235 * x84 + &
&x244 * x381 - x28 * x333 * x443 + x28 * x385 - 52.5_rp * x28 * x428 - 15.0_rp * x28 * x446 - x291 * (-x158 - 12.0_rp * x2 * x4 + &
&x288 * x426 + 3.0_rp) + 2.0_rp * x314 * x435 - x345 * x421 - x348 * x426 - x364 * x453 - x366 * x436 + x376 * x455 - x378 * x457&
& - x382 * x456
            T(i, 7) = -x582 * (7.5_rp * x1 * x15 * x179 * x59 - x1 * x209 * x475 - x1 * x209 * x518 + 85.5_rp * x1 * x347 + x1 *&
& x500 * x59 + x1 * x547 * x548 - 315.0_rp * x161 - 630.0_rp * x163 - 22.5_rp * x165 * x61 - 67.5_rp * x165 * x63 - x172 * x557 + 1&
&57.5_rp * x175 * x292 - 67.5_rp * x178 - 31.5_rp * x180 + x185 * x462 + x185 * x475 - 4.5_rp * x192 * x59 - x209 * x332 * x336 + x&
&209 * x477 + x209 * x494 - 67.5_rp * x270 * x375 + x287 * x460 + x287 * x567 - x313 * x477 - x313 * x494 + x321 * x502 * x59 - &
&x326 * x462 - x326 * x475 - x332 * x405 * x59 + x362 * x557 - 4.5_rp * x367 * x40 - 27.0_rp * x383 - x413 * (-35.0_rp * x1 * x151&
& * x174 + 19.0_rp * x138 * x287 - x24 * x287 * x408 + x261 + x263 - x268 + x287 * x408 - x302 - x303 - x304 + x305 - 153.0_rp * &
&x350 + 507.0_rp * x379 + x401 - 33.0_rp * x478 + 48.0_rp * x480 + 192.0_rp * x482 - 192.0_rp * x483 + x485 - x488 + x489 + x490) - &
&4.5_rp * x42 * x59 - x473 * x474 + x477 * x545 * x63 + x496 * x61 + x505 * x506 + 27.0_rp * x509 + 45.0_rp * x511 - x514 * x516 -&
& x516 * x522 + x527 * x61 + x551 * x552 + x560 * x61 + 270.0_rp * x566 + 67.5_rp * x569 * x63 + 157.5_rp * x571 * x61 + 157.5_rp *&
& x571 * x63 - x575 * x61)
            T(i, 8) = x622 * (-7.5_rp * x1 * x15 * x179 * x59 - 18.0_rp * x1 * x188 * x393 + x1 * x209 * x475 + x1 * x209 * x518&
& - 85.5_rp * x1 * x347 - x1 * x393 * x604 * x63 + x1 * x397 * x603 + x1 * x397 * x604 - x1 * x500 * x59 - x1 * x547 * x548 + x1&
&5 * x297 * x45 * (-x107 + x108 - x261 - x263 + x268 + x302 + x303 + x304 - x305 + x307 - x400 - x401) + x157 * (57.0_rp * x138 &
&* x287 - x151 * x308 + x24 * x583 - x24 * x585 + x261 + x263 - x268 - x302 - x303 - x304 + x305 - 459.0_rp * x350 + 1521.0_rp * &
&x379 + x401 - 99.0_rp * x478 + 144.0_rp * x480 + 576.0_rp * x482 - 576.0_rp * x483 + x485 + x490 - x583 + x585) + 105.0_rp * x161 +&
& 210.0_rp * x163 + 22.5_rp * x165 * x63 + x167 * x61 + x169 * x61 + x172 * x518 * x63 - 157.5_rp * x175 * x292 + 22.5_rp * x178 + &
&10.5_rp * x180 - x185 * x518 - x185 * x598 + x191 * x332 * x59 - x209 * x233 - x209 * x235 + x209 * x332 * x336 + x233 * x313 +&
& x235 * x313 + 67.5_rp * x270 * x375 - x287 * x460 - x287 * x567 + x295 * x404 + x30 * x332 * x59 - x321 * x502 * x59 + x326 * &
&x518 + x326 * x598 - x333 * x362 * x61 - x342 * x528 * x61 - x362 * x518 * x63 + x367 * x547 + 9.0_rp * x383 + x406 * x59 + x47&
&7 * x596 - x477 * x611 * x63 - x496 * x61 - x505 * x506 - 9.0_rp * x509 - 15.0_rp * x511 + x514 * x516 + x516 * x522 - x551 * x5&
&52 - 300.0_rp * x566 - 67.5_rp * x569 * x63 - 157.5_rp * x571 * x61 - 157.5_rp * x571 * x63 - 7.5_rp * x605 * x61)
            T(i, 9) = x582 * (-7.5_rp * x15 * x177 * x443 * x61 - x15 * x179 * x2 * x551 + x15 * x430 * x45 * (-x107 + x108 - x&
&261 - x263 + x268 + x302 + x303 + x304 - x305 + x307 - x400 - x401) + x157 * (-x130 * x24 * x426 + x130 * x426 + 57.0_rp * x138&
& * x426 + x24 * x623 + x261 + x263 - x268 - x302 - x303 - x304 + x305 + x401 - x432 * x45 - 459.0_rp * x445 + 1521.0_rp * x456 +&
& x485 + x490 - x623 - 99.0_rp * x624 + 144.0_rp * x628 + 576.0_rp * x629 - 576.0_rp * x630) + 105.0_rp * x161 + 210.0_rp * x163 + 22&
&.5_rp * x165 * x63 + x167 * x61 + x169 * x61 + x172 * x518 * x63 - 157.5_rp * x175 * x427 + 22.5_rp * x178 + 10.5_rp * x180 - x185&
& * x518 - x185 * x598 - 18.0_rp * x188 * x2 * x393 + x191 * x332 * x59 + x2 * x209 * x475 + x2 * x209 * x518 - 85.5_rp * x2 * x3&
&47 - x2 * x393 * x604 * x63 - x2 * x547 * x548 - x209 * x233 - x209 * x235 + x209 * x332 * x444 + x233 * x435 + x235 * x435 + &
&x30 * x332 * x59 - x333 * x453 * x61 - 67.5_rp * x342 * x416 * x63 - x342 * x61 * x639 + x367 * x547 + 67.5_rp * x375 * x416 + 9&
&.0_rp * x383 + x404 * x429 + x406 * x59 - 18.0_rp * x416 * x505 - x426 * x460 - x426 * x567 - x436 * x502 * x59 + x440 * x518 + &
&x440 * x598 - x453 * x518 * x63 - x477 * x63 * x671 + x477 * x643 + x514 * x680 + x522 * x680 - x551 * x657 - x59 * x674 + x60&
&3 * x654 + x604 * x654 - x61 * x633 - x61 * x688 - x63 * x688 - 15.0_rp * x650 - 9.0_rp * x652 - 300.0_rp * x685)
            T(i, 10) = -x622 * (x15 * x179 * x2 * x551 - 315.0_rp * x161 - 630.0_rp * x163 - 22.5_rp * x165 * x61 - 67.5_rp * x165&
& * x63 - x172 * x557 + 157.5_rp * x175 * x427 - 67.5_rp * x178 - 31.5_rp * x180 + x185 * x462 + x185 * x475 - 4.5_rp * x192 * x59 &
&- x2 * x209 * x475 - x2 * x209 * x518 + 85.5_rp * x2 * x347 + x2 * x547 * x548 - x209 * x332 * x444 + x209 * x477 + x209 * x494&
& - x332 * x405 * x59 + 67.5_rp * x342 * x416 * x63 - 4.5_rp * x367 * x40 - 67.5_rp * x375 * x416 - 27.0_rp * x383 - x413 * (19.0_rp&
& * x138 * x426 - 35.0_rp * x151 * x174 * x2 - x24 * x408 * x426 + x24 * x701 + x261 + x263 - x268 - x302 - x303 - x304 + x305 +&
& x401 + x408 * x426 - 153.0_rp * x445 + 507.0_rp * x456 + x485 + x490 - 33.0_rp * x624 + 48.0_rp * x628 + 192.0_rp * x629 - 192.0_rp&
& * x630 - x701) + 18.0_rp * x416 * x505 - 4.5_rp * x42 * x59 + x426 * x460 + x426 * x567 - x435 * x477 - x435 * x494 + x436 * x5&
&02 * x59 - x440 * x462 - x440 * x475 + x453 * x557 - x474 * x697 + x477 * x63 * x708 - x514 * x680 - x522 * x680 + x551 * x657&
& - x575 * x61 + x59 * x674 + x61 * x633 + x61 * x688 + x61 * x707 + x61 * x709 + x63 * x688 + 45.0_rp * x650 + 27.0_rp * x652 + &
&270.0_rp * x685)
            T(i, 11) = ai * x1 * x205 * x327 - x0 * x710 * x790 + 45.0_rp * x1 * x180 - x1 * x214 * x508 + x1 * x326 * x769 + x&
&1 * x327 * x46 * x72 + 180.0_rp * x1 * x383 + 3.0_rp * x1 * x40 * x499 - 180.0_rp * x1 * x511 - 52.5_rp * x10 * x40 * x725 - x127 &
&* x327 - 3780.0_rp * x162 * x3**(-6.5_rp) * x710 - 90.0_rp * x174 * x25 + x175 * x282 - 45.0_rp * x18 - 12.0_rp * x188 * x797 - 4.5&
&_rp * x191 * x40 + 30.0_rp * x21 * x515 * x725 - x23 * x710 * x778 - x233 * x63 * x786 + x25 * x583 * x63 - x27 * x63 - x282 * x&
&510 * x63 + 1890.0_rp * x287 * x299 - x287 * x373 * x545 - 6.0_rp * x295 * x754 + 6.0_rp * x297 * x473 * x766 - 472.5_rp * x299 * &
&x63 * x710 - 4.5_rp * x30 * x40 - x314 * (-x1 * x24 * x733 + 180.0_rp * x1 * x68 + 108.0_rp * x1 * x78 - 36.0_rp * x1 * x79 + x108&
& * x725 - 135.0_rp * x23 * x725 + x24 * x745 + 24.0_rp * x538 - 6.0_rp * x540 - x710 * x744 + x720 * x749 - 210.0_rp * x720 - x736&
& + x738 - x739 + x740 + x742 - 5.0_rp * x745 + x747 * x748) - x327 * x712 + x331 * x604 * x63 + 90.0_rp * x351 * x63 - x363 * x6&
&3 * x797 - x374 * x375 + x374 * x505 + 990.0_rp * x380 + x389 * x510 - x405 * x547 + x413 * (-x0 * x24 * x408 * x710 + x0 * x40&
&8 * x710 + 443.0_rp * x0 * x727 - 306.0_rp * x1 * x12 - 66.0_rp * x1 * x250 + 96.0_rp * x1 * x252 - x1 * x407 + x1 * x411 + 93.0_rp&
& * x174 * x725 + 35.0_rp * x24 * x717 - x24 * x719 + x397 + x536 - 35.0_rp * x717 + x719 - 35.0_rp * x720 - 128.0_rp * x727 * x77)&
& + 6.0_rp * x473 * x764 + x508 * x67 + x515 * x710 * x771 + x515 * x710 * x772 + 2.0_rp * x515 * x758 * x94 - x516 * x753 - 420.&
&0_rp * x565 * x710 - x710 * x757 - x710 * x789 - x717 * x776 - x727 * x794
            T(i, 12) = -x414 * (x0 * x282 * x839 - x1 * x397 * x475 - x1 * x397 * x518 + 5.0_rp * x1 * x40 * x548 + x1 * x789 +&
& x1 * x838 - x130 * x162 - 495.0_rp * x161 - 45.0_rp * x165 * x63 - x172 * x518 * x63 - 22.5_rp * x178 - 22.5_rp * x180 + x185 * x&
&604 + x187 * x52 * x545 + 3.0_rp * x188 * x526 - x205 * x40 * x72 + x233 * x545 * x63 + x270 * x778 + x287 * x790 + x287 * x802&
& + x292 * x776 + x314 * x827 - x326 * x475 - x326 * x518 - x326 * x807 - x332 * x336 * x397 + x362 * x518 * x63 + 9.0_rp * x371&
& - x375 * x719 - 99.0_rp * x383 - x393 * x835 + x397 * x477 + x397 * x494 - x413 * (-x146 - x24 * x287 * x408 + x24 * x804 + x2&
&87 * x408 - 35.0_rp * x350 + 443.0_rp * x379 + 93.0_rp * x482 - 128.0_rp * x483 + x485 - x589 - x591 + x593 + x64 - x804 + x81) + &
&x46 * x754 - x473 * x477 - x473 * x494 - x477 * x596 - x494 * x596 + x496 * x63 - x499 * x547 + 54.0_rp * x509 + 90.0_rp * x511 &
&- x516 * x771 - x516 * x830 - x516 * x831 + x527 * x63 + x552 * x63 * x843 + x560 * x63 + 510.0_rp * x566 + 210.0_rp * x571 * x6&
&3 - x575 * x63 + x611 * x835)
            T(i, 13) = 7.5_rp * x1 * x180 + x1 * x333 * x440 - 3.0_rp * x1 * x371 + 30.0_rp * x1 * x383 + x1 * x39 * x500 + x1 * &
&x440 * x598 - 30.0_rp * x1 * x650 - 21.0_rp * x1 * x652 + x155 * x326 * x444 + x155 * x336 * x440 + x157 * (-x1 * x589 - x1 * x5&
&91 + x1 * x593 - x2 * x24 * x585 - x2 * x307 + 1329.0_rp * x2 * x379 - x2 * x488 + x2 * x585 - x2 * x589 - x2 * x591 + x2 * x59&
&3 + x24 * x308 - x308 - 384.0_rp * x379 * x872 + x397 + x432 + 279.0_rp * x481 * x849 + x542 + x635 + x639 - x641 - 105.0_rp * x8&
&71 + x875 + x876) + 7.5_rp * x175 * x270 + 7.5_rp * x175 * x416 - 15.0_rp * x18 + 7.5_rp * x180 * x2 - 60.0_rp * x188 * x393 * x849&
& + 6.0_rp * x192 * x397 * x849 - x2 * x270 * x778 - x2 * x287 * x790 - x2 * x287 * x802 + x2 * x326 * x333 + x2 * x326 * x598 -&
& x2 * x342 * x63 * x719 + x2 * x375 * x719 + 30.0_rp * x2 * x383 - x2 * x496 * x63 - 21.0_rp * x2 * x509 - 30.0_rp * x2 * x511 - &
&x233 * x326 - x233 * x440 - x233 * x63 * (-x1 * x24 * x49 + x1 * x44 + x1 * x47 + x172 + x2 * x44 + x2 * x47 - x2 * x50 - x2 *&
& x528 - x2 * x536 - x2 * x537 - x2 * x539 + x2 * x540 + x2 * x543 + x321 - x322 - x359 + x436 - x437 - x450 + x80 * x849 + x86&
&7) - x235 * x326 - x235 * x440 + x25 * x287 * x408 + x25 * x408 * x426 - 7.5_rp * x270 * x63 * x649 - x287 * x363 * x671 - x291&
& + 52.5_rp * x293 * x63 + 2.0_rp * x295 * x502 * x654 - x295 * x754 + x297 * x643 * x766 - x30 * x547 - x314 * (-x1 * x80 - x2 *&
& x488 - x2 * x80 - x2 * x809 - x2 * x810 - x2 * x811 - x2 * x813 + x2 * x814 + x2 * x816 + x2 * x817 + x2 * x818 - x24 * x635 &
&+ x393 + x536 + x537 + x539 - x540 - x543 + x639 - x641 + x663 + x664 + 4.0_rp * x665 - x667 + x875 + x876) - x321 * x375 + x32&
&1 * x505 - x321 * x63 * x651 + x331 * x333 * x63 + x333 * x443 * x63 + 15.0_rp * x351 * x63 + x363 * x397 * x849 - x363 * x426 &
&* x611 - x371 * x859 - x375 * x436 + 165.0_rp * x380 - x385 * x63 - x39 * x406 - x39 * x598 + x39 * x674 - 7.5_rp * x416 * x510 &
&* x63 + 52.5_rp * x428 * x63 - x429 * x754 + x436 * x505 - x436 * x508 * x63 + 15.0_rp * x446 * x63 + 165.0_rp * x457 + x508 * x5&
&2 + 3.0_rp * x510 * x88 + x515 * x771 * x849 + x515 * x830 * x849 + x515 * x831 * x849 + x52 * x651 - 540.0_rp * x565 * x849 + x&
&596 * x863 + x596 * x865 + x643 * x764 + 3.0_rp * x649 * x88 - x757 * x849 - x776 * x856 - x789 * x849 - x794 * x850 - x838 * x&
&849
            T(i, 14) = -x414 * (x0 * x421 * x839 - x130 * x162 - 495.0_rp * x161 - 45.0_rp * x165 * x63 - x172 * x518 * x63 - 22&
&.5_rp * x178 - 22.5_rp * x180 + x185 * x604 + x187 * x52 * x708 + 3.0_rp * x188 * x706 + x2 * x757 + x2 * x789 + x2 * x838 - x205&
& * x40 * x72 + x233 * x63 * x708 + x314 * x902 - x332 * x397 * x444 + 9.0_rp * x371 - x375 * x904 - 99.0_rp * x383 - x393 * x835&
& + x397 * x477 + x397 * x494 - x413 * (-x146 - x24 * x408 * x426 + 35.0_rp * x24 * x427 + x408 * x426 - 35.0_rp * x427 - 35.0_rp &
&* x445 + 443.0_rp * x456 + x485 - x589 - x591 + x593 + 93.0_rp * x629 - 128.0_rp * x630 + x64 + x81) + x416 * x778 + x426 * x790 &
&+ x426 * x802 + x427 * x776 - x440 * x475 - x440 * x518 - x440 * x807 + x453 * x518 * x63 + x46 * x754 - x475 * x654 - x477 * &
&x643 - x477 * x697 - x494 * x643 - x494 * x697 - x499 * x547 - x518 * x654 - x575 * x63 + x63 * x633 + x63 * x657 * x843 + 210&
&.0_rp * x63 * x687 + x63 * x707 + x63 * x709 + 90.0_rp * x650 + 54.0_rp * x652 + x671 * x835 - x680 * x771 - x680 * x830 - x680 *&
& x831 + 510.0_rp * x685)
            T(i, 15) = -x0 * x790 * x907 - 52.5_rp * x10 * x40 * x909 - x127 * x15 * x440 - x134 * x794 * x907 - x15 * x440 * x&
&712 - 3780.0_rp * x162 * x3**(-6.5_rp) * x907 - 90.0_rp * x174 * x25 + x175 * x421 - 45.0_rp * x18 + 45.0_rp * x180 * x2 - 12.0_rp *&
& x188 * x2 * x706 - 4.5_rp * x191 * x40 - x2 * x214 * x651 - x2 * x363 * x63 * x706 + 180.0_rp * x2 * x383 + x2 * x440 * x46 * x&
&502 + x2 * x440 * x604 + x2 * x440 * x769 - 180.0_rp * x2 * x650 + 30.0_rp * x21 * x515 * x909 - x23 * x778 * x907 - x233 * x63 &
&* x929 + x25 * x623 * x63 - x27 * x63 + 1890.0_rp * x299 * x426 - 472.5_rp * x299 * x63 * x907 - 4.5_rp * x30 * x40 - x314 * x922&
& - x373 * x426 * x708 - x375 * x455 + x389 * x649 + x40 * x499 * x859 - x405 * x547 + x413 * (443.0_rp * x0 * x134 * x907 - x0 &
&* x24 * x408 * x907 + x0 * x408 * x907 - 306.0_rp * x12 * x2 - 128.0_rp * x134 * x77 * x907 + 93.0_rp * x174 * x909 - 66.0_rp * x2&
& * x250 + 96.0_rp * x2 * x252 - x2 * x407 + x2 * x411 - x24 * x904 + 35.0_rp * x24 * x911 + x397 + x663 + x904 - 35.0_rp * x911 -&
& 35.0_rp * x913) - x421 * x63 * x649 - 6.0_rp * x429 * x754 + x443 * x604 * x63 + 90.0_rp * x446 * x63 + x455 * x505 + 990.0_rp * &
&x457 + x515 * x771 * x907 + x515 * x772 * x907 + 2.0_rp * x515 * x923 * x94 - 420.0_rp * x565 * x907 + x651 * x67 - x680 * x753 &
&+ 6.0_rp * x697 * x863 + 6.0_rp * x697 * x865 - x757 * x907 - x776 * x911 - x789 * x907
            T(i, 16) = -x582 * (x1 * x186 * x545 - x1 * x253 * x342 - x115 * x508 + x157 * (-x1 * x24 * x241 - x135 * x710 - x&
&136 * x717 - x141 * x710 + x143 * x24 * x710 - 6.0_rp * x143 * x710 + x146 * x710 * x94 - x150 + x24 * x253 - 450.0_rp * x24 * x&
&478 - 100.0_rp * x24 * x815 - x24 * x935 + 15.0_rp * x24 * x97 - x253 + 225.0_rp * x256 + 2100.0_rp * x350 + 1350.0_rp * x478 - 105&
&0.0_rp * x480 - 450.0_rp * x68 - x710 * x942 + x717 * x946 + x749 * x945 + 90.0_rp * x79 + 50.0_rp * x812 + 400.0_rp * x815 + x935 &
&- x941 + x944 + 945.0_rp * x947 - 60.0_rp * x97) + 225.0_rp * x165 + 787.5_rp * x168 + x170 * x786 - x186 * x526 - x210 * x545 + x&
&233 * x958 - x253 * x510 + 75.0_rp * x270 * x508 + 525.0_rp * x292 * x510 + 262.5_rp * x342 * x717 + x42 * x843 - 4725.0_rp * x495&
& + 75.0_rp * x526 * x552 - 2100.0_rp * x571 + x710 * x961 + x710 * x963 + 2.5_rp * x86 * (-x1 * x24 * x733 + 180.0_rp * x1 * x68 +&
& 108.0_rp * x1 * x78 - 36.0_rp * x1 * x79 + x108 * x725 - 135.0_rp * x23 * x725 + x24 * x745 + 24.0_rp * x538 - 6.0_rp * x540 - x71&
&0 * x744 + x720 * x749 - 210.0_rp * x720 - x736 + x738 - x739 + x740 + x742 - 5.0_rp * x745 + x747 * x748))
            T(i, 17) = -x622 * (-x1 * x253 * x342 + x1 * x373 * x545 + x1 * x604 * x611 - x107 * x510 - x128 * x611 + x15 * x3&
&93 * x710 * x772 + x157 * (-x104 * x815 - x135 * x710 - x136 * x717 - x141 * x710 + x143 * x24 * x710 - 6.0_rp * x143 * x710 + &
&x146 * x710 * x94 - 630.0_rp * x24 * x350 - 270.0_rp * x24 * x478 - x24 * x966 + 1260.0_rp * x350 + 810.0_rp * x478 - x710 * x942 &
&+ x717 * x946 + x749 * x945 + 30.0_rp * x812 - 6.0_rp * x814 + 240.0_rp * x815 - x941 + x944 + 945.0_rp * x947 + x966 + x969) + x3&
&33 * x786 + 472.5_rp * x342 * x717 - x374 * x839 - 2835.0_rp * x495 - x508 * x66 + x508 * x738 + x510 * x583 - 1260.0_rp * x571 +&
& x710 * x961 + x710 * x963 + 2.0_rp * x764 * x827 + x86 * (-x1 * x24 * x733 + 180.0_rp * x1 * x68 + 108.0_rp * x1 * x78 - 36.0_rp &
&* x1 * x79 + x108 * x725 - 135.0_rp * x23 * x725 + x24 * x745 + 24.0_rp * x538 - 6.0_rp * x540 - x710 * x744 + x720 * x749 - 210.&
&0_rp * x720 - x736 + x738 - x739 + x740 + x742 - 5.0_rp * x745 + x747 * x748) + x965 + x975)
            T(i, 18) = -x582 * (x1 * x518 * x671 + x157 * (x130 * x24 * x849 - x130 * x849 - x135 * x849 - x136 * x856 + x143 &
&* x24 * x849 - 6.0_rp * x143 * x849 + x149 * x24 * x849 + x153 * x24 * x849 - 75.0_rp * x2 * x270 * x94 - x223 * x898 - x24 * x6&
&23 - x24 * x893 - x437 * x94 + 630.0_rp * x445 + x488 - x489 + x623 + 405.0_rp * x624 - 315.0_rp * x628 + x809 + x810 + x811 + x8&
&13 - x814 - x816 - x817 - x818 - x849 * x942 + x856 * x946 + x867 * x94 + 15.0_rp * x895 + 120.0_rp * x898 + x969) + x2 * x282 *&
& x839 - x2 * x393 * x604 + x2 * x604 * x611 - x233 * x545 + x233 * (-x1 * x80 - x2 * x488 + x2 * x733 - x2 * x809 - x2 * x810 &
&- x2 * x811 - x2 * x813 + x2 * x814 + x2 * x816 + x2 * x817 + x2 * x818 + x2 * x819 + x2 * x820 - x2 * x821 - x24 * x664 + x42&
&1 - x422 + x536 + x537 + x539 - x540 - x543 - 3.0_rp * x667 + x742 + x875 + x876) - 112.5_rp * x342 * x416 + 577.5_rp * x342 * x8&
&56 - x362 * x518 + 67.5_rp * x416 * x508 + 157.5_rp * x427 * x510 - x453 * x604 - x477 * x671 + x477 * (-x1 * x80 - x2 * x488 - &
&x2 * x80 - x2 * x809 - x2 * x810 - x2 * x811 - x2 * x813 + x2 * x814 + x2 * x816 + x2 * x817 + x2 * x818 - x24 * x635 + x393 +&
& x536 + x537 + x539 - x540 - x543 + x639 - x641 + x663 + x664 + 4.0_rp * x665 - x667 + x875 + x876) - x496 + x518 * (-x1 * x24 &
&* x49 + x1 * x44 + x1 * x47 + x172 + x2 * x44 + x2 * x47 - x2 * x50 - x2 * x528 - x2 * x536 - x2 * x537 - x2 * x539 + x2 * x54&
&0 + x2 * x543 + x321 - x322 - x359 + x436 - x437 - x450 + x80 * x849 + x867) + 7.5_rp * x526 * x657 - x527 + x545 * x983 - x552&
& * x843 - x560 - 210.0_rp * x571 - 1417.5_rp * x632 + x651 * x719 - 630.0_rp * x687 - x709 + x827 * x863 + x849 * x961 + x849 * x&
&963 + x975 + 52.5_rp * x984)
            T(i, 19) = -x622 * (-x1 * x393 * x604 + x1 * x518 * x708 + x1 * x604 * x671 + x157 * (x130 * x24 * x849 - x130 * x&
&849 - x135 * x849 - x136 * x856 + x143 * x24 * x849 - 6.0_rp * x143 * x849 + x149 * x24 * x849 + x153 * x24 * x849 - 75.0_rp * x&
&2 * x270 * x94 + x2 * x744 - x223 * x815 - x24 * x583 - x24 * x701 - x24 * x810 - x322 * x94 + 630.0_rp * x350 - x422 * x45 + 2&
&10.0_rp * x445 + 405.0_rp * x478 - 315.0_rp * x480 + x583 - 105.0_rp * x628 + x701 - x748 * x898 + 15.0_rp * x812 + 120.0_rp * x815 &
&- x849 * x942 + x856 * x946 + x867 * x94 + x893 + 5.0_rp * x895 - x897 + x969) + x2 * x282 * x839 - x233 * x708 + x233 * (x1 * &
&x733 + x1 * x819 + x1 * x820 - x1 * x821 - x2 * x488 - x2 * x80 - x2 * x809 - x2 * x810 - x2 * x811 - x2 * x813 + x2 * x814 + &
&x2 * x816 + x2 * x817 + x2 * x818 - x24 * x537 - x24 * x635 + x282 - x283 - 3.0_rp * x540 + x639 - x641 + x663 + x664 + 4.0_rp *&
& x665 - x667 + x742 + x875) - 112.5_rp * x270 * x342 + 67.5_rp * x270 * x651 + 577.5_rp * x342 * x856 - x362 * x604 + 52.5_rp * x4&
&27 * x510 - x453 * x518 - x477 * x611 + x477 * (-x1 * x80 - x2 * x488 - x2 * x80 - x2 * x809 - x2 * x810 - x2 * x811 - x2 * x8&
&13 + x2 * x814 + x2 * x816 + x2 * x817 + x2 * x818 - x24 * x635 + x393 + x536 + x537 + x539 - x540 - x543 + x639 - x641 + x663&
& + x664 + 4.0_rp * x665 - x667 + x875 + x876) - 1417.5_rp * x495 + x508 * x904 + x518 * (-x1 * x24 * x49 + x1 * x44 + x1 * x47 +&
& x172 + x2 * x44 + x2 * x47 - x2 * x50 - x2 * x528 - x2 * x536 - x2 * x537 - x2 * x539 + x2 * x540 + x2 * x543 + x321 - x322 -&
& x359 + x436 - x437 - x450 + x80 * x849 + x867) + 7.5_rp * x552 * x706 - x560 - 630.0_rp * x571 + x611 * x983 - x633 - x657 * x8&
&43 - 210.0_rp * x687 - x707 - x709 + x764 * x902 + x849 * x961 + x849 * x963 + x975 + 157.5_rp * x984)
            T(i, 20) = -x582 * (-x107 * x649 - x128 * x671 + x15 * x393 * x772 * x907 + x157 * (-x104 * x898 + x130 * x24 * x9&
&07 - x130 * x907 - x135 * x907 - x136 * x911 - x141 * x907 + x143 * x24 * x907 - 6.0_rp * x143 * x907 + x146 * x907 * x94 + x14&
&9 * x24 * x907 + x153 * x24 * x907 - 630.0_rp * x24 * x427 - 630.0_rp * x24 * x445 - 270.0_rp * x24 * x624 + 630.0_rp * x427 + 126&
&0.0_rp * x445 + 810.0_rp * x624 + 30.0_rp * x895 - 6.0_rp * x897 + 240.0_rp * x898 - x907 * x942 + x911 * x946 + x969) + x187 * x70&
&6 * x904 - x2 * x253 * x342 + x2 * x373 * x708 + x2 * x604 * x671 + x333 * x929 + 472.5_rp * x342 * x911 + 90.0_rp * x416 * x651&
& - x455 * x839 + x623 * x649 - 2835.0_rp * x632 - x651 * x66 - 1260.0_rp * x687 + x86 * x922 + 2.0_rp * x863 * x902 + x907 * x961&
& + x907 * x963 + x975)
            T(i, 21) = -x622 * (-x115 * x651 + x157 * (x130 * x24 * x907 - x130 * x907 - x135 * x907 - x136 * x911 - x141 * x9&
&07 + x143 * x24 * x907 - 6.0_rp * x143 * x907 + x146 * x907 * x94 + x149 * x24 * x907 - x150 + x153 * x24 * x907 + x24 * x253 -&
& 1050.0_rp * x24 * x427 - 450.0_rp * x24 * x624 - 100.0_rp * x24 * x898 + 15.0_rp * x24 * x97 - x241 * x872 - x253 + 225.0_rp * x25&
&6 + 1050.0_rp * x427 + 2100.0_rp * x445 + 1350.0_rp * x624 - 1050.0_rp * x628 - 450.0_rp * x68 + 90.0_rp * x79 + 50.0_rp * x895 + 400&
&.0_rp * x898 - x907 * x942 + x911 * x946 - 60.0_rp * x97) + 225.0_rp * x165 + 787.5_rp * x168 + x170 * x929 + x186 * x2 * x708 - x&
&186 * x706 - x2 * x253 * x342 - x210 * x708 + x233 * (x108 * x909 + x125 - 150.0_rp * x2 * x256 + 300.0_rp * x2 * x68 - 60.0_rp *&
& x2 * x79 + x2 * x99 - 135.0_rp * x23 * x909 - 150.0_rp * x24 * x416 + 105.0_rp * x24 * x911 + x24 * x918 + 150.0_rp * x416 + 40.0&
&_rp * x665 - 10.0_rp * x667 - x744 * x907 + x748 * x907 * x92 + x749 * x913 - 105.0_rp * x911 - 210.0_rp * x913 - 5.0_rp * x918) - &
&x253 * x649 + 262.5_rp * x342 * x911 + 75.0_rp * x416 * x651 + x42 * x843 + 525.0_rp * x427 * x649 - 4725.0_rp * x632 + 75.0_rp * x&
&657 * x706 - 2100.0_rp * x687 + 2.5_rp * x86 * x922 + x907 * x961 + x907 * x963)
            T(i, 22) = -x1 * x128 * x958 - 22.5_rp * x1011 * x331 * x42 - 30.0_rp * x1013 * x42 * x797 - x1015 * x786 * x87 - x1&
&57 * (x1 * x150 * x24 + 420.0_rp * x1010 * x110 * x151 + x1010 * x130 * x24 - x1010 * x130 - x1010 * x135 - x1010 * x137 - x101&
&0 * x141 + x1010 * x143 * x24 - 6.0_rp * x1010 * x143 + x1010 * x147 + x1010 * x149 * x24 + x1010 * x153 * x24 - x1010 * x942 -&
& 1350.0_rp * x1019 - x102 * x747 - 810.0_rp * x1020 - x1022 * x145 + 2025.0_rp * x1022 + 675.0_rp * x1026 + x154 + 675.0_rp * x24 *&
& x270 - 1575.0_rp * x24 * x717 - 1575.0_rp * x24 * x720 - 15.0_rp * x24 * x745 - 675.0_rp * x270 - 180.0_rp * x538 + 45.0_rp * x540 &
&+ 1575.0_rp * x717 + 3150.0_rp * x720 + 75.0_rp * x745 + 600.0_rp * x747) - x27 * (x1010 * x22 + 105.0_rp * x297 - 315.0_rp * x758 -&
& 5.0_rp) - 45.0_rp * x351 * (-70.0_rp * x297 + 63.0_rp * x758 + 15.0_rp)
            T(i, 23) = -x414 * (37.5_rp * x1011 * x165 + x1011 * x170 * x172 - x1013 * x1040 * x186 + 15.0_rp * x1013 * x179 * x&
&558 + 75.0_rp * x1013 * x605 - x1015 * x1045 * x210 + x1015 * x186 * x526 + 5.0_rp * x1015 * x192 * x526 - x157 * (-x1 * x241 - &
&x104 * x945 - x119 * x45 + 35.0_rp * x140 * x710 + x143 * x710 - 315.0_rp * x151 * x717 - x24 * x248 * x710 * x94 - x24 * x253 +&
& 300.0_rp * x24 * x478 + x24 * x811 + x24 * x935 + x253 - 180.0_rp * x256 - 1050.0_rp * x270 * x45 - 1950.0_rp * x350 + 1050.0_rp *&
& x45 * x717 + 900.0_rp * x480 + 405.0_rp * x68 + 1785.0_rp * x727 - 220.0_rp * x815 - x935 + x941 - x944 + 285.0_rp * x945 - 840.0_rp * x947 + 15.0_rp * x97 + x99) + x167 * (-70.0_rp * x297 + 63.0_rp * x758 + 15.0_rp) + x169 * (-70.0_rp * x297 + 63.0_rp * x758 + 1&
&5.0_rp) + x170 * x786 + 2.5_rp * x192 * x786 + x233 * x958 + x235 * x958 + 150.0_rp * x292 * x510 + 210.0_rp * x495 * (9.0_rp * x29&
&7 - 5.0_rp) + 150.0_rp * x571 * (7.0_rp * x297 - 3.0_rp) - 2.5_rp * x86 * (-162.0_rp * x1019 - 72.0_rp * x1020 - x1022 * x223 + 105.0&
&_rp * x1022 + 72.0_rp * x1026 + x230 + x24 * x537 - 90.0_rp * x24 * x720 - 4.0_rp * x24 * x747 - 6.0_rp * x538 + 195.0_rp * x720 + x&
&736 - x738 + x739 - x740 + x745 + 22.0_rp * x747) + x965)
            T(i, 24) = -x1 * x1013 * x363 * x453 - 36.0_rp * x10 * x331 * x766 * x849 - x1011 * x333 * x443 + x1011 * x385 - 52&
&.5_rp * x1011 * x428 - 15.0_rp * x1011 * x446 - 60.0_rp * x1013 * x2 * x569 + 30.0_rp * x1013 * x351 - x1013 * x793 * x850 + x1015&
& * x1040 * x187 * x455 + 6.0_rp * x1015 * x1040 * x192 * x2 + x1015 * x1060 * x128 - 21.0_rp * x1015 * x179 * x2 * x558 - x1015 &
&* x331 * x436 * x502 + x1015 * x331 * x603 + x1015 * x331 * x604 - x1015 * x421 * x510 + 36.0_rp * x1040 * x2 * x552 + x1045 * &
&x369 * x849 + x1045 * x373 * x849 - x15 * x353 * x526 * x849 - x15 * x368 * x526 * x849 - x155 * x444 * x786 + x157 * (162.0_rp&
& * x1019 + 72.0_rp * x1020 + x1022 * x223 - 105.0_rp * x1022 - 72.0_rp * x1026 - 432.0_rp * x1067 + 432.0_rp * x1069 + x2 * x24 * x&
&966 + x2 * x261 + x2 * x263 - x2 * x268 + x2 * x45 * x739 + 1665.0_rp * x2 * x727 - 42.0_rp * x2 * x815 + x2 * x941 - x2 * x944 &
&+ 162.0_rp * x2 * x945 - 720.0_rp * x2 * x947 - x2 * x966 + x208 - 210.0_rp * x24 * x427 * x725 - x24 * x537 + 90.0_rp * x24 * x72&
&0 + 4.0_rp * x24 * x747 + 9.0_rp * x416 * x710 * x94 + x421 - x422 + 825.0_rp * x427 * x725 + 6.0_rp * x538 - 195.0_rp * x720 - x73&
&6 + x738 - x739 + x740 - x745 - 22.0_rp * x747 - 1062.0_rp * x871 - 24.0_rp * x872 * x945) + 2.0_rp * x192 * x797 - x2 * x333 * x7&
&86 - 840.0_rp * x2 * x338 * x710 - 300.0_rp * x2 * x495 * (7.0_rp * x297 - 3.0_rp) - x2 * x598 * x786 - x2 * x965 + x233 * x786 + &
&x235 * x786 - 30.0_rp * x293 * (-12.0_rp * x2 * x4 + x288 * x849 - 7.0_rp * x297 + 3.0_rp) - 120.0_rp * x342 * x427 * x710 + x363 *&
& x797 + x381 * x727 - x381 * x850 * (7.0_rp * x297 - 3.0_rp) + x506 * x510 - 180.0_rp * x510 * x856 + 2.0_rp * x764 * (-27.0_rp * x&
&1019 - 12.0_rp * x1020 + 12.0_rp * x1026 - x2 * x261 - x2 * x263 + x2 * x268 + 177.0_rp * x2 * x350 + 72.0_rp * x2 * x478 - 72.0_rp&
& * x2 * x480 + x2 * x488 + 7.0_rp * x2 * x815 + x228 - x24 * x44 + x322 * x45 - x421 + x422 + 15.0_rp * x43 - x45 * x867 + x47 -&
& x528 + x529 - x538 - x875) + x863 * (-162.0_rp * x1019 - 72.0_rp * x1020 - x1022 * x223 + 105.0_rp * x1022 + 72.0_rp * x1026 + x2&
&30 + x24 * x537 - 90.0_rp * x24 * x720 - 4.0_rp * x24 * x747 - 6.0_rp * x538 + 195.0_rp * x720 + x736 - x738 + x739 - x740 + x745 &
&+ 22.0_rp * x747) + x865 * (-162.0_rp * x1019 - 72.0_rp * x1020 - x1022 * x223 + 105.0_rp * x1022 + 72.0_rp * x1026 + x230 + x24 * &
&x537 - 90.0_rp * x24 * x720 - 4.0_rp * x24 * x747 - 6.0_rp * x538 + 195.0_rp * x720 + x736 - x738 + x739 - x740 + x745 + 22.0_rp * &
&x747)
            T(i, 25) = -x414 * (85.5_rp * x10 * x2 * x331 * x766 - 22.5_rp * x1013 * x165 - x1013 * x575 + x1013 * x633 + x1013 &
&* x688 + x1013 * x707 + x1013 * x709 - 67.5_rp * x1015 * x165 - x1015 * x172 * x462 + 67.5_rp * x1015 * x342 * x416 + x1015 * x4&
&53 * x462 + x1015 * x477 * x708 + x1015 * x688 - x1040 * x15 * x2 * x514 - x1040 * x15 * x2 * x522 + x1040 * x462 + x1040 * x4&
&75 - 67.5_rp * x1040 * x657 - x1045 * x2 * x475 - x1045 * x332 * x444 + x1045 * x477 + x1045 * x494 - x1045 * x983 - x1060 * x4&
&62 - x1060 * x475 + 7.5_rp * x15 * x179 * x2 * x526 + x2 * x331 * x332 * x548 + 27.0_rp * x270 * x651 + 18.0_rp * x331 * x416 * x&
&502 + 270.0_rp * x342 * x856 - x413 * (19.0_rp * x138 * x849 - x2 * x24 * x45 * x804 - x24 * x408 * x849 + x24 * x701 - 192.0_rp &
&* x24 * x850 + x261 + x263 - x268 - 177.0_rp * x350 + x408 * x849 - 153.0_rp * x445 + x45 * x529 + 192.0_rp * x45 * x856 - 72.0_rp&
& * x478 + 72.0_rp * x480 - x488 + x489 + x490 - 33.0_rp * x624 + 48.0_rp * x628 - x701 - 7.0_rp * x815 + 507.0_rp * x850) + 157.5_rp&
& * x427 * x510 + x436 * x502 * x526 + x460 * x849 - x474 * (-x1 * x317 - x1 * x318 + x1 * x410 + 11.0_rp * x1067 - 16.0_rp * x10&
&69 + x2 * x24 * x248 - x2 * x24 * x804 - x2 * x248 - 5.0_rp * x2 * x68 + x2 * x804 + x49 + x595 + 51.0_rp * x871) - x475 * x526 &
&- x477 * (-27.0_rp * x1019 - 12.0_rp * x1020 + 12.0_rp * x1026 - x2 * x261 - x2 * x263 + x2 * x268 + 177.0_rp * x2 * x350 + 72.0_rp&
& * x2 * x478 - 72.0_rp * x2 * x480 + x2 * x488 + 7.0_rp * x2 * x815 + x228 - x24 * x44 + x322 * x45 - x421 + x422 + 15.0_rp * x43&
& - x45 * x867 + x47 - x528 + x529 - x538 - x875) - x494 * (-27.0_rp * x1019 - 12.0_rp * x1020 + 12.0_rp * x1026 - x2 * x261 - x2 &
&* x263 + x2 * x268 + 177.0_rp * x2 * x350 + 72.0_rp * x2 * x478 - 72.0_rp * x2 * x480 + x2 * x488 + 7.0_rp * x2 * x815 + x228 - x2&
&4 * x44 + x322 * x45 - x421 + x422 + 15.0_rp * x43 - x45 * x867 + x47 - x528 + x529 - x538 - x875) - 630.0_rp * x495 - x514 * x5&
&58 - x518 * x526 - x522 * x558 + 7.5_rp * x526 * x657 + x526 * x674 - x526 * x807 + x567 * x849 - 27.0_rp * x569 - 315.0_rp * x57&
&1 - 67.5_rp * x605 + 45.0_rp * x984)
            T(i, 26) = -ai**5 * x155 * x29 * x331 * x907 - 52.5_rp * x10 * x558 * x909 - x1015 * x134 * x793 * x907 - x1015 * x&
&2 * x363 * x706 - x1015 * x233 * x929 + x1015 * x25 * x623 - x1015 * x27 - 472.5_rp * x1015 * x299 * x907 - x1015 * x421 * x649&
& + x1015 * x443 * x604 + 90.0_rp * x1015 * x446 - x1040 * x15 * x2 * x753 + 30.0_rp * x1040 * x15 * x21 * x909 + x1040 * x15 * x&
&771 * x907 + x1040 * x15 * x772 * x907 + 2.0_rp * x1040 * x15 * x923 * x94 - x1040 * x187 * x455 - 6.0_rp * x1040 * x429 * x502 &
&- x1060 * x128 - x1060 * x15 * x712 + x1060 * x2 * x46 * x502 + x1060 * x2 * x604 + x1060 * x2 * x769 - x1081 * x790 - x1081 *&
& x802 + 990.0_rp * x17 * x850 + 45.0_rp * x179 * x2 * x558 + 1890.0_rp * x2 * x495 - 12.0_rp * x2 * x552 * x706 + 180.0_rp * x2 * x&
&569 - x2 * x651 * x738 - 22.5_rp * x23 * x331 * x502 * x907 - 420.0_rp * x292 * x342 * x907 - 90.0_rp * x293 + x331 * x455 * x502&
& - x331 * x475 - x331 * x518 - x331 * x807 - 45.0_rp * x351 - x373 * x708 * x849 + x374 * x651 + x413 * (x1 * x317 + x1 * x318 &
&- x1 * x410 - 66.0_rp * x1067 + 96.0_rp * x1069 - 128.0_rp * x1081 * x134 * x24 + 443.0_rp * x1081 * x134 - x1081 * x24 * x408 + x&
&1081 * x408 + 210.0_rp * x24 * x856 - x24 * x904 + 35.0_rp * x24 * x911 + 93.0_rp * x292 * x909 + x391 - x49 - x51 + x663 - 210.0&
&_rp * x856 - 306.0_rp * x871 + x876 + x904 - 35.0_rp * x911 - 35.0_rp * x913) + x421 * x510 + x499 * x558 * x859 + x506 * x649 - 5&
&2.5_rp * x510 * x911 - 5.0_rp * x548 * x558 * x907 - 180.0_rp * x649 * x856 - x764 * x922 + 6.0_rp * x863 * (-x1 * x317 - x1 * x31&
&8 + x1 * x410 + 11.0_rp * x1067 - 16.0_rp * x1069 + x2 * x24 * x248 - x2 * x24 * x804 - x2 * x248 - 5.0_rp * x2 * x68 + x2 * x804&
& + x49 + x595 + 51.0_rp * x871) + 6.0_rp * x865 * (-x1 * x317 - x1 * x318 + x1 * x410 + 11.0_rp * x1067 - 16.0_rp * x1069 + x2 * x&
&24 * x248 - x2 * x24 * x804 - x2 * x248 - 5.0_rp * x2 * x68 + x2 * x804 + x49 + x595 + 51.0_rp * x871)
            T(i, 27) = -x414 * (-x115 * x651 + x157 * (x130 * x24 * x907 - x130 * x907 - x135 * x907 - x136 * x911 - x141 * x9&
&07 + x143 * x24 * x907 - 6.0_rp * x143 * x907 + x146 * x907 * x94 + x149 * x24 * x907 - x150 + x153 * x24 * x907 + x24 * x253 -&
& 1050.0_rp * x24 * x427 - 450.0_rp * x24 * x624 - 100.0_rp * x24 * x898 + 15.0_rp * x24 * x97 - x241 * x872 - x253 + 225.0_rp * x25&
&6 + 1050.0_rp * x427 + 2100.0_rp * x445 + 1350.0_rp * x624 - 1050.0_rp * x628 - 450.0_rp * x68 + 90.0_rp * x79 + 50.0_rp * x895 + 400&
&.0_rp * x898 - x907 * x942 + x911 * x946 - 60.0_rp * x97) + 225.0_rp * x165 + 787.5_rp * x168 + x170 * x929 + x186 * x2 * x708 - x&
&186 * x706 - x2 * x253 * x342 - x210 * x708 + x233 * (x108 * x909 + x125 - 150.0_rp * x2 * x256 + 300.0_rp * x2 * x68 - 60.0_rp *&
& x2 * x79 + x2 * x99 - 135.0_rp * x23 * x909 - 150.0_rp * x24 * x416 + 105.0_rp * x24 * x911 + x24 * x918 + 150.0_rp * x416 + 40.0&
&_rp * x665 - 10.0_rp * x667 - x744 * x907 + x748 * x907 * x92 + x749 * x913 - 105.0_rp * x911 - 210.0_rp * x913 - 5.0_rp * x918) - &
&x253 * x649 + 262.5_rp * x342 * x911 + 75.0_rp * x416 * x651 + x42 * x843 + 525.0_rp * x427 * x649 - 4725.0_rp * x632 + 75.0_rp * x&
&657 * x706 - 2100.0_rp * x687 + 2.5_rp * x86 * x922 + x907 * x961 + x907 * x963)
            T(i, 28) = -x128 * x2 * (x108 * x909 + x125 - 150.0_rp * x2 * x256 + 300.0_rp * x2 * x68 - 60.0_rp * x2 * x79 + x2 * &
&x99 - 135.0_rp * x23 * x909 - 150.0_rp * x24 * x416 + 105.0_rp * x24 * x911 + x24 * x918 + 150.0_rp * x416 + 40.0_rp * x665 - 10.0_rp * x667 - x744 * x907 + x748 * x907 * x92 + x749 * x913 - 105.0_rp * x911 - 210.0_rp * x913 - 5.0_rp * x918) - x157 * (-x102 * x&
&907 * x92 + 420.0_rp * x1087 * x110 * x151 + x1087 * x130 * x24 - x1087 * x130 - x1087 * x135 - x1087 * x137 - x1087 * x141 + x&
&1087 * x143 * x24 - 6.0_rp * x1087 * x143 + x1087 * x147 + x1087 * x149 * x24 + x1087 * x153 * x24 - x1087 * x942 - x145 * x23 &
&* x909 + x150 * x872 + x154 + 675.0_rp * x2 * x256 - 1350.0_rp * x2 * x68 - 810.0_rp * x2 * x78 + 2025.0_rp * x23 * x909 + 675.0_rp&
& * x24 * x416 - 1575.0_rp * x24 * x911 - 1575.0_rp * x24 * x913 - 15.0_rp * x24 * x918 - 675.0_rp * x416 - 180.0_rp * x665 + 45.0_rp&
& * x667 + 600.0_rp * x907 * x92 + 1575.0_rp * x911 + 3150.0_rp * x913 + 75.0_rp * x918) - 30.0_rp * x2 * x42 * x706 * (5.0_rp * x430&
& - 3.0_rp) - x27 * (x1087 * x22 + 105.0_rp * x430 - 315.0_rp * x923 - 5.0_rp) - 22.5_rp * x42 * x443 * (-30.0_rp * x430 + 35.0_rp * x&
&923 + 3.0_rp) - 45.0_rp * x446 * (-70.0_rp * x430 + 63.0_rp * x923 + 15.0_rp) - x87 * x929 * (3.0_rp * x430 - 1.0_rp)
        end do
    end subroutine T6_damp_thole
end module T_tensor_damp_thole
