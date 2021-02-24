module tensors_recursive
    implicit none
    private

    public Tn_recursive
    real(8), dimension(:, :, :), allocatable, save :: Cnij
    integer, save :: current_Cnij_order = 0
contains
    function Txyz(kx, ky, kz, x_powers, y_powers, z_powers, inv_R_power) result(T)
        integer, intent(in) :: kx, ky, kz
        real(8), intent(in) :: x_powers(:, :), y_powers(:, :), z_powers(:, :), inv_R_power(:)
        real(8), allocatable :: T(:)
        integer :: a, b, c
        real(8), allocatable :: Cx(:), Cy(:), Cz(:)
        real(8) :: Cnij_factor

        allocate (T(size(inv_R_power)))
        T(:) = 0.0

        allocate (Cx(size(inv_R_power)))
        allocate (Cy(size(inv_R_power)))
        allocate (Cz(size(inv_R_power)))

        do a = 0, kx
            Cnij_factor = Cnij(1, kx, a)
            if (Cnij_factor == 0) cycle
            Cx = Cnij_factor*x_powers(:, a + 1)
            do b = 0, ky
                Cnij_factor = Cnij(a + kx + 1, ky, b)
                if (Cnij_factor == 0) cycle
                Cy = Cx*Cnij_factor*y_powers(:, b + 1)
                do c = 0, kz
                    Cnij_factor = Cnij(a + kx + b + ky + 1, kz, c)
                    if (Cnij_factor == 0) cycle
                    Cz = Cy*Cnij_factor*z_powers(:, c + 1)
                    T = T + Cz
                end do
            end do
        end do

        T = T*inv_R_power
    end function Txyz

    subroutine set_Cnij_coefficients(max_order)
        integer, intent(in) :: max_order
        integer :: i, j, k, n
        if (allocated(Cnij)) deallocate (Cnij)
        allocate (Cnij(2*max_order + 3, 0:max_order + 1, 0:max_order + 1))
        Cnij = 0
        Cnij(:, 0, 0) = 1
        do n = 1, 2*max_order + 3
            if (mod(n, 2) == 0) cycle
            do i = 1, max_order + 1
                if (mod(i, 2) /= 0) then
                    k = i - 1
                else if (mod(i, 2) == 0) then
                    k = i
                end if
                do j = 0, i
                    if (mod(i + j, 2) /= 0) cycle
                    if (j == 0) then
                        Cnij(n, i, j) = Cnij(n, i - 1, j + 1)
                    else if (j /= i) then
                        Cnij(n, i, j) = (j + 1)*Cnij(n, i - 1, j + 1)
                        Cnij(n, i, j) = Cnij(n, i, j) - (n + k)*Cnij(n, i - 1, j - 1)
                        k = k + 2
                    else if (j == i) then
                        Cnij(n, i, j) = -(n + k)*Cnij(n, i - 1, j - 1)
                    end if
                end do
            end do
        end do
    end subroutine set_Cnij_coefficients

    subroutine Tn_recursive(n, x, y, z, T)
        integer, intent(in) :: n
        real(8), intent(in) :: x(:), y(size(x)), z(size(x))
        real(8), intent(inout) :: T(size(x), (n + 1)*(n + 2)*(n + 3)/6)
        integer :: i, kx, ky, kz
        real(8), allocatable :: x_powers(:, :), y_powers(:, :), z_powers(:, :) ! powers (x/R)**k
        real(8), allocatable :: R(:), inv_R_power(:)                  ! R, 1/R**n

        if (current_Cnij_order < n) then
            call set_Cnij_coefficients(n)
            current_Cnij_order = n
        end if

        allocate (R(size(x)))
        allocate (inv_R_power(size(x)))
        allocate (x_powers(size(x), 1:n + 1))
        allocate (y_powers(size(x), 1:n + 1))
        allocate (z_powers(size(x), 1:n + 1))

        R = sqrt(x**2 + y**2 + z**2)
        inv_R_power = 1/R**(n + 1)
        x_powers(:, 1) = 1.0 ! (x/R)**0 = 1.0
        y_powers(:, 1) = 1.0 ! (y/R)**0 = 1.0
        z_powers(:, 1) = 1.0 ! (z/R)**0 = 1.0
        do i = 2, n + 1
            x_powers(:, i) = x_powers(:, i - 1)*(x/R)
            y_powers(:, i) = y_powers(:, i - 1)*(y/R)
            z_powers(:, i) = z_powers(:, i - 1)*(z/R)
        end do

        i = 1
        do kx = n, 0, -1
            do ky = n - kx, 0, -1
                kz = n - kx - ky
                T(:, i) = Txyz(kx, ky, kz, x_powers, y_powers, z_powers, inv_R_power)
                i = i + 1
            end do
        end do
    end subroutine Tn_recursive
end module tensors_recursive
