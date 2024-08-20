module fmm
    use precision
    use t_tensor
    use t_tensor_damp_thole
    use t_tensor_damp_amoeba
    use t_tensor_damp_erf
    implicit none
    private

    public field_direct
    public field_fmm
    public get_multipoles, get_boxes, interpolate_electron_fields

    type settings_type
        integer(ip) :: ncrit
        real(rp) :: theta
        integer(ip) :: expansion_order
        integer(ip) :: multipole_size
        logical :: is_periodic
        integer :: num_interaction_lists ! non-periodic: 1, periodic: 27
        integer, allocatable, dimension(:,:) :: aijk
        real(rp), dimension(3,3) :: box_matrix
    end type settings_type

    type list_type
        integer(ip) :: head = 1
        integer(ip), allocatable :: elements(:)
    end type list_type

    type node_type
        ! tree things
        integer(ip) :: idx
        integer(ip) :: depth
        integer(ip) :: nleaf
        integer(ip), allocatable :: leaf(:)
        real(rp) :: center(3)
        real(rp) :: r
        real(rp) :: rmax
        type(node_type), pointer :: parent
        type(node_type), pointer :: children(:)
        logical(lp) :: occupied(0:7) = .FALSE.
        logical(lp) :: is_terminal = .TRUE.
        ! interaction list - "list-of-lists" or "ragged array" 
        ! for nonperiodic: just a list (with the 000 cell)
        ! for periodic: 27 elements in 3x3x3 supercell:
        !    1   2   3   4   5   6   7   8   9  10  11  12  13  14
        ! (---,--0,--+,-0-,-00,-0+,-+-,-+0,-++,0--,0-0,0-+,00-,000,00+,0+-,0+0,0++,+--,+-0,+-+,+0-,+00,+0+,++-,++0,+++)
        type(list_type), allocatable :: particle_interactions(:)
        type(list_type), allocatable :: cell_interactions(:)
    end type node_type
    !             parent
    !               |
    !              node
    !             /    \
    !     child(0) ... child(7)
    !                   | | | |
    !                   1 2 3 4 (particles)

    type node_list_type
        type(node_type), pointer :: node
    end type node_list_type

    type tree_type
        real(rp), pointer :: coordinates(:, :)
        real(rp), pointer :: source_multipoles(:, :)
        type(node_type), pointer :: root_node
        integer(ip) :: num_nodes
        integer(ip) :: max_depth
        type(node_list_type), allocatable :: node_list(:)
        real(rp), allocatable :: centers(:, :)
        real(rp), allocatable :: cell_multipoles(:, :)
        real(rp), allocatable :: local_expansion(:, :)
    end type tree_type

    type(settings_type) settings
    type(tree_type) tree

contains

    pure function xyz2idx(x, y, z) result(idx)
        implicit none
        integer(ip), intent(in) :: x, y, z
        integer(ip) :: k
        integer(ip) :: idx
        k = x + y + z
        ! number of components before order k is k*(k+1)*(k+2)/6
        ! (y**2 + 2*y*z + y + z**2 + 3*z)/2 + 1 is the symmetry-packed index of the current slice
        idx = k * (k + 1) * (k + 2) / 6 + (y**2 + 2 * y * z + y + z**2 + 3 * z) / 2 + 1
    end function xyz2idx

    pure function factorial(n)
        integer(ip), intent(in) :: n
        integer(ip) :: i
        real(rp) :: factorial
        factorial = 1
        do i = n, 1, -1
            factorial = factorial * real(i, rp)
        end do
    end function

    pure function trinom(i, j, k)
        integer(ip), intent(in) :: i, j, k
        real(rp) :: trinom
        trinom = factorial(i + j + k) / (factorial(i) * factorial(j) * factorial(k))
    end function trinom

    pure function binom(n, k)
        integer(ip), intent(in) :: n, k
        real(rp) :: binom
        binom = factorial(n) / (factorial(k) * factorial(n - k))
    end function binom

    subroutine list_append(list, element)
        type(list_type), intent(inout) :: list
        integer(ip), intent(in) :: element
        integer(ip), allocatable :: swap(:)
        if (list%head > size(list%elements)) then
            allocate (swap(list%head))
            call move_alloc(list%elements, swap)
            allocate (list%elements(2 * list%head))
            list%elements(1:size(swap)) = swap(:)
            deallocate (swap)
        end if
        list%elements(list%head) = element
        list%head = list%head + 1
    end subroutine list_append

    subroutine list_trim(list)
        type(list_type), intent(inout) :: list
        integer(ip), allocatable :: swap(:)
        call move_alloc(list%elements, swap)
        allocate (list%elements(1:list%head - 1))
        list%elements(1:list%head - 1) = swap(1:list%head - 1)
        deallocate (swap)
    end subroutine list_trim

    subroutine octree_build(ncrit, expansion_order, theta, coordinates, multipoles)
        integer(ip), intent(in) :: ncrit
        integer(ip), intent(in) :: expansion_order
        real(rp), intent(in) :: theta
        real(rp), intent(in), pointer :: coordinates(:, :)
        real(rp), intent(in), pointer :: multipoles(:, :)
        integer(ip) :: i, ai, aj, ak, idx
        integer(ip) :: octant
        integer(ip) :: multipole_size
        integer(ip) :: n
        type(node_type), pointer :: node

        settings%ncrit = ncrit
        settings%expansion_order = expansion_order
        settings%theta = theta
        n = settings%expansion_order
        multipole_size = (n + 1) * (n + 2) * (n + 3) / 6
        settings%multipole_size = multipole_size
        ! periodic data
        if (settings%is_periodic) then 
            settings%num_interaction_lists = 27
            allocate(settings%aijk(27, 3)) 
            idx = 1
            do ai=-1,1
                do aj=-1,1
                    do ak=-1,1
                        settings%aijk(idx, 1) = ai
                        settings%aijk(idx, 2) = aj
                        settings%aijk(idx, 3) = ak
                        idx = idx + 1
                    end do
                end do
            end do
        else
            settings%num_interaction_lists = 1
            allocate(settings%aijk(1, 3))
            settings%box_matrix = 0.0_rp
            settings%aijk(1, :) = [0, 0, 0]
        end if

        if (.not. associated(tree%root_node)) allocate (tree%root_node)
        tree%root_node%nleaf = 0
        allocate (tree%root_node%leaf(settings%ncrit))
        allocate (tree%root_node%children(0:7))
        nullify (tree%root_node%parent)
        tree%num_nodes = 1
        tree%root_node%depth = 1
        tree%root_node%center = sum(coordinates, 1) / real(size(coordinates, 1), rp)
        tree%root_node%r = max( &
                           maxval(abs(coordinates(:, 1) - tree%root_node%center(1))), &
                           maxval(abs(coordinates(:, 2) - tree%root_node%center(2))), &
                           maxval(abs(coordinates(:, 3) - tree%root_node%center(3)))) * 1.0001_rp
        tree%root_node%rmax = sqrt(0.5 * 3 * tree%root_node%r * tree%root_node%r)
        tree%coordinates => coordinates
        tree%source_multipoles => multipoles

        do i = 1, size(coordinates, 1)
            node => tree%root_node
            do while (node%nleaf >= settings%ncrit)
                node%nleaf = node%nleaf + 1
                octant = 0
                if (coordinates(i, 1) > node%center(1)) then
                    octant = octant + 1
                end if
                if (coordinates(i, 2) > node%center(2)) then
                    octant = octant + 2
                end if
                if (coordinates(i, 3) > node%center(3)) then
                    octant = octant + 4
                end if
                if (.not. node%occupied(octant)) then
                    call add_child(node, octant)
                end if
                node => node%children(octant)
            end do
            node%nleaf = node%nleaf + 1
            node%leaf(node%nleaf) = i
            if (node%nleaf >= settings%ncrit) then
                call split_node(node)
            end if
        end do
        call create_node_list()
    end subroutine octree_build

    subroutine add_child(node, octant)
        type(node_type), intent(inout), target :: node
        integer(ip), intent(in) :: octant
        real(rp) :: r
        real(rp) :: center(3)

        tree%num_nodes = tree%num_nodes + 1
        node%occupied(octant) = .TRUE.
        node%is_terminal = .FALSE.

        r = node%r * 0.5
        center = node%center
        center(1) = center(1) + r * real((IAND(octant, 1) * 2 - 1), rp)
        center(2) = center(2) + r * real((IAND(octant, 2) - 1), rp)
        center(3) = center(3) + r * real((IAND(octant, 4) / 2 - 1), rp)

        node%children(octant)%r = r
        node%children(octant)%rmax = sqrt(0.5 * 3 * r * r)
        node%children(octant)%center = center
        node%children(octant)%parent => node
        node%children(octant)%depth = node%depth + 1
        node%children(octant)%nleaf = 0
        allocate (node%children(octant)%children(0:7))
        allocate (node%children(octant)%leaf(settings%ncrit))
    end subroutine add_child

    recursive subroutine split_node(node)
        type(node_type), intent(inout), target :: node
        integer(ip) :: i
        integer(ip) :: leaf
        integer(ip) :: octant
        do i = 1, node%nleaf
            leaf = node%leaf(i)
            octant = 0
            if (tree%coordinates(leaf, 1) > node%center(1)) then
                octant = octant + 1
            end if
            if (tree%coordinates(leaf, 2) > node%center(2)) then
                octant = octant + 2
            end if
            if (tree%coordinates(leaf, 3) > node%center(3)) then
                octant = octant + 4
            end if
            if (.not. node%occupied(octant)) then
                call add_child(node, octant)
            end if
            node%children(octant)%nleaf = node%children(octant)%nleaf + 1
            node%children(octant)%leaf(node%children(octant)%nleaf) = leaf
            if (node%children(octant)%nleaf >= settings%ncrit) then
                call split_node(node%children(octant))
            end if
        end do
    end subroutine split_node

    subroutine create_node_list
        integer(ip) :: pos, i
        allocate (tree%node_list(tree%num_nodes))
        pos = 1
        call downward_pass(tree%root_node, pos)
        allocate (tree%centers(tree%num_nodes, 3))
        do i = 1, tree%num_nodes
            tree%centers(i, :) = tree%node_list(i)%node%center(:)
        end do
    end subroutine create_node_list

    recursive subroutine downward_pass(node, pos)
        type(node_type), intent(in), target :: node
        integer(ip), intent(inout) :: pos
        integer(ip) :: octant
        tree%max_depth = max(tree%max_depth, node%depth)
        tree%node_list(pos)%node => node
        tree%node_list(pos)%node%idx = pos
        pos = pos + 1
        do octant = 0, 7
            if (node%occupied(octant)) then
                call downward_pass(node%children(octant), pos)
            end if
        end do
    end subroutine downward_pass

    subroutine finalize
        call clean_node(tree%root_node)
        if (associated(tree%root_node)) deallocate (tree%root_node)
        if (allocated(tree%node_list)) deallocate (tree%node_list)
        if (allocated(tree%centers)) deallocate (tree%centers)
        if (allocated(tree%cell_multipoles)) deallocate (tree%cell_multipoles)
        if (allocated(tree%local_expansion)) deallocate (tree%local_expansion)
        if (allocated(settings%aijk)) deallocate (settings%aijk)
    end subroutine finalize

    recursive subroutine clean_node(node)
        type(node_type), intent(inout) :: node
        integer(ip) octant
        do octant = 0, 7
            if (node%occupied(octant)) then
                call clean_node(node%children(octant))
            end if
        end do
        deallocate (node%children)
        deallocate (node%particle_interactions)
        deallocate (node%cell_interactions)
        deallocate (node%leaf)
    end subroutine clean_node

    subroutine build_interaction_lists_fmm
        integer(ip) :: i, j, k
        k = 1
        do i = 1, tree%num_nodes
            ! just temporary allocation, will be resized to fit later
            allocate (tree%node_list(i)%node%particle_interactions(settings%num_interaction_lists))
            allocate (tree%node_list(i)%node%cell_interactions(settings%num_interaction_lists))
            do j=1, settings%num_interaction_lists
                allocate (tree%node_list(i)%node%particle_interactions(j)%elements(64))
                allocate (tree%node_list(i)%node%cell_interactions(j)%elements(64))
            end do
        end do
        if (settings%is_periodic) then
            do k=1, settings%num_interaction_lists
                call interact_fmm(tree%root_node, tree%root_node, settings%theta, k)
            end do
        else
            call interact_fmm(tree%root_node, tree%root_node, settings%theta, 1)
        end if
        do i = 1, tree%num_nodes
            do j=1, settings%num_interaction_lists
                call list_trim(tree%node_list(i)%node%particle_interactions(j))
                call list_trim(tree%node_list(i)%node%cell_interactions(j))
            end do
        end do
    end subroutine build_interaction_lists_fmm

    recursive subroutine interact_fmm(node_i, node_j, theta, k)
        type(node_type) :: node_i, node_j
        integer, intent(in) :: k
        integer(ip) :: octant
        real(rp) :: theta
        real(rp) :: delta(3), displacement(3)
        real(rp) :: r

        displacement = matmul(settings%box_matrix, real(settings%aijk(k, :), rp))

        delta = node_i%center - (node_j%center + displacement)
        r = norm2(delta)
        if (r * theta > node_i%rmax + node_j%rmax) then
            ! far field
            call list_append(node_i%cell_interactions(k), node_j%idx)
        else if (node_i%is_terminal .and. node_j%is_terminal) then
            ! both are leaf nodes, and should interact directly
            call list_append(node_i%particle_interactions(k), node_j%idx)
        else if (node_j%is_terminal .or. ((node_i%rmax > node_j%rmax) .and. (.not. node_i%is_terminal))) then
            ! too big, descend i
            do octant = 0, 7
                if (node_i%occupied(octant)) then
                    call interact_fmm(node_i%children(octant), node_j, theta, k)
                end if
            end do
        else
            ! too big, descend j
            do octant = 0, 7
                if (node_j%occupied(octant)) then
                    call interact_fmm(node_i, node_j%children(octant), theta, k)
                end if
            end do
        end if
    end subroutine interact_fmm

    subroutine multipole_accumulate_kernel(delta, source_multipoles, target_multipoles)
        real(rp), intent(in) :: delta(3)
        real(rp), dimension(:), intent(in) :: source_multipoles
        real(rp), dimension(:), intent(inout) :: target_multipoles
        integer(ip) :: source_index, source_order, sx, sy, sz
        integer(ip) :: target_index, target_order, tx, ty, tz
        real(rp) :: symfac
        integer(ip) :: max_multipole_order
        real(rp) :: monomial
        max_multipole_order = NINT((real(size(tree%cell_multipoles, 2), rp) * 6_rp)**(1.0_rp/3.0_rp)) - 2_ip
        do source_order = 0, max_multipole_order
            do target_order = source_order, settings%expansion_order
                do sx = source_order, 0, -1
                    do sy = source_order, 0, -1
                        do sz = source_order, 0, -1
                            if (sx + sy + sz /= source_order) cycle
                            source_index = xyz2idx(sx, sy, sz)
                            do tx = target_order, sx, -1
                                do ty = target_order, sy, -1
                                    do tz = target_order, sz, -1
                                        if (tx + ty + tz /= target_order) cycle
                                        target_index = xyz2idx(tx, ty, tz)
                                        ! unsure about this ....
                                        ! maybe needs symmetry factor?
                                        symfac = binom(tx, sx) * binom(ty, sy) * binom(tz, sz)
                                        monomial = delta(1)**(tx - sx) * delta(2)**(ty - sy) * delta(3)**(tz - sz)
                                        target_multipoles(target_index) = target_multipoles(target_index) &
                                                                        + symfac * monomial * source_multipoles(source_index)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine multipole_accumulate_kernel

    recursive subroutine multipole_accumulate(source_node_index, target_node_index)
        integer, intent(in) :: source_node_index, target_node_index
        real(rp) :: delta(3)

        if (tree%node_list(target_node_index)%node%idx /= 1) then
            call multipole_accumulate(source_node_index, tree%node_list(target_node_index)%node%parent%idx)
        end if
        delta = tree%centers(source_node_index, :) - tree%centers(target_node_index, :)
        call multipole_accumulate_kernel(delta, tree%cell_multipoles(source_node_index, :), &
                                                tree%cell_multipoles(target_node_index, :))
    end subroutine multipole_accumulate

    subroutine shift_multipole_vec(source_multipole, target_multipole, dx, dy, dz)
        real(rp), intent(in) :: source_multipole(:, :)
        real(rp), intent(inout) :: target_multipole(:)
        real(rp), intent(in) :: dx(:), dy(:), dz(:)
        integer(ip) :: source_index, source_order, sx, sy, sz
        integer(ip) :: target_index, target_order, tx, ty, tz
        real(rp) :: symfac
        integer(ip) :: max_multipole_order
        real(rp) :: monomial(size(dx))

        max_multipole_order = NINT((real(size(source_multipole, 2), rp) * 6_rp)**(1.0_rp/3.0_rp)) - 2_ip
        do source_order = 0, max_multipole_order
            do target_order = source_order, settings%expansion_order
                do sx = source_order, 0, -1
                    do sy = source_order, 0, -1
                        do sz = source_order, 0, -1
                            if (sx + sy + sz /= source_order) cycle
                            source_index = xyz2idx(sx, sy, sz)
                            do tx = target_order, sx, -1
                                do ty = target_order, sy, -1
                                    do tz = target_order, sz, -1
                                        if (tx + ty + tz /= target_order) cycle
                                        target_index = xyz2idx(tx, ty, tz)
                                        symfac = binom(tx, sx) * binom(ty, sy) * binom(tz, sz)
                                        monomial = dx**(tx - sx) * dy**(ty - sy) * dz**(tz - sz)
                                        target_multipole(target_index) = target_multipole(target_index) &
                                                                         + symfac * sum(monomial * source_multipole(:, source_index))
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine shift_multipole_vec

    subroutine multipole_expansion(comm)
        ! particle_to_multipole and multipole_to_multipole combined
#ifdef VAR_MPI
#if defined(USE_MPI_MOD_F90)
        use mpi
#else
        include 'mpif.h'
#endif
#endif
        integer(ip), intent(in) :: comm
        integer(ip) :: i, num_terminal
        real(rp), allocatable :: dx(:), dy(:), dz(:)
        type(node_type) :: node
        integer(ip) :: mpi_size, mpi_rank, ierr, work_size, work_start, work_stop, node_idx
        integer(ip), allocatable :: terminal_node_list(:)

        allocate (tree%cell_multipoles(tree%num_nodes, settings%multipole_size))
        tree%cell_multipoles = 0.0_rp
#ifdef VAR_MPI
        call mpi_comm_size(comm, mpi_size, ierr)
        call mpi_comm_rank(comm, mpi_rank, ierr)
#else
        mpi_size = 1
        mpi_rank = 0
#endif

        num_terminal = 0
        do i = 1, tree%num_nodes
            node = tree%node_list(i)%node
            if (.not. node%is_terminal) cycle
            num_terminal = num_terminal + 1
        end do
        allocate (terminal_node_list(num_terminal))
        num_terminal = 0
        do i = 1, tree%num_nodes
            node = tree%node_list(i)%node
            if (.not. node%is_terminal) cycle
            num_terminal = num_terminal + 1
            terminal_node_list(num_terminal) = node%idx
        end do

        work_size = num_terminal / mpi_size
        work_start = 1 + mpi_rank * work_size
        work_stop = (mpi_rank + 1) * work_size
        if (mpi_rank == mpi_size - 1) work_stop = num_terminal

        ! particle_to_multipole
        do i = work_start, work_stop
            node_idx = terminal_node_list(i)
            node = tree%node_list(node_idx)%node
            allocate (dx(node%nleaf))
            allocate (dy(node%nleaf))
            allocate (dz(node%nleaf))
            dx(:) = tree%coordinates(node%leaf(1:node%nleaf), 1) - node%center(1)
            dy(:) = tree%coordinates(node%leaf(1:node%nleaf), 2) - node%center(2)
            dz(:) = tree%coordinates(node%leaf(1:node%nleaf), 3) - node%center(3)
            call shift_multipole_vec(tree%source_multipoles(node%leaf(1:node%nleaf), :), tree%cell_multipoles(node_idx, :), dx, dy, dz)
            deallocate (dx)
            deallocate (dy)
            deallocate (dz)
        end do

        do i = work_start, work_stop
            node_idx = terminal_node_list(i)
            node = tree%node_list(node_idx)%node
            call multipole_accumulate(node%idx, node%parent%idx)
        end do

        deallocate (terminal_node_list)
#ifdef VAR_MPI
        call mpi_allreduce(MPI_IN_PLACE, tree%cell_multipoles(1, 1), size(tree%cell_multipoles), MPI_REAL8, MPI_SUM, comm, ierr)
#endif
    end subroutine multipole_expansion

    function multipole_field(multipole, x, y, z, max_field_order, damp_type, damping_factors) result(F)
        ! field up to and including field_order
        ! compute field from N multipoles on one point
        real(rp), allocatable :: F(:)
        real(rp), intent(in) :: multipole(:, :)
        real(rp), intent(in) :: x(:), y(:), z(:)
        integer(ip), intent(in) :: max_field_order
        character(len=*), intent(in) :: damp_type
        real(rp), intent(in) :: damping_factors(:)
        integer(ip) :: max_multipole_order, multipole_order, sx, sy, sz
        integer(ip) :: field_order, tx, ty, tz
        integer(ip) :: multipole_index, field_index, tensor_index
        integer(ip) :: max_tensor_order
        real(rp) :: symfac
        real(rp) :: taylor
        integer(ip) :: length, N
        real(rp), allocatable :: T(:, :)

        N = size(x)

        length = (max_field_order + 1) * (max_field_order + 2) * (max_field_order + 3) / 6
        allocate (F(length))
        F = 0.0_rp
        max_multipole_order = NINT((real(size(multipole, 2), rp) * 6_rp)**(1.0_rp/3.0_rp)) - 2_ip

        max_tensor_order = max_field_order + max_multipole_order

        length = (max_tensor_order + 1) * (max_tensor_order + 2) * (max_tensor_order + 3) / 6
        allocate (T(N, length))
        T = 0.0_rp
        ! T(i, 1) -> T0
        ! T(i, 2) -> T1_x
        ! T(i, 3) -> T1_y
        ! T(i, 4) -> T1_z
        ! T(i, 5) -> T2_xx
        ! T(i, 6) -> T2_xy
        ! ...
        if (trim(damp_type) == 'AMOEBA') then
            call Tn_damp_amoeba(max_tensor_order, x, y, z, damping_factors, T)
        else if (trim(damp_type) == 'THOLE') then
            call Tn_damp_thole(max_tensor_order, x, y, z, damping_factors, T)
        else if (trim(damp_type) == 'ERF') then
            call Tn_damp_erf(max_tensor_order, x, y, z, damping_factors, T)
        else if (trim(damp_type) == '') then
            call Tn(max_tensor_order, x, y, z, T)
        else
            error stop "Unknown damp type: " // damp_type
        end if
        do field_order = 0, max_field_order
            do multipole_order = 0, max_multipole_order
                taylor = (-1.0)**(multipole_order + 1) / factorial(multipole_order)
                do sx = multipole_order, 0, -1
                    do sy = multipole_order, 0, -1
                        do sz = multipole_order, 0, -1
                            if (sx + sy + sz /= multipole_order) cycle
                            multipole_index = xyz2idx(sx, sy, sz)
                            do tx = field_order, 0, -1
                                do ty = field_order, 0, -1
                                    do tz = field_order, 0, -1
                                        if (field_order + multipole_order > max_tensor_order) cycle
                                        if (tx + ty + tz /= field_order) cycle
                                        field_index = xyz2idx(tx, ty, tz)
                                        tensor_index = xyz2idx(sx + tx, sy + ty, sz + tz)
                                        symfac = trinom(sx, sy, sz)
                                        F(field_index) = F(field_index) + taylor * symfac * SUM(multipole(:, multipole_index) * T(:, tensor_index))
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end function multipole_field

    subroutine particle_fields(comm, field, exclusions, max_field_order, damp_type, damping_factors)
#ifdef VAR_MPI
#if defined(USE_MPI_MOD_F90)
        use mpi
#else
        include 'mpif.h'
#endif
#endif
        integer(ip), intent(in) :: comm
        real(rp), intent(inout) :: field(:, :)
        integer(ip), intent(in) :: exclusions(:, :)
        integer(ip), intent(in) :: max_field_order
        character(len=*), intent(in) :: damp_type
        real(rp), intent(in) :: damping_factors(:)
        real(rp), allocatable :: full_damping_factors(:)

        integer(ip) :: j, ci, cj, leaf_i, leaf_j
        real(rp)  :: dx(settings%ncrit), dy(settings%ncrit), dz(settings%ncrit)
        type(list_type) :: particle_list
        type(node_type) :: node_i, node_j
        integer(ip) :: mpi_size, mpi_rank, ierr
        integer(ip), allocatable :: work_node_list(:), work_leaf_list(:)
        integer(ip) :: work_start, work_stop, work_idx, N, total_work
        integer(ip) :: ideal_work_size, ideal_work_start, ideal_work_stop
        logical(lp) :: set_work_start, set_work_stop
        integer(ip) :: periodic_index
        logical(lp) :: central
        real(rp), dimension(3) :: displacement
#ifdef VAR_MPI
        call mpi_comm_size(comm, mpi_size, ierr)
        call mpi_comm_rank(comm, mpi_rank, ierr)
#else
        mpi_size = 1
        mpi_rank = 0
#endif
        N = size(field, 1)
        allocate (work_node_list(N))
        allocate (work_leaf_list(N))
        ! get total work estimate
        do periodic_index = 1, settings%num_interaction_lists
            displacement = matmul(settings%box_matrix, real(settings%aijk(periodic_index, :), rp))
            central = .false.
            if (all(settings%aijk(periodic_index, :) == 0)) central = .true.
            total_work = 0
            do ci = 1, tree%num_nodes
                node_i = tree%node_list(ci)%node
                if (.not. node_i%is_terminal) cycle
                total_work = total_work + node_i%nleaf * (node_i%particle_interactions(periodic_index)%head - 1)
            end do
            if (total_work == 0) cycle

            ! start/stop estimates
            ideal_work_size = total_work / mpi_size
            ideal_work_start = mpi_rank * ideal_work_size
            ideal_work_stop = (mpi_rank + 1) * ideal_work_size
            work_start = total_work

            total_work = 0
            work_idx = 1
            set_work_start = .true.
            set_work_stop = .true.
            do ci = 1, tree%num_nodes
                node_i = tree%node_list(ci)%node
                if (.not. node_i%is_terminal) cycle
                do j = 1, node_i%nleaf
                    if (total_work >= ideal_work_start .and. set_work_start) then
                        work_start = work_idx
                        set_work_start = .false.
                    end if
                    if (total_work >= ideal_work_stop .and. set_work_stop) then
                        work_stop = work_idx - 1
                        set_work_stop = .false.
                    end if
                    work_leaf_list(work_idx) = node_i%leaf(j)
                    work_node_list(work_idx) = node_i%idx
                    work_idx = work_idx + 1
                    total_work = total_work + (node_i%particle_interactions(periodic_index)%head - 1)
                end do
            end do
            if (mpi_rank == mpi_size - 1) work_stop = N

            do ci = work_start, work_stop
                node_i = tree%node_list(work_node_list(ci))%node
                particle_list = node_i%particle_interactions(periodic_index)
                leaf_i = work_leaf_list(ci)
                if (size(particle_list%elements) == 0) cycle
                do cj = 1, size(particle_list%elements)
                    node_j = tree%node_list(particle_list%elements(cj))%node
                    dx(1:node_j%nleaf) = tree%coordinates(leaf_i, 1) - (tree%coordinates(node_j%leaf(1:node_j%nleaf), 1) + displacement(1))
                    dy(1:node_j%nleaf) = tree%coordinates(leaf_i, 2) - (tree%coordinates(node_j%leaf(1:node_j%nleaf), 2) + displacement(2))
                    dz(1:node_j%nleaf) = tree%coordinates(leaf_i, 3) - (tree%coordinates(node_j%leaf(1:node_j%nleaf), 3) + displacement(3))
                    do j = 1, node_j%nleaf
                        leaf_j = node_j%leaf(j)
                        if (any(exclusions(leaf_i, :) == leaf_j) .and. central) then
                            dx(j) = 1e10_rp ! T(x,y,z) -> 0. as x,y,z -> large numbers
                            dy(j) = 1e10_rp ! a bit of a hacky way to do exclusions
                            dz(j) = 1e10_rp ! but makes it easy to pass the entire vector into multipole_field
                        end if
                    end do
                    if (damp_type /= '') then
                        allocate (full_damping_factors(node_j%nleaf))
                        do j = 1, node_j%nleaf
                            leaf_j = node_j%leaf(j)
                            full_damping_factors(j) = damping_factors(leaf_i) * damping_factors(leaf_j)
                        end do
                    else
                        allocate (full_damping_factors(0))
                    end if
                    field(leaf_i, :) = field(leaf_i, :) + multipole_field(tree%source_multipoles(node_j%leaf(1:node_j%nleaf), :), &
                                                                          dx(1:node_j%nleaf), dy(1:node_j%nleaf), dz(1:node_j%nleaf), max_field_order, &
                                                                          damp_type, full_damping_factors)
                if (allocated(full_damping_factors)) deallocate (full_damping_factors)
                end do
            end do
        end do
        deallocate (work_node_list)
        deallocate (work_leaf_list)
#ifdef VAR_MPI
        if (mpi_rank == 0) then
            call mpi_reduce(MPI_IN_PLACE, field(1, 1), size(field), MPI_REAL8, MPI_SUM, 0, comm, ierr)
        else
            call mpi_reduce(field(1, 1), field(1, 1), size(field), MPI_REAL8, MPI_SUM, 0, comm, ierr)
        end if
#endif
    end subroutine particle_fields

    subroutine particle_fields_direct(comm, field, coordinates, multipoles, exclusions, max_field_order, damp_type, damping_factors)
        !use pelib_options, only: box_matrix, box_matrix_inv, pelib_mic
#ifdef VAR_MPI
#if defined(USE_MPI_MOD_F90)
        use mpi
#else
        include 'mpif.h'
#endif
#endif
        integer(ip), intent(in) :: comm
        real(rp), intent(inout) :: field(:, :)
        real(rp), intent(in) :: coordinates(:, :)
        real(rp), intent(in) :: multipoles(:, :)
        integer(ip), intent(in) :: exclusions(:, :)
        integer(ip), intent(in) :: max_field_order
        character(len=*), intent(in) :: damp_type
        real(rp), intent(in) :: damping_factors(:)
        real(rp), allocatable :: full_damping_factors(:)
        integer(ip) :: i, j, exclusion
        integer(ip) :: N
        integer(ip) :: work_start, work_stop, work_size
        integer(ip) :: mpi_rank, mpi_size, ierr
        real(rp), allocatable :: Rji(:, :)
        real(rp) :: s_i(3)
#ifdef VAR_MPI
        call mpi_comm_size(comm, mpi_size, ierr)
        call mpi_comm_rank(comm, mpi_rank, ierr)
#else
        mpi_size = 1
        mpi_rank = 0
#endif
        N = size(coordinates, 1)
        work_size = N / mpi_size
        work_start = 1 + mpi_rank * work_size
        work_stop = (mpi_rank + 1) * work_size
        if (mpi_rank == mpi_size - 1) work_stop = N

        allocate(Rji(N,3))

        do i = work_start, work_stop
            ! Works for any triclinic box through the box matrix
            !if (pelib_mic) then
            !    s_i = matmul(box_matrix_inv, coordinates(i, :))
            !    do j = 1, N
            !        Rji(j, :) = s_i - matmul(box_matrix_inv, coordinates(j, :))
            !        Rji(j, :) = Rji(j, :) - anint(Rji(j, :))
            !        Rji(j, :) = matmul(box_matrix, Rji(j, :))
            !    end do
            !else
                Rji(:, 1) = coordinates(i, 1) - coordinates(:, 1)
                Rji(:, 2) = coordinates(i, 2) - coordinates(:, 2)
                Rji(:, 3) = coordinates(i, 3) - coordinates(:, 3)
            !end if
            do j = 1, size(exclusions, 2)
                exclusion = exclusions(i, j)
                if (exclusion == 0) cycle
                Rji(exclusion, :) = 1e10_rp
            end do
            if (damp_type /= '') then
                allocate (full_damping_factors(N))
                do j = 1, N
                    full_damping_factors(j) = damping_factors(i) * damping_factors(j)
                end do
            end if
            field(i, :) = field(i, :) + multipole_field(multipoles, Rji(:, 1), Rji(:, 2), Rji(:, 3), max_field_order, damp_type, full_damping_factors)
            if (allocated(full_damping_factors)) deallocate (full_damping_factors)
        end do
#ifdef VAR_MPI
        if (mpi_rank == 0) then
            call mpi_reduce(MPI_IN_PLACE, field(1, 1), size(field), MPI_REAL8, MPI_SUM, 0, comm, ierr)
        else
            call mpi_reduce(field(1, 1), field(1, 1), size(field), MPI_REAL8, MPI_SUM, 0, comm, ierr)
        end if
#endif
        deallocate(Rji)
    end subroutine particle_fields_direct

    subroutine multipole_to_local(comm)
#ifdef VAR_MPI
#if defined(USE_MPI_MOD_F90)
        use mpi
#else
        include 'mpif.h'
#endif
#endif
        integer(ip), intent(in) :: comm
        integer(ip) :: ci, N
        real(rp), allocatable :: dx(:), dy(:), dz(:)
        type(node_type) :: node_i
        type(list_type) :: cell_list
        integer(ip) :: work_start, work_stop, work_size
        integer(ip) :: mpi_rank, mpi_size, ierr
        integer(ip) :: periodic_index
        real(rp) :: displacement(3)
        real(rp) :: dummy_damping_factors(0)
        allocate (tree%local_expansion(tree%num_nodes, settings%multipole_size))
        tree%local_expansion = 0.0
#ifdef VAR_MPI
        call mpi_comm_size(comm, mpi_size, ierr)
        call mpi_comm_rank(comm, mpi_rank, ierr)
#else
        mpi_size = 1
        mpi_rank = 0
#endif
        N = tree%num_nodes
        work_size = N / mpi_size
        work_start = 1 + mpi_rank * work_size
        work_stop = (mpi_rank + 1) * work_size
        if (mpi_rank == mpi_size - 1) work_stop = N

        do periodic_index=1, settings%num_interaction_lists
            displacement = matmul(settings%box_matrix, real(settings%aijk(periodic_index, :), rp))
            do ci = work_start, work_stop
                node_i = tree%node_list(ci)%node
                cell_list = node_i%cell_interactions(periodic_index)
                N = size(cell_list%elements)
                allocate (dx(N))
                allocate (dy(N))
                allocate (dz(N))
                dx(:) = node_i%center(1) - (tree%centers(cell_list%elements, 1) + displacement(1))
                dy(:) = node_i%center(2) - (tree%centers(cell_list%elements, 2) + displacement(2))
                dz(:) = node_i%center(3) - (tree%centers(cell_list%elements, 3) + displacement(3))
                tree%local_expansion(node_i%idx, :) = tree%local_expansion(node_i%idx, :) + &
                                                      multipole_field(tree%cell_multipoles(cell_list%elements, :), dx, dy, dz, &
                                                                      settings%expansion_order, '', dummy_damping_factors)
                deallocate (dx)
                deallocate (dy)
                deallocate (dz)
            end do
        end do
#ifdef VAR_MPI
        call mpi_allreduce(MPI_IN_PLACE, tree%local_expansion(1, 1), size(tree%local_expansion), MPI_REAL8, MPI_SUM, comm, ierr)
#endif
    end subroutine multipole_to_local

    subroutine periodic_sum
        integer(ip) :: ai, aj, ak, bi, bj, bk, multipole_size, k, level
        real(rp) :: aijk(3), bijk(3)
        real(rp), allocatable :: cell_multipoles(:,:), cell_multipoles_scaled(:), contribution(:)
        real(rp) :: d_aijk(3), d_bijk(3)
        real(rp), dimension(26*27) :: dx, dy, dz
        real(rp), parameter :: tol_periodic = 1.0e-8_rp
        integer(ip), parameter :: limit = 50
        real(rp) :: dummy_damping_factors(0)
        ! scale up to next 3x3x3 supercell
        ! loop over 26 neighbors to the central cell and get contributions from 27 subcells onto central cell
        multipole_size = size(tree%cell_multipoles, 2)
        ! force charge neutrality
        tree%cell_multipoles(1,1) = 0.0_rp
        allocate(cell_multipoles(26*27, multipole_size))
        allocate(cell_multipoles_scaled(multipole_size))
        allocate(contribution(size(tree%local_expansion, 2)))
        do k=1, 26*27
            cell_multipoles(k, :) = tree%cell_multipoles(tree%root_node%idx, :)
        end do
        cell_multipoles_scaled(:) = cell_multipoles(1, :)
        do level=1, limit
            !if (level == limit) error stop "periodic sum did not converge"
            k = 1
            do ai=-1,1
                do aj=-1,1
                    do ak=-1,1
                        if (ai == 0 .and. aj == 0 .and. ak == 0) cycle
                        ! accumulate for next scaled multipole
                        aijk = [real(ai,rp)*3.0_rp**level, real(aj, rp)*3.0_rp**level, real(ak, rp)*3.0_rp**level]
                        d_aijk = matmul(settings%box_matrix, aijk)
                        call multipole_accumulate_kernel(d_aijk/3.0_rp, &
                                                         cell_multipoles(1,:), &
                                                         cell_multipoles_scaled)
                        ! get positions of the 27 subcells
                        do bi=-1,1
                            do bj=-1,1
                                do bk=-1,1
                                    bijk = [real(bi, rp)*3.0_rp**(level-1), real(bj, rp)*3.0_rp**(level-1), real(bk, rp)*3.0_rp**(level-1)]
                                    d_bijk = matmul(settings%box_matrix, bijk)
                                    dx(k) = d_aijk(1) + d_bijk(1)
                                    dy(k) = d_aijk(2) + d_bijk(2)
                                    dz(k) = d_aijk(3) + d_bijk(3)
                                    k = k + 1
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            ! alternative would be
            ! m2p barnes-hut like step
            ! dx, dy, dz are distances from root cell to image cells
            ! do i=1, size(field, 1)
            !     delta = tree%coordinates(i,:) - tree%root_node%center
            !     contribution = multipole_field(cell_multipoles, dx+delta(1), dy+delta(2), dz+delta(3), max_field_order)
            !     field(i, :) = field(i,:) + contribution
            !     contribution_tracker = contribution_tracker + abs(contribution)
            ! end do
            contribution(:) = multipole_field(cell_multipoles, dx, dy, dz, settings%expansion_order, '', dummy_damping_factors)
            tree%local_expansion(tree%root_node%idx, :) = tree%local_expansion(tree%root_node%idx, :) + contribution

            do k=1, 26*27
                cell_multipoles(k, :) = cell_multipoles_scaled(:)
            end do
            ! print *, 'periodic sum: level, norm2(contribution)', level, norm2(contribution(2:))
            if (norm2(contribution(2:)) < tol_periodic) exit
        end do
    end subroutine periodic_sum

    subroutine dipole_correction(field, max_field_order)
        real(rp), intent(inout) :: field(:, :)
        integer, intent(in) :: max_field_order
        real(rp) :: volume, dipole(3), quadrupole_trace
        real(rp) :: h(3,3)
        integer :: i
        h(:,:) = settings%box_matrix
        volume = h(1,1)*(h(2,2)*h(3,3) - h(2,3)*h(3,2)) &
               - h(1,2)*(h(2,1)*h(3,3) - h(2,3)*h(3,1)) &
               + h(1,3)*(h(1,2)*h(3,2) - h(2,2)*h(3,1))
        dipole = tree%cell_multipoles(tree%root_node%idx, 2:4)
        quadrupole_trace = tree%cell_multipoles(tree%root_node%idx, 5) &
                         + tree%cell_multipoles(tree%root_node%idx, 8) &
                         + tree%cell_multipoles(tree%root_node%idx, 10)
        ! correction to the potential
        
        do i=1, size(field, 1)
            field(i, 1) = field(i, 1) + 4.0_rp * 3.14159265359_rp * sum(dipole*tree%coordinates(i,:)) / (3.0_rp * volume) &
                        + 2.0_rp*3.14159265359_rp*quadrupole_trace/(3.0_rp * volume) !? 
        end do
        if (max_field_order >= 1) then
            do i=1, size(field, 1)
                field(i, 2:4) = field(i, 2:4) + 4.0_rp * 3.14159265359_rp * dipole / (3.0_rp * volume)
            end do
        end if
    end subroutine dipole_correction

    subroutine local_to_local
        integer(ip) :: i, num_nodes_depth, depth
        integer(ip) :: source_order, source_index, sx, sy, sz
        integer(ip) :: target_order, target_index, tx, ty, tz
        real(rp) :: prefactor
        real(rp), allocatable :: dx(:), dy(:), dz(:)
        integer(ip), allocatable :: source_cell_indices(:), target_cell_indices(:)
        type(node_type) :: node
        ! important in case of MPI: can only parallelize from level-1 -> level
        ! might need to make separate lists for this
        do depth = 2, tree%max_depth
            allocate (source_cell_indices(8**depth))
            allocate (target_cell_indices(8**depth))
            num_nodes_depth = 0
            ! illegal to parallelize over depth
            ! can parallelize in one depth-layer
            ! need to do slices of 1:num_nodes_depth
            do i = 1, tree%num_nodes
                node = tree%node_list(i)%node
                if (node%depth == depth) then
                    num_nodes_depth = num_nodes_depth + 1
                    source_cell_indices(num_nodes_depth) = node%parent%idx
                    target_cell_indices(num_nodes_depth) = node%idx
                end if
            end do
            allocate (dx(num_nodes_depth))
            allocate (dy(num_nodes_depth))
            allocate (dz(num_nodes_depth))
            dx(:) = tree%centers(target_cell_indices(1:num_nodes_depth), 1) - tree%centers(source_cell_indices(1:num_nodes_depth), 1)
            dy(:) = tree%centers(target_cell_indices(1:num_nodes_depth), 2) - tree%centers(source_cell_indices(1:num_nodes_depth), 2)
            dz(:) = tree%centers(target_cell_indices(1:num_nodes_depth), 3) - tree%centers(source_cell_indices(1:num_nodes_depth), 3)
            do target_order = 0, settings%expansion_order
                do tx = target_order, 0, -1
                    do ty = target_order, 0, -1
                        do tz = target_order, 0, -1
                            if (tx + ty + tz /= target_order) cycle
                            target_index = xyz2idx(tx, ty, tz)
                            do source_order = target_order, settings%expansion_order
                                do sx = source_order, tx, -1
                                    do sy = source_order, ty, -1
                                        do sz = source_order, tz, -1
                                            if (sx + sy + sz /= source_order) cycle
                                            source_index = xyz2idx(sx, sy, sz)
                                            prefactor = 1.0 / (factorial(sx - tx) * factorial(sy - ty) * factorial(sz - tz))
                                            tree%local_expansion(target_cell_indices(1:num_nodes_depth), target_index) = &
                                                tree%local_expansion(target_cell_indices(1:num_nodes_depth), target_index) &
                                                + prefactor &
                                                * tree%local_expansion(source_cell_indices(1:num_nodes_depth), source_index) &
                                                * (dx**(sx - tx) * dy**(sy - ty) * dz**(sz - tz))
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            deallocate (source_cell_indices)
            deallocate (target_cell_indices)
            deallocate (dx)
            deallocate (dy)
            deallocate (dz)
        end do
    end subroutine local_to_local

    subroutine local_to_particle(comm, field, max_field_order)
#ifdef VAR_MPI
#if defined(USE_MPI_MOD_F90)
        use mpi
#else
        include 'mpif.h'
#endif
#endif
        integer(ip), intent(in) :: comm
        real(rp), intent(inout) :: field(:, :)
        integer(ip), intent(in) :: max_field_order
        type(node_type) :: node
        integer(ip) :: source_order, source_index, sx, sy, sz
        integer(ip) :: target_order, target_index, tx, ty, tz
        real(rp) :: prefactor
        real(rp) :: dx(settings%ncrit), dy(settings%ncrit), dz(settings%ncrit)
        integer(ip) :: work_start, work_stop, work_size, i, num_terminal
        integer(ip) :: mpi_rank, mpi_size, ierr
        integer(ip), allocatable :: terminal_node_list(:)
#ifdef VAR_MPI
        call mpi_comm_size(comm, mpi_size, ierr)
        call mpi_comm_rank(comm, mpi_rank, ierr)
#else
        mpi_size = 1
        mpi_rank = 0
#endif

        if (mpi_rank /= 0) field(:, :) = 0.0

        num_terminal = 0
        do i = 1, tree%num_nodes
            node = tree%node_list(i)%node
            if (.not. node%is_terminal) cycle
            num_terminal = num_terminal + 1
        end do
        allocate (terminal_node_list(num_terminal))
        num_terminal = 0
        do i = 1, tree%num_nodes
            node = tree%node_list(i)%node
            if (.not. node%is_terminal) cycle
            num_terminal = num_terminal + 1
            terminal_node_list(num_terminal) = node%idx
        end do

        work_size = num_terminal / mpi_size
        work_start = 1 + mpi_rank * work_size
        work_stop = (mpi_rank + 1) * work_size
        if (mpi_rank == mpi_size - 1) work_stop = num_terminal

        do i = work_start, work_stop
            node = tree%node_list(terminal_node_list(i))%node
            dx(1:node%nleaf) = tree%coordinates(node%leaf(1:node%nleaf), 1) - node%center(1)
            dy(1:node%nleaf) = tree%coordinates(node%leaf(1:node%nleaf), 2) - node%center(2)
            dz(1:node%nleaf) = tree%coordinates(node%leaf(1:node%nleaf), 3) - node%center(3)
            do target_order = 0, max_field_order
                do tx = target_order, 0, -1
                    do ty = target_order, 0, -1
                        do tz = target_order, 0, -1
                            if (tx + ty + tz /= target_order) cycle
                            target_index = xyz2idx(tx, ty, tz)
                            do source_order = target_order, settings%expansion_order
                                do sx = source_order, tx, -1
                                    do sy = source_order, ty, -1
                                        do sz = source_order, tz, -1
                                            if (sx + sy + sz /= source_order) cycle
                                            source_index = xyz2idx(sx, sy, sz)
                                            prefactor = 1.0 / (factorial(sx - tx) * factorial(sy - ty) * factorial(sz - tz))
                                            field(node%leaf(1:node%nleaf), target_index) = field(node%leaf(1:node%nleaf), target_index) &
                                                                                           + prefactor * tree%local_expansion(node%idx, source_index) * dx**(sx - tx) * dy**(sy - ty) * dz**(sz - tz)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
        deallocate (terminal_node_list)
#ifdef VAR_MPI
        if (mpi_rank == 0) then
            call mpi_reduce(MPI_IN_PLACE, field(1, 1), size(field), MPI_REAL8, MPI_SUM, 0, comm, ierr)
        else
            call mpi_reduce(field(1, 1), field(1, 1), size(field), MPI_REAL8, MPI_SUM, 0, comm, ierr)
        end if
#endif
    end subroutine local_to_particle

    recursive subroutine classify_box(target_coordinates, theta, node, interaction_type)
        real(rp), intent(in) :: target_coordinates(:, :)
        real(rp), intent(in) :: theta
        type(node_type), intent(inout) :: node
        integer(ip), intent(inout) :: interaction_type(:)
        logical(lp) :: too_close
        integer(ip) :: i, octant
        real(rp) :: delta(3), r
        ! 0 -> don't use node (initial value)
        ! 1 -> use node, box
        ! 2 -> use node, particles

        too_close = .false.
        ! any(r*theta < box%r)?
        do i = 1, size(target_coordinates, 1)
            delta = target_coordinates(i, :) - node%center
            r = norm2(delta)
            if (r * theta <= node%rmax) then
                too_close = .true.
            end if
        end do

        if (too_close .and. .not. node%is_terminal) then
            ! node too close and not terminal: visit every child
            do octant = 0, 7
                if (node%occupied(octant)) then
                    call classify_box(target_coordinates, theta, node%children(octant), interaction_type)
                end if
            end do
        else if (too_close .and. node%is_terminal) then
            ! node too close and terminal, use particles
            interaction_type(node%idx) = 2
        else
            ! node is far away, OK to use box multipole
            interaction_type(node%idx) = 1
        end if
    end subroutine classify_box

    subroutine get_multipoles(comm, coordinates, multipoles, target_coordinates, theta, ncrit, expansion_order, idx_near_field, box_coordinates, box_multipoles)
        integer(ip), intent(in) :: comm
        real(rp), intent(in), target :: coordinates(:, :)
        real(rp), intent(in), target :: multipoles(:, :, :)
        real(rp), intent(in) :: target_coordinates(:, :)
        real(rp), intent(in) :: theta
        integer(ip), intent(in) :: ncrit
        integer(ip), intent(in) :: expansion_order
        real(rp), intent(out), allocatable :: box_coordinates(:, :)
        real(rp), intent(out), allocatable :: box_multipoles(:, :, :)
        integer(ip), intent(out), allocatable :: idx_near_field(:)
        real(rp), pointer :: p_coordinates(:, :)
        real(rp), pointer :: p_multipoles(:, :)
        integer(ip) :: i, j, k, d, ndens, leaf, num_near_field, num_far_field
        integer(ip), allocatable :: interaction_type(:)
        p_coordinates => coordinates(:,:)
        p_multipoles => multipoles(:,:,1)
        ndens = size(multipoles, 3)
        call octree_build(ncrit, expansion_order, theta, p_coordinates, p_multipoles)
        allocate (interaction_type(tree%num_nodes), source=0)
        call classify_box(target_coordinates, theta, tree%root_node, interaction_type)
        ! find dimension of output coordinates/multipoles
        ! 0 -> don't use node
        ! 1 -> use node, box
        ! 2 -> use node, particles
        num_near_field = 0
        num_far_field = 0
        do i = 1, tree%num_nodes
            if (interaction_type(i) == 0) then
                continue
            else if (interaction_type(i) == 1) then
                num_far_field = num_far_field + 1
            else if (interaction_type(i) == 2) then
                num_near_field = num_near_field + tree%node_list(i)%node%nleaf
            else
                error stop "Wrong interaction_type"
            end if
        end do

        if (allocated(box_coordinates)) deallocate (box_coordinates)
        if (allocated(box_multipoles)) deallocate (box_multipoles)
        if (allocated(idx_near_field)) deallocate (idx_near_field)
        allocate(box_coordinates(num_far_field, 3))
        allocate(box_multipoles(num_far_field, (expansion_order + 1) * (expansion_order + 2) * (expansion_order + 3) / 6, ndens))
        allocate(idx_near_field(num_near_field))
        ! output: near-field
        k = 1
        do i = 1, tree%num_nodes
            if (interaction_type(i) == 2) then
                do j = 1, tree%node_list(i)%node%nleaf
                    leaf = tree%node_list(i)%node%leaf(j)
                    idx_near_field(k) = leaf
                    k = k + 1
                end do
            end if
        end do
        ! output: far-field - place resulting multipoles in output arrays
        do d = 1, ndens
            call finalize
            p_coordinates => coordinates(:,:)
            p_multipoles => multipoles(:,:,d)
            call octree_build(ncrit, expansion_order, theta, p_coordinates, p_multipoles)
            call multipole_expansion(comm)
            k = 1
            do i = 1, tree%num_nodes
                if (interaction_type(i) == 1) then
                    box_coordinates(k, :) = tree%centers(i, :)
                    box_multipoles(k, :, d) = tree%cell_multipoles(i, :)
                    k = k + 1
                end if 
            end do
        end do
        call finalize
    end subroutine get_multipoles

    subroutine get_boxes(coordinates, target_coordinates, theta, ncrit, expansion_order, idx_near_field, box_coordinates, box_indices)
        real(rp), intent(in), target :: coordinates(:, :)
        real(rp), intent(in) :: target_coordinates(:, :)
        real(rp), intent(in) :: theta
        integer(ip), intent(in) :: ncrit
        integer(ip), intent(in) :: expansion_order
        real(rp), intent(out), allocatable :: box_coordinates(:, :)
        integer(ip), intent(out), allocatable :: idx_near_field(:), box_indices(:)
        real(rp), pointer :: p_coordinates(:, :)
        real(rp), pointer :: p_multipoles(:, :) ! dummy
        integer(ip) :: i, j, k, leaf, num_near_field, num_far_field
        integer(ip), allocatable :: interaction_type(:)
        p_coordinates => coordinates(:,:)
        call octree_build(ncrit, expansion_order, theta, p_coordinates, p_multipoles)
        allocate (interaction_type(tree%num_nodes), source=0)
        call classify_box(target_coordinates, theta, tree%root_node, interaction_type)
        ! find dimension of output coordinates/multipoles
        ! 0 -> don't use node
        ! 1 -> use node, box
        ! 2 -> use node, particles
        num_near_field = 0
        num_far_field = 0
        do i = 1, tree%num_nodes
            if (interaction_type(i) == 0) then
                continue
            else if (interaction_type(i) == 1) then
                num_far_field = num_far_field + 1
            else if (interaction_type(i) == 2) then
                num_near_field = num_near_field + tree%node_list(i)%node%nleaf
            else
                error stop "Wrong interaction_type"
            end if
        end do
        if (allocated(box_coordinates)) deallocate (box_coordinates)
        if (allocated(box_indices)) deallocate (box_indices)
        if (allocated(idx_near_field)) deallocate (idx_near_field)
        allocate(box_coordinates(num_far_field, 3))
        allocate(box_indices(num_far_field))
        allocate(idx_near_field(num_near_field))
        ! output: near-field
        k = 1
        do i = 1, tree%num_nodes
            if (interaction_type(i) == 2) then
                do j = 1, tree%node_list(i)%node%nleaf
                    leaf = tree%node_list(i)%node%leaf(j)
                    idx_near_field(k) = leaf
                    k = k + 1
                end do
            end if
        end do
        ! output: far-field
        k = 1
        do i = 1, tree%num_nodes
            if (interaction_type(i) == 1) then
                box_coordinates(k, :) = tree%centers(i, :)
                box_indices(k) = i
                k = k + 1
            end if 
        end do
        call finalize
    end subroutine get_boxes

    subroutine interpolate_electron_fields(comm, coordinates, local_expansion, box_indices, theta, ncrit, expansion_order, field_order, field)
        integer(ip), intent(in) :: comm
        real(rp), intent(in), target :: coordinates(:, :)
        real(rp), intent(in) :: local_expansion(:, :, :)
        integer(ip), intent(in) :: box_indices(:)
        real(rp), intent(in) :: theta
        integer(ip), intent(in) :: ncrit
        integer(ip), intent(in) :: expansion_order
        integer(ip), intent(in) :: field_order
        real(rp), intent(out) :: field(size(coordinates, 1), (field_order + 1) * (field_order + 2) * (field_order + 3) / 6, size(local_expansion, 3))
        real(rp), pointer :: p_coordinates(:, :)
        real(rp), pointer :: p_multipoles(:, :)
        integer(ip) :: i, k, box_index, ndens
        p_coordinates => coordinates(:,:)
        call octree_build(ncrit, expansion_order, theta, p_coordinates, p_multipoles)
        ndens = size(local_expansion, 3)
        allocate (tree%local_expansion(tree%num_nodes, settings%multipole_size), source=0.0_rp)
        do k = 1, ndens
            tree%local_expansion(:,:) = 0.0_rp
            do i = 1, size(box_indices)
                box_index = box_indices(i)
                tree%local_expansion(box_index, :) = local_expansion(i, :, k)
            end do
            call local_to_local
            call local_to_particle(comm, field(:,:,k), field_order)
        end do
        call finalize
    end subroutine interpolate_electron_fields

    subroutine field_direct(comm, coordinates, multipoles, exclusions, field_order, field, damp_type, damping_factors)
        integer(ip), intent(in) :: comm
        real(rp), intent(in), target :: coordinates(:, :)
        real(rp), intent(in), target :: multipoles(:, :)
        integer(ip), intent(in) :: exclusions(:, :)
        integer(ip), intent(in) :: field_order
        character(len=*), intent(in) :: damp_type
        real(rp), intent(in) :: damping_factors(:)
        real(rp), intent(out) :: field(size(coordinates, 1), (field_order + 1) * (field_order + 2) * (field_order + 3) / 6)
        real(rp), pointer :: p_coordinates(:, :)
        real(rp), pointer :: p_multipoles(:, :)
        p_coordinates => coordinates
        p_multipoles => multipoles
        call particle_fields_direct(comm, field, coordinates, multipoles, exclusions, field_order, damp_type, damping_factors)
    end subroutine field_direct

    subroutine field_fmm(comm, coordinates, multipoles, exclusions, theta, ncrit, expansion_order, field_order, field,&
                         damp_type, damping_factors, is_periodic, box_matrix)
        integer(ip), intent(in) :: comm
        real(rp), intent(in), target :: coordinates(:, :)
        real(rp), intent(in), target :: multipoles(:, :)
        integer(ip), intent(in) :: exclusions(:, :)
        real(rp), intent(in) :: theta
        integer(ip), intent(in) :: ncrit
        integer(ip), intent(in) :: expansion_order
        integer(ip), intent(in) :: field_order
        character(len=*), intent(in) :: damp_type ! 'THOLE', 'AMOEBA', 'ERF' or '' (off)
        real(rp), intent(in) :: damping_factors(:)
        logical(lp), intent(in) :: is_periodic
        real(rp), intent(in), dimension(3,3) :: box_matrix
        real(rp), intent(out) :: field(size(coordinates, 1), (field_order + 1) * (field_order + 2) * (field_order + 3) / 6)
        real(rp), pointer :: p_coordinates(:, :)
        real(rp), pointer :: p_multipoles(:, :)

        if (.not. (trim(damp_type) == 'THOLE' .or. &
                   trim(damp_type) == 'AMOEBA' .or. &
                   trim(damp_type) == 'ERF' .or. &
                   trim(damp_type) == '')) then
            error stop "Unknown damp type: " // damp_type
        end if
        settings%is_periodic = is_periodic
        settings%box_matrix = box_matrix

        if (size(coordinates, 1) <= ncrit) then
            if (settings%is_periodic) error stop "Not enough coordinates for periodic FMM, please decrease ncrit"
            call field_direct(comm, coordinates, multipoles, exclusions, field_order, field, damp_type, damping_factors)
        else
            p_coordinates => coordinates
            p_multipoles => multipoles
            call octree_build(ncrit, expansion_order, theta, p_coordinates, p_multipoles)
            call multipole_expansion(comm)
            call build_interaction_lists_fmm
            call particle_fields(comm, field, exclusions, field_order, damp_type, damping_factors)
            call multipole_to_local(comm)
            if (settings%is_periodic) then
                call periodic_sum
                call dipole_correction(field, field_order)
            end if
            call local_to_local
            call local_to_particle(comm, field, field_order)
            call finalize
        end if
    end subroutine field_fmm
end module fmm
