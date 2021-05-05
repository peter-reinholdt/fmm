module octree
    use t_tensor
    use t_tensor_damp_thole
    use t_tensor_damp_amoeba
    use t_tensor_damp_erf
    implicit none
    private

    public field_direct
    public field_fmm

    type settings_type
        integer :: ncrit
        real(8) :: theta
        integer :: expansion_order
        integer :: multipole_size
    end type settings_type

    type list_type
        integer :: head = 1
        integer, allocatable :: elements(:)
    end type list_type

    type node_type
        ! tree things
        integer :: idx
        integer :: depth
        integer :: nleaf
        integer, allocatable :: leaf(:)
        real(8) :: center(3)
        real(8) :: r
        real(8) :: rmax
        type(node_type), pointer :: parent
        type(node_type), pointer :: children(:)
        logical :: occupied(0:7) = .FALSE.
        logical :: is_terminal = .TRUE.
        ! data
        type(list_type) :: particle_interactions
        type(list_type) :: cell_interactions
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
        real(8), pointer :: coordinates(:, :)
        real(8), pointer :: source_multipoles(:, :)
        type(node_type), pointer :: root_node
        integer :: num_nodes
        integer :: max_depth
        type(node_list_type), allocatable :: node_list(:)
        real(8), allocatable :: centers(:, :)
        real(8), allocatable :: cell_multipoles(:, :)
        real(8), allocatable :: local_expansion(:, :)
    end type tree_type

    type(settings_type) settings
    type(tree_type) tree

contains

    pure function xyz2idx(x, y, z) result(idx)
        implicit none
        integer, intent(in) :: x, y, z
        integer :: k
        integer :: idx
        k = x + y + z
        ! number of components before order k is k*(k+1)*(k+2)/6
        ! (y**2 + 2*y*z + y + z**2 + 3*z)/2 + 1 is the symmetry-packed index of the current slice
        idx = k * (k + 1) * (k + 2) / 6 + (y**2 + 2 * y * z + y + z**2 + 3 * z) / 2 + 1
    end function xyz2idx

    pure function factorial(n)
        integer, intent(in) :: n
        integer :: i
        real(8) :: factorial
        factorial = 1
        do i = n, 1, -1
            factorial = factorial * i
        end do
    end function

    pure function trinom(i, j, k)
        integer, intent(in) :: i, j, k
        real(8) :: trinom
        trinom = factorial(i + j + k) / (factorial(i) * factorial(j) * factorial(k))
    end function trinom

    pure function binom(n, k)
        integer, intent(in) :: n, k
        real(8) :: binom
        binom = factorial(n) / (factorial(k) * factorial(n - k))
    end function binom

    subroutine list_append(list, element)
        type(list_type), intent(inout) :: list
        integer, intent(in) :: element
        integer, allocatable :: swap(:)
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
        integer, allocatable :: swap(:)
        call move_alloc(list%elements, swap)
        allocate (list%elements(1:list%head - 1))
        list%elements(1:list%head - 1) = swap(1:list%head - 1)
        deallocate (swap)
    end subroutine list_trim

    subroutine octree_build(ncrit, expansion_order, theta, coordinates, multipoles)
        integer, intent(in) :: ncrit
        integer, intent(in) :: expansion_order
        real(8), intent(in) :: theta
        real(8), intent(in), pointer :: coordinates(:, :)
        real(8), intent(in), pointer :: multipoles(:, :)
        integer :: i
        integer :: octant
        integer :: multipole_size
        integer :: n
        type(node_type), pointer :: node

        settings%ncrit = ncrit
        settings%expansion_order = expansion_order
        settings%theta = theta
        n = settings%expansion_order
        multipole_size = (n + 1) * (n + 2) * (n + 3) / 6
        settings%multipole_size = multipole_size

        if (.not. associated(tree%root_node)) allocate (tree%root_node)
        tree%root_node%nleaf = 0
        allocate (tree%root_node%leaf(settings%ncrit))
        allocate (tree%root_node%children(0:7))
        nullify (tree%root_node%parent)
        tree%num_nodes = 1
        tree%root_node%depth = 1
        tree%root_node%center = sum(coordinates, 1) / size(coordinates, 1)
        tree%root_node%r = max( &
                           maxval(abs(coordinates(:, 1) - tree%root_node%center(1))), &
                           maxval(abs(coordinates(:, 2) - tree%root_node%center(2))), &
                           maxval(abs(coordinates(:, 3) - tree%root_node%center(3)))) * 1.0001
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
        integer, intent(in) :: octant
        real(8) :: r
        real(8) :: center(3)

        tree%num_nodes = tree%num_nodes + 1
        node%occupied(octant) = .TRUE.
        node%is_terminal = .FALSE.

        r = node%r * 0.5
        center = node%center
        center(1) = center(1) + r * (IAND(octant, 1) * 2 - 1)
        center(2) = center(2) + r * (IAND(octant, 2) - 1)
        center(3) = center(3) + r * (IAND(octant, 4) / 2 - 1)

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
        integer :: i
        integer :: leaf
        integer :: octant
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
        integer :: pos, i
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
        integer, intent(inout) :: pos
        integer :: octant
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
    end subroutine finalize

    recursive subroutine clean_node(node)
        type(node_type), intent(inout) :: node
        integer octant
        do octant = 0, 7
            if (node%occupied(octant)) then
                call clean_node(node%children(octant))
            end if
        end do
        deallocate (node%children)
    end subroutine clean_node

    subroutine build_interaction_lists_fmm
        integer :: i
        do i = 1, tree%num_nodes
            ! just temporary allocation, will be resized to fit later
            allocate (tree%node_list(i)%node%particle_interactions%elements(64))
            allocate (tree%node_list(i)%node%cell_interactions%elements(64))
        end do
        call interact_fmm(tree%root_node, tree%root_node, settings%theta)
        do i = 1, tree%num_nodes
            call list_trim(tree%node_list(i)%node%particle_interactions)
            call list_trim(tree%node_list(i)%node%cell_interactions)
        end do
    end subroutine build_interaction_lists_fmm

    recursive subroutine interact_fmm(node_i, node_j, theta)
        type(node_type) :: node_i, node_j
        integer :: octant
        real(8) :: theta
        real(8) :: delta(3)
        real(8) :: r
        delta = node_i%center - node_j%center
        r = norm2(delta)
        if (r * theta > node_i%rmax + node_j%rmax) then
            ! far field
            call list_append(node_i%cell_interactions, node_j%idx)
        else if (node_i%is_terminal .and. node_j%is_terminal) then
            ! both are leaf nodes, and should interact directly
            call list_append(node_i%particle_interactions, node_j%idx)
        else if (node_j%is_terminal .or. ((node_i%rmax > node_j%rmax) .and. (.not. node_i%is_terminal))) then
            ! too big, descend i
            do octant = 0, 7
                if (node_i%occupied(octant)) then
                    call interact_fmm(node_i%children(octant), node_j, theta)
                end if
            end do
        else
            ! too big, descend j
            do octant = 0, 7
                if (node_j%occupied(octant)) then
                    call interact_fmm(node_i, node_j%children(octant), theta)
                end if
            end do
        end if
    end subroutine interact_fmm

    recursive subroutine multipole_accumulate(source_node_index, target_node_index)
        integer, intent(in) :: source_node_index, target_node_index
        real(8) :: delta(3)
        integer :: source_index, source_order, sx, sy, sz
        integer :: target_index, target_order, tx, ty, tz
        real(8) :: symfac
        integer :: max_multipole_order
        real(8) :: monomial

        if (tree%node_list(target_node_index)%node%idx /= 1) then
            call multipole_accumulate(source_node_index, tree%node_list(target_node_index)%node%parent%idx)
        end if
        delta = tree%centers(source_node_index, :) - tree%centers(target_node_index, :)
        max_multipole_order = NINT((size(tree%cell_multipoles, 2) * 6)**(1./3.)) - 2

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
                                        tree%cell_multipoles(target_node_index, target_index) = tree%cell_multipoles(target_node_index, target_index) &
                                                                                                + symfac * monomial * tree%cell_multipoles(source_node_index, source_index)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine multipole_accumulate

    subroutine shift_multipole_vec(source_multipole, target_multipole, dx, dy, dz)
        real(8), intent(in) :: source_multipole(:, :)
        real(8), intent(inout) :: target_multipole(:)
        real(8), intent(in) :: dx(:), dy(:), dz(:)
        integer :: source_index, source_order, sx, sy, sz
        integer :: target_index, target_order, tx, ty, tz
        real(8) :: symfac
        integer :: max_multipole_order
        real(8) :: monomial(size(dx))

        max_multipole_order = NINT((size(source_multipole, 2) * 6)**(1./3.)) - 2
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
        integer, intent(in) :: comm
        integer :: i, num_terminal
        real(8), allocatable :: dx(:), dy(:), dz(:)
        type(node_type) :: node
        integer :: mpi_size, mpi_rank, ierr, work_size, work_start, work_stop, node_idx
        integer, allocatable :: terminal_node_list(:)

        allocate (tree%cell_multipoles(tree%num_nodes, settings%multipole_size), source=0.0D0)
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
            dx = tree%coordinates(node%leaf(1:node%nleaf), 1) - node%center(1)
            dy = tree%coordinates(node%leaf(1:node%nleaf), 2) - node%center(2)
            dz = tree%coordinates(node%leaf(1:node%nleaf), 3) - node%center(3)
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
        real(8), allocatable :: F(:)
        real(8), intent(in) :: multipole(:, :)
        real(8), intent(in) :: x(:), y(:), z(:)
        integer, intent(in) :: max_field_order
        character(len=6), intent(in), optional :: damp_type
        real(8), intent(in), optional :: damping_factors(:)
        integer :: max_multipole_order, multipole_order, sx, sy, sz
        integer :: field_order, tx, ty, tz
        integer :: multipole_index, field_index, tensor_index
        integer :: max_tensor_order
        real(8) :: symfac
        real(8) :: taylor
        integer :: length, N
        real(8), allocatable :: T(:, :)

        N = size(x)

        length = (max_field_order + 1) * (max_field_order + 2) * (max_field_order + 3) / 6
        allocate (F(length), source=0.0D0)
        max_multipole_order = NINT((size(multipole, 2) * 6)**(1./3.)) - 2

        max_tensor_order = max_field_order + max_multipole_order

        length = (max_tensor_order + 1) * (max_tensor_order + 2) * (max_tensor_order + 3) / 6
        allocate (T(N, length), source=0.0D0)

        ! T(i, 1) -> T0
        ! T(i, 2) -> T1_x
        ! T(i, 3) -> T1_y
        ! T(i, 4) -> T1_z
        ! T(i, 5) -> T2_xx
        ! T(i, 6) -> T2_xy
        ! ...
        if (present(damp_type)) then
            if (damp_type == 'AMOEBA') then
                call Tn_damp_amoeba(max_tensor_order, x, y, z, damping_factors, T)
            else if (damp_type == 'THOLE ') then
                call Tn_damp_thole(max_tensor_order, x, y, z, damping_factors, T)
            else if (damp_type == 'ERF   ') then
                call Tn_damp_erf(max_tensor_order, x, y, z, damping_factors, T)
            end if
        else
            call Tn(max_tensor_order, x, y, z, T)
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
        integer, intent(in) :: comm
        real(8), intent(inout) :: field(:, :)
        integer, intent(in) :: exclusions(:, :)
        integer, intent(in) :: max_field_order
        character(len=6), intent(in), optional :: damp_type
        real(8), intent(in), optional :: damping_factors(:)
        real(8), allocatable :: full_damping_factors(:)

        integer :: j, ci, cj, leaf_i, leaf_j
        real(8)  :: dx(settings%ncrit), dy(settings%ncrit), dz(settings%ncrit)
        type(list_type) :: particle_list
        type(node_type) :: node_i, node_j
        integer :: mpi_size, mpi_rank, ierr
        integer, allocatable :: work_node_list(:), work_leaf_list(:)
        integer :: work_start, work_stop, work_idx, N, total_work
        integer :: ideal_work_size, ideal_work_start, ideal_work_stop
        logical :: set_work_start, set_work_stop
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
        total_work = 0
        do ci = 1, tree%num_nodes
            node_i = tree%node_list(ci)%node
            if (.not. node_i%is_terminal) cycle
            total_work = total_work + node_i%nleaf * node_i%particle_interactions%head
        end do

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
                total_work = total_work + node_i%particle_interactions%head
            end do
        end do
        if (mpi_rank == mpi_size - 1) work_stop = N

        do ci = work_start, work_stop
            node_i = tree%node_list(work_node_list(ci))%node
            particle_list = node_i%particle_interactions
            leaf_i = work_leaf_list(ci)
            do cj = 1, size(particle_list%elements)
                node_j = tree%node_list(particle_list%elements(cj))%node
                dx(1:node_j%nleaf) = tree%coordinates(leaf_i, 1) - tree%coordinates(node_j%leaf(1:node_j%nleaf), 1)
                dy(1:node_j%nleaf) = tree%coordinates(leaf_i, 2) - tree%coordinates(node_j%leaf(1:node_j%nleaf), 2)
                dz(1:node_j%nleaf) = tree%coordinates(leaf_i, 3) - tree%coordinates(node_j%leaf(1:node_j%nleaf), 3)
                do j = 1, node_j%nleaf
                    leaf_j = node_j%leaf(j)
                    if (any(exclusions(leaf_i, :) == leaf_j)) then
                        dx(j) = 1d10 ! T(x,y,z) -> 0. as x,y,z -> large numbers
                        dy(j) = 1d10 ! a bit of a hacky way to do exclusions
                        dz(j) = 1d10 ! but makes it easy to pass the entire vector into multipole_field
                    end if
                end do
                if (present(damp_type)) then
                    allocate (full_damping_factors(node_j%nleaf))
                    do j = 1, node_j%nleaf
                        leaf_j = node_j%leaf(j)
                        full_damping_factors(j) = damping_factors(leaf_i) * damping_factors(leaf_j)
                    end do
                    field(leaf_i, :) = field(leaf_i, :) + multipole_field(tree%source_multipoles(node_j%leaf(1:node_j%nleaf), :), &
                                                                          dx(1:node_j%nleaf), dy(1:node_j%nleaf), dz(1:node_j%nleaf), max_field_order, &
                                                                          damp_type, full_damping_factors)
                    deallocate (full_damping_factors)
                else
                    field(leaf_i, :) = field(leaf_i, :) + multipole_field(tree%source_multipoles(node_j%leaf(1:node_j%nleaf), :), &
                                                                          dx(1:node_j%nleaf), dy(1:node_j%nleaf), dz(1:node_j%nleaf), max_field_order)
                end if
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
#ifdef VAR_MPI
#if defined(USE_MPI_MOD_F90)
        use mpi
#else
        include 'mpif.h'
#endif
#endif
        integer, intent(in) :: comm
        real(8), intent(inout) :: field(:, :)
        real(8), intent(in) :: coordinates(:, :)
        real(8), intent(in) :: multipoles(:, :)
        integer, intent(in) :: exclusions(:, :)
        integer, intent(in) :: max_field_order
        character(len=6), intent(in), optional :: damp_type
        real(8), intent(in), optional :: damping_factors(:)
        real(8), allocatable :: full_damping_factors(:)
        integer :: i, j, exclusion
        integer :: N
        integer :: work_start, work_stop, work_size
        integer :: mpi_rank, mpi_size, ierr
        real(8) :: dx(size(coordinates, 1)), dy(size(coordinates, 1)), dz(size(coordinates, 1))
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

        do i = work_start, work_stop
            dx = coordinates(i, 1) - coordinates(:, 1)
            dy = coordinates(i, 2) - coordinates(:, 2)
            dz = coordinates(i, 3) - coordinates(:, 3)
            do j = 1, size(exclusions, 2)
                exclusion = exclusions(i, j)
                if (exclusion == 0) cycle
                dx(exclusion) = 1d10
                dy(exclusion) = 1d10
                dz(exclusion) = 1d10
            end do
            if (present(damp_type)) then
                allocate (full_damping_factors(N))
                do j = 1, N
                    full_damping_factors(j) = damping_factors(i) * damping_factors(j)
                end do
                field(i, :) = field(i, :) + multipole_field(multipoles, dx, dy, dz, max_field_order, damp_type, full_damping_factors)
                deallocate (full_damping_factors)
            else
                field(i, :) = field(i, :) + multipole_field(multipoles, dx, dy, dz, max_field_order)
            end if
        end do
#ifdef VAR_MPI
        if (mpi_rank == 0) then
            call mpi_reduce(MPI_IN_PLACE, field(1, 1), size(field), MPI_REAL8, MPI_SUM, 0, comm, ierr)
        else
            call mpi_reduce(field(1, 1), field(1, 1), size(field), MPI_REAL8, MPI_SUM, 0, comm, ierr)
        end if
#endif
    end subroutine particle_fields_direct

    subroutine multipole_to_local(comm)
#ifdef VAR_MPI
#if defined(USE_MPI_MOD_F90)
        use mpi
#else
        include 'mpif.h'
#endif
#endif
        integer, intent(in) :: comm
        integer :: ci, N
        real(8), allocatable :: dx(:), dy(:), dz(:)
        type(node_type) :: node_i
        type(list_type) :: cell_list
        integer :: work_start, work_stop, work_size
        integer :: mpi_rank, mpi_size, ierr
        allocate (tree%local_expansion(tree%num_nodes, settings%multipole_size), source=0.0D0)
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

        do ci = work_start, work_stop
            node_i = tree%node_list(ci)%node
            cell_list = node_i%cell_interactions
            N = size(cell_list%elements)
            allocate (dx(N))
            allocate (dy(N))
            allocate (dz(N))
            dx = node_i%center(1) - tree%centers(cell_list%elements, 1)
            dy = node_i%center(2) - tree%centers(cell_list%elements, 2)
            dz = node_i%center(3) - tree%centers(cell_list%elements, 3)
            tree%local_expansion(node_i%idx, :) = tree%local_expansion(node_i%idx, :) + &
                                                  multipole_field(tree%cell_multipoles(cell_list%elements, :), dx, dy, dz, settings%expansion_order)
            deallocate (dx)
            deallocate (dy)
            deallocate (dz)
        end do
#ifdef VAR_MPI
        call mpi_allreduce(MPI_IN_PLACE, tree%local_expansion(1, 1), size(tree%local_expansion), MPI_REAL8, MPI_SUM, comm, ierr)
#endif
    end subroutine multipole_to_local

    subroutine local_to_local
        integer :: i, num_nodes_depth, depth
        integer :: source_order, source_index, sx, sy, sz
        integer :: target_order, target_index, tx, ty, tz
        real(8) :: prefactor
        real(8), allocatable :: dx(:), dy(:), dz(:)
        integer, allocatable :: source_cell_indices(:), target_cell_indices(:)
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
            dx = tree%centers(target_cell_indices(1:num_nodes_depth), 1) - tree%centers(source_cell_indices(1:num_nodes_depth), 1)
            dy = tree%centers(target_cell_indices(1:num_nodes_depth), 2) - tree%centers(source_cell_indices(1:num_nodes_depth), 2)
            dz = tree%centers(target_cell_indices(1:num_nodes_depth), 3) - tree%centers(source_cell_indices(1:num_nodes_depth), 3)
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
                                            prefactor = 1.0 / real(factorial(sx - tx) * factorial(sy - ty) * factorial(sz - tz))
                                            tree%local_expansion(target_cell_indices(1:num_nodes_depth), target_index) = tree%local_expansion(target_cell_indices(1:num_nodes_depth), target_index) &
                                                                                                                         + prefactor * tree%local_expansion(source_cell_indices(1:num_nodes_depth), source_index) * (dx**(sx - tx) * dy**(sy - ty) * dz**(sz - tz))
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
        integer, intent(in) :: comm
        real(8), intent(inout) :: field(:, :)
        integer, intent(in) :: max_field_order
        type(node_type) :: node
        integer :: source_order, source_index, sx, sy, sz
        integer :: target_order, target_index, tx, ty, tz
        real(8) :: prefactor
        real(8) :: dx(settings%ncrit), dy(settings%ncrit), dz(settings%ncrit)
        integer :: work_start, work_stop, work_size, i, num_terminal
        integer :: mpi_rank, mpi_size, ierr
        integer, allocatable :: terminal_node_list(:)
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
                                            prefactor = 1.0 / real(factorial(sx - tx) * factorial(sy - ty) * factorial(sz - tz))
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
        real(8), intent(in) :: target_coordinates(:, :)
        real(8), intent(in) :: theta
        type(node_type), intent(inout) :: node
        integer, intent(inout) :: interaction_type(:)
        logical :: too_close
        integer :: i, octant
        real(8) :: delta(3), r
        ! 0 -> don't use node (initial value)
        ! 1 -> use node, box
        ! 2 -> use node, particles

        too_close = .false.
        ! any(r*theta < box%r)?
        do i = 1, size(target_coordinates, 1)
            delta = target_coordinates(i, :) - node%center
            delta(3) = 0.0d0
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

    subroutine get_multipoles(comm, coordinates, multipoles, target_coordinates, theta, ncrit, expansion_order, result_coordinates, result_multipoles, result_sizes)
        integer, intent(in) :: comm
        real(8), intent(in), target :: coordinates(:, :)
        real(8), intent(in), target :: multipoles(:, :)
        real(8), intent(in) :: target_coordinates(:, :)
        real(8), intent(in) :: theta
        integer, intent(in) :: ncrit
        integer, intent(in) :: expansion_order
        real(8), intent(out), allocatable :: result_coordinates(:, :)
        real(8), intent(out), allocatable :: result_multipoles(:, :)
        real(8), intent(out), allocatable :: result_sizes(:)
        real(8), pointer :: p_coordinates(:, :)
        real(8), pointer :: p_multipoles(:, :)
        real(8), pointer :: p_target_coordinates(:, :)
        integer :: i, j, k, leaf, output_size
        integer, allocatable :: interaction_type(:)
        p_coordinates => coordinates
        p_multipoles => multipoles
        call octree_build(ncrit, expansion_order, theta, p_coordinates, p_multipoles)
        call multipole_expansion(comm)
        allocate (interaction_type(tree%num_nodes), source=0)
        call classify_box(target_coordinates, theta, tree%root_node, interaction_type)
        ! find dimension of output coordinates/multipoles
        ! 0 -> don't use node
        ! 1 -> use node, box
        ! 2 -> use node, particles
        output_size = 0
        do i = 1, tree%num_nodes
            if (interaction_type(i) == 0) then
                continue
            else if (interaction_type(i) == 1) then
                output_size = output_size + 1
            else if (interaction_type(i) == 2) then
                output_size = output_size + tree%node_list(i)%node%nleaf
            else
                error stop "Wrong interaction_type"
            end if
        end do

        if (allocated(result_coordinates)) deallocate (result_coordinates)
        if (allocated(result_multipoles)) deallocate (result_multipoles)
        if (allocated(result_sizes)) deallocate (result_sizes)

        allocate (result_coordinates(output_size, 3))
        allocate (result_multipoles(output_size, (expansion_order + 1) * (expansion_order + 2) * (expansion_order + 3) / 6))
        allocate (result_sizes(output_size))
        ! place resulting multipoles in output arrays
        k = 1
        do i = 1, tree%num_nodes
            if (interaction_type(i) == 0) then
                continue
            else if (interaction_type(i) == 1) then
                result_coordinates(k, :) = tree%centers(i, :)
                result_multipoles(k, :) = tree%cell_multipoles(i, :)
                result_sizes(k) = tree%node_list(i)%node%rmax 
                k = k + 1
            else if (interaction_type(i) == 2) then
                do j = 1, tree%node_list(i)%node%nleaf
                    leaf = tree%node_list(i)%node%leaf(j)
                    result_coordinates(k, :) = coordinates(leaf, :)
                    result_multipoles(k, :) = multipoles(leaf, :)
                    result_sizes(k) = 0.0d0
                    k = k + 1
                end do
            else
                error stop 'Wrong interaction_type'
            end if
        end do
    end subroutine get_multipoles

    subroutine field_direct(comm, coordinates, multipoles, exclusions, field_order, field, damp_type, damping_factors)
        integer, intent(in) :: comm
        real(8), intent(in), target :: coordinates(:, :)
        real(8), intent(in), target :: multipoles(:, :)
        integer, intent(in) :: exclusions(:, :)
        integer, intent(in) :: field_order
        character(len=6), intent(in), optional :: damp_type
        real(8), intent(in), optional :: damping_factors(:)
        real(8), intent(out) :: field(size(coordinates, 1), (field_order + 1) * (field_order + 2) * (field_order + 3) / 6)
        real(8), pointer :: p_coordinates(:, :)
        real(8), pointer :: p_multipoles(:, :)
        p_coordinates => coordinates
        p_multipoles => multipoles
        if (present(damp_type)) then
            call particle_fields_direct(comm, field, coordinates, multipoles, exclusions, field_order, damp_type, damping_factors)
        else
            call particle_fields_direct(comm, field, coordinates, multipoles, exclusions, field_order)
        end if
    end subroutine field_direct

    subroutine field_fmm(comm, coordinates, multipoles, exclusions, theta, ncrit, expansion_order, field_order, field, damp_type, damping_factors)
        integer, intent(in) :: comm
        real(8), intent(in), target :: coordinates(:, :)
        real(8), intent(in), target :: multipoles(:, :)
        integer, intent(in) :: exclusions(:, :)
        real(8), intent(in) :: theta
        integer, intent(in) :: ncrit
        integer, intent(in) :: expansion_order
        integer, intent(in) :: field_order
        character(len=6), intent(in), optional :: damp_type
        real(8), intent(in), optional :: damping_factors(:)
        real(8), intent(out) :: field(size(coordinates, 1), (field_order + 1) * (field_order + 2) * (field_order + 3) / 6)
        real(8), pointer :: p_coordinates(:, :)
        real(8), pointer :: p_multipoles(:, :)
        if (size(coordinates, 1) <= ncrit) then
            if (present(damp_type)) then
                call field_direct(comm, coordinates, multipoles, exclusions, field_order, field, damp_type, damping_factors)
            else
                call field_direct(comm, coordinates, multipoles, exclusions, field_order, field)
            end if
        else
            p_coordinates => coordinates
            p_multipoles => multipoles
            call octree_build(ncrit, expansion_order, theta, p_coordinates, p_multipoles)
            call multipole_expansion(comm)
            call build_interaction_lists_fmm
            if (present(damp_type)) then
                if (damp_type /= '') then
                    call particle_fields(comm, field, exclusions, field_order, damp_type, damping_factors)
                else
                    call particle_fields(comm, field, exclusions, field_order)
                end if
            else
                call particle_fields(comm, field, exclusions, field_order)
            end if
            call multipole_to_local(comm)
            call local_to_local
            call local_to_particle(comm, field, field_order)
            call finalize
        end if
    end subroutine field_fmm
end module octree
