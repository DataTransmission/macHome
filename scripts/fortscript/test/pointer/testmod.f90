
module testmod 
    type physics_state
        integer :: &
            lchnk
        real, dimension(:), allocatable    :: &
            lat, lon
        real, dimension(:,:,:), allocatable :: &
            t, u, v
    end type
    interface alloc
        module procedure physics_type_alloc
    end interface alloc
    contains
    subroutine physics_type_alloc(phys_state, state)
        type(physics_state), pointer :: phys_state(:)
        type(physics_state), pointer :: state
        
        allocate(phys_state(1:5),state)
        allocate(state%lat(3))
    end subroutine physics_type_alloc
end module testmod
