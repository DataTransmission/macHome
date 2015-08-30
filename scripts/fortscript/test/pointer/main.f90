

program main
! See how the subroutine testmod with two different structure variables allocated (1d
! structure array allocated and a single structure variable allocated)
! were imported into the main program. 
    use testmod
    implicit none
    type(physics_state), pointer :: phys_state(:) => null() ! define a dissociated pointer
    type(physics_state), pointer :: state => null() 
    type(physics_state), target  :: state_tgt
    call physics_type_alloc(phys_state, state) 
    state_tgt%lchnk=1.1
    allocate(state_tgt%lat(3)) ! need to allocate here sine state_tgt
                               ! wasn't called from physics_type_alloc
    state_tgt%lat=(/3,2,1/)
    state => state_tgt
    phys_state(:)%lchnk=2.2

    print*, state
    print*, phys_state
end program main
