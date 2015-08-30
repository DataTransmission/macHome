
      module point_mod
      implicit none
          !private
          type buffer_field_int
              !private
              integer, pointer :: data => null()
          end type buffer_field_int
      interface buffer
          module procedure buff
      end interface buffer
      public buffer_field_int
      contains
      subroutine buff(iptr)
          type(buffer_field_int) :: iptr 
          integer, target :: itar = 1
          if (associated(iptr%data)) then
              print*, itar
          end if
          iptr%data => itar
      end subroutine buff
      end module point_mod
