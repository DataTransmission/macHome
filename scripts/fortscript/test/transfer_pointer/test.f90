      program main 

      use point_mod 
      implicit none
      type buffer_field_default_type
          !private
          real, pointer :: data => null()
      end type buffer_field_default_type
!      type buffer_field_int    ! this is already defined in point_mod
          !private
!          integer, pointer :: data => null()
!      end type buffer_field_int
      
      !character(len=10) :: tom
      type(buffer_field_default_type) :: r
      real, target :: r2 
      type(buffer_field_int) :: i
   
      ! r, i are both pointers, r2 is a target
      call buff(i)
      if( associated(i%data)) then
          print*, 'i is associated, and the value is',i%data
      end if
      !r2=1.11
      !r=>r2 ! associate the pointer r to the target r2
      r=transfer(i,r) ! transfer the bit value of i to r and use the type (real) of r
      print*,r
      ! r is associated with i here, i is treated as a target that
      ! changes value with r in a "bit-valued" way
      r%data=1.11
      print*, i%data, r%data, r2 ! if r%data is changed, then i%data changes in a bit-valued way
      print*, transfer(i%data,r%data) ! check whether r2 actually saves the bit values of i
      !tom=transfer(r,tom) 
      !print*,r,tom
      !tom=transfer(i,tom) 
      !print*,i,tom
      !r2 = transfer(i,r2)
      !print*,i,r2
      end program main



 
