module findminval
implicit none
! compare whether the centered value is the minval in a boxarea 
! if not then see if it's the second minval in the boxarea, if not
! then see if it's the third minval in the box area and so on. 
! if runned over nmin times, the centered value still isn't the minval,
! then condz will equal nmin, else it will be smaller than nmin
real,allocatable :: zvmat(:), zv(:,:,:,:,:)
integer :: ni,nj,nz,nii,nmin,condz,boxarea,hfbox
integer :: i,j,k,t,e,r
ni=;nj=;nz=;nmin=;hfbox=3
boxarea=(hfbox*2+1)**2
allocate(zvmat(boxarea))
condz = 0

do e=
do t=
do i=
do j=

do k=1,nz
  nii=0
!======save the boxarea matrix in zvmat(:)=======
  do ni=i-hfbox,i+hfbox
   do nj=j-hfbox,j+hfbox
       nii=nii+1
       zvmat(nii) = zv(ni,nj,k,t,e)
    enddo
  enddo

!=====compare nmin times whether the centered value is a r-th minval in the boxarea========
  do r=1,nmin
    !write(*,'(2f15.8)'),zv(i,j,1,t,e),minval(zvmat)
    if (minval(zvmat) == zv(i,j,k,t,e)) then
    else
      call minimum(zvmat,boxarea,condz)!find the nmin-th minimum value for minimum zv value
    endif
  enddo
enddo
if ( condz < nz*nmin ) then
!detemine what variable to save here
endif

enddo
enddo
enddo
enddo

contains

subroutine minimum(zvmat,boxarea,condz)
implicit none          
integer,intent(inout) :: condz, boxarea
real,intent(inout) :: zvmat(boxarea)
integer :: ni
real :: minzvmat
minzvmat=minval(zvmat)
do ni=1,boxarea
  if (zvmat(ni) /= minzvmat) then
  else
    zvmat(ni) = 10e5 
  endif
enddo
condz = condz + 1
return
end subroutine minimum
end module 
