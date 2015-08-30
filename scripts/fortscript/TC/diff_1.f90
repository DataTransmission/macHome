module !1st derivative module
implicit none
contains

subroutine diffx1(ux,u,dx,nx,ny,nz)
implicit none
real*8,  intent(out) :: ux(nx,ny,nz)
real*8,  intent(in)  :: u(nx,ny,nz), dx
integer, intent(in)  :: nx,ny,nz
integer :: i
do k = 1,nz
   ux(1,:,k) = (u(2,:,k) - u(1,:,k))/dx !FD
     do i = 2,nx-1
        ux(i,:,k) = (u(i+1,:,k) - u(i-1,:,k))/(2*dx) !CD
     enddo
   ux(nx,:,k) = (u(nx,:,k) - u(nx-1,:,k))/dx !BD
enddo
return
end subroutine 

subroutine diffy1(vy,v,dy,nx,ny,nz)
implicit none
real*8,  intent(out) :: vy(nx,ny,nz)
real*8,  intent(in)  :: v(nx,ny,nz), dy
integer, intent(in)  :: nx,ny,nz
integer :: i
do k = 1,nz
   vy(:,1,k) = (v(:,2,k) - v(:,1,k))/dy !FD
     do i = 2,ny-1
        vy(i,:,k) = (v(:,i+1,k) - v(:,i-1,k))/(2*dy) !CD 
     enddo
   vy(:,ny,k) = (v(:,ny,k) - v(;,ny-1,k))/dy !BD
enddo
return
end subroutine 

subroutine diffx2(uxx, u, dx, N)
implicit none
integer, intent(in)   :: N
real*8, intent(in)  :: dx, u(N)
real*8, intent(out) :: uxx(N)
real*8 :: dx2
integer :: i
dx2 = dx*dx
uxx(1)   = (-u(N-2) + 16.0*u(N-1) - 30.0*u(1)   + 16.0*u(2) - u(3)) / (12.0*dx2)
uxx(2)   = (-u(N-1) + 16.0*u(1)   - 30.0*u(2)   + 16.0*u(3) - u(4)) / (12.0*dx2)
uxx(N-1) = (-u(N-3) + 16.0*u(N-2) - 30.0*u(N-1) + 16.0*u(N) - u(2)) / (12.0*dx2)
uxx(N)   = (-u(N-2) + 16.0*u(N-1) - 30.0*u(N)   + 16.0*u(2) - u(3)) / (12.0*dx2)

do i = 3,(N-2)
uxx(i)= (-u(i-2) + 16.0*u(i-1) - 30.0*u(i) + 16.0*u(i+1) - u(i+2)) / (12.0*dx2)
enddo

print*,'uxx = [', uxx, ']'
!do i = 1,N
!write(6,*) uxx(i)
!enddo

return
end subroutine diffx2

end module
