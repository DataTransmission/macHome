program vortexinit 
!module vortexinit
!use kinds
use diff

implicit none

real :: a, omg, rd, g, zt, qt, zq1, zq2
real :: gam, dp, rp, zp, phic, lamc, fc, e, eps, rad
real :: d1, d2, d
real :: r, ps, dpdz, f, tcradius, ztc1, vttc !tangential velocity
real :: dlam, dphi
real, allocatable, dimension(:,:,:) :: xxx, t, logp0 
real, allocatable, dimension(:,:)  :: q0b,t0b, dtdp
real, allocatable, dimension(:,:)  :: p0, q0, t0, tv0, tvt, pt
real, allocatable, dimension(:)    :: pmod, lam1, phi1, lam, phi
real, allocatable, dimension(:,:,:):: ztc, plev, qtc, ttc, utc, vtc, pstc 
real, allocatable, dimension(:,:,:):: divtc, vortc, ux, uy, vx, vy 
!real(r8),allocatable :: z(:,:,:), d(:,:,:), t(:,:,:), w(:,:,:)
!real(r8),allocatable :: ztcini(:,:,:), plevini(:,:,:)
integer :: i, j, k, ii, jj
integer :: tcxgrid, tcygrid
integer :: xc, yc, cond, ncond, irec
integer :: nx, ny, nz 
real, parameter :: pi =  3.141592653589793
!In a p-coor model, we use the algorithm to calc the background surface pressure
!distribution do an initial TC. Then by using this ps(r), we could get the z distribution
!at each pmod(k) level. Since the analytical conversion of p to z is inaccurate in the TC,
!(due to warm core and pressure perturbation), so need to use Newton's method to make 
!the plev(k) converge to pmod(k).

!============READIN BCKGRD VARIABLES=======================================
!==========================================================================
tcxgrid=5                                !the x-dir/y-dir +-grids you have in your domain 
tcygrid=2
nx=48
ny=40
nz=5
phic=-15.0                           !phic(lamc)=center lat(lon) of initial vortex, 
lamc=150.0
irec=0
close(111)
open(111,file='t_lastday.dat',form='unformatted',access='direct',recl=nx*ny)
allocate(t(nx,ny,nz),logp0(nx,ny,1)) !z(nx,ny,nz),d(nx,ny,nz),t(nx,ny,nz),w(nx,ny,nz),logp0(nx,ny,1))
!!!!!!!!!!!!INDEX SHOULD BE CONSISTENT WITH OTHER DIMS OF VAR
!call readdat(xxx   ,nx,ny,nz,irec,0) !could be calc from grads and output binary
!call readdat(xxx   ,nx,ny,nz,irec,0) !0 = noread; 1 = read
call readdat(t,nx,ny,nz,irec,1)
call readdat(logp0,nx,ny,1 ,irec,1) !this is called q in Ben's model which is logp0 = sqrt(2)*log(p0)
close(111)
!call readdat(q0b,nx,ny,1,irec)
!call readdat(t0b,nx,ny,1,irec)
!=====suppose to be from model=======
allocate(pmod(nz))
allocate(q0b(nx,ny),t0b(nx,ny),dtdp(nx,ny))
pmod = (/ 900.0, 700.0, 500.0, 300.0, 100.0 /)
q0b=0.021                               !q0=sfc spec hum, 
t0b=302.15 
!do i=1,nx
!   do j=1,ny
!      dtdp(i,j)=(log(t(i,j,2)) - log(t(i,j,1)))/(log(pmod(2))-log(pmod(1))) !
!      t0b (i,j) = exp(log(t(i,j,1)) - dtdp(i,j)*(log(pmod(1)) - 1000*logp0(i,j,1)/2**.5))
!print*,log(pmod(2)), 1000*logp0(i,j,1)/2**.5
!   enddo
!enddo
!print*, t0b
!stop
!==========================================================================
!==========================Genesis Zone====================================
!==========================================================================
allocate( lam1(nx), phi1(ny) ) 
allocate( lam(2*tcxgrid+1), phi(2*tcygrid+1) )
allocate( p0(size(lam),size(phi)),q0(size(lam),size(phi)),t0(size(lam),size(phi)) )
rad = dble(pi/180)
dlam = 360.0/dble(nx-1)
dphi = 180.0/dble(ny-1)
do i = 1,nx                              !i grids do hi res
   lam1(i) = (i-1)*dlam
   if (abs(lam1(i) - lamc) < dlam) then
     xc = i
   endif
enddo
do i = 1,ny 
   phi1(i) = -90.0 + (i-1)*dphi
   if ( abs(phi1(i) - phic) < dphi) then
      yc = i
   endif
enddo
ii=1
do i = xc-tcxgrid,xc+tcxgrid
   if (i <= 0) then
      lam(ii) = lam1(i+nx)*rad
   else
      lam(ii) = lam1(i)*rad
   endif
ii = ii+1
enddo
      phi = phi1(yc-tcygrid:yc+tcygrid)*rad
ii = 1
do i = xc-tcxgrid,xc+tcxgrid
   if (i <= 0) then
      p0(ii,1:2*tcygrid+1)=1000.0*exp(logp0(i+nx,yc-tcygrid:yc+tcygrid,1)/2.0**(0.5)) ! sqrt(2) is a normalize factor, 
      q0(ii,1:2*tcygrid+1)=           q0b(i+nx,yc-tcygrid:yc+tcygrid)
      t0(ii,1:2*tcygrid+1)=           t0b(i+nx,yc-tcygrid:yc+tcygrid)
   else
      p0(ii,1:2*tcygrid+1)=1000.0*exp(logp0(i,yc-tcygrid:yc+tcygrid,1)/2.0**(0.5)) ! sqrt(2) is a normalize factor, 
      q0(ii,1:2*tcygrid+1)=           q0b(i,yc-tcygrid:yc+tcygrid)
      t0(ii,1:2*tcygrid+1)=           t0b(i,yc-tcygrid:yc+tcygrid)
   endif
ii=ii+1
enddo
print*,phic,lamc!,phi/rad,lam/rad
phic=phic*rad                        !phic(lamc)=center lat(lon) of initial vortex, 
lamc=lamc*rad
!====================================================================
!=====================END OF GENESIS ZONE============================
!====================================================================

   
!restrict sfc var to the area of TC initialize zone
!p0(ii,jj)=1015                        !p0=bckgrd sfc p, 
!==========================================================================
!============PARAMETERS====================================================
!==========================================================================
a=6.37122D6                !a=radius of earth, 
omg=7.292115D-5           !omg=rot speed of earth, 
rd=287.04                       !rd=dry air gas ctc, 
g=9.80616                       !g=gravity, 
zt=1.5D4                       !zt=tropopaus ht,
e= 1D-30                       !e=small ctc to avoid division by zero, 
eps=2D-13                   !eps=convergence limit do fixed point iterations
qt=1D-8*1D-3                !qt=upper atm spec hum, 
zq1=3D3                        !zq1 zq2=ctc for spec hum profile, 
zq2=8D3
gam=0.007                       !gam=tv lapse rate, 
allocate(tv0(size(lam),size(phi)),  tvt(size(lam),size(phi)),  pt(size(lam),size(phi)))
do i = 1,size(lam)
   do j = 1,size(phi)
      tv0(i,j)=t0(i,j)*(1.0+0.608*q0(i,j))          !tv0=sfc bckgrd tv(virtual temp), 
      tvt(i,j)=tv0(i,j) - gam*zt                 !tvt=upper atm t set to a ctc,
      pt(i,j) =p0(i,j)*((tvt(i,j)/tv0(i,j))**(g/(rd*gam)))
   enddo
enddo
dp=11.15                       !dp=p0 - ps, sfc p diff btwn bckgrd sfc p0 and ps at ini vortex center, 
rp=282D3                    !rp & zp=ctc for p fit, 
zp=7D3
!------specify genesis domain---------
fc=2.0*omg*sin(phic)            !fc=coriolis param, 
ncond=2
tcradius = 1D8
!==========================================================================
!=============End of Parameter=============================================
!==========================================================================

allocate(utc(size(lam),size(phi),nz),   vtc(size(lam),size(phi),nz), &
         ztc(size(lam),size(phi),nz),   ttc(size(lam),size(phi),nz), &
         qtc(size(lam),size(phi),nz),  pstc(size(lam),size(phi),1),  &
         plev(size(lam),size(phi),nz) )
!==========================================================================
!================CALCULATE VAR=============================================
!==========================================================================
do i=1,size(lam) !1:15
      do j=1,size(phi) !1:8
         r  = a*acos(sin(phic)*sin(phi(j)) + cos(phic)*cos(phi(j))*cos(lam(i)-lamc))
         d1 = sin(phic)*cos(phi(j)) - cos(phic)*sin(phi(j))*cos(lam(i) - lamc)
         d2 = cos(phic)*sin(lam(i) - lamc)                                             !opposite sign in SH
         d  = max(e, (d1**(2.0) + d2**(2.0))**(0.5))
         ps = p0(i,j) - dp*exp(-(r/rp)**(1.5))  !ps is the sfc-p to initial TC, p0 is the surface p from the model output
         do k=1,nz
            if (pmod(k) .GE. pt(i,j) .AND. ps .GE. pmod(k)) then 
              cond = 1
              ztc(i,j,k) = tv0(i,j)/gam*(1.0 - (pmod(k)/ps)**(rd*gam/g)) !ztc at each pmod level 
!              ztcini(i,j,k)=ztc(i,j,k);
!              plevini(i,j,k) = (p0(i,j) - dp*exp(-(r/rp)**1.5)*exp(-(ztc(i,j,k)/zp)**2))*((tv0 -  &
!                      gam*ztc(i,j,k))/tv0)**(g/(rd*gam)); ! pressure coor + perturbation 
              if ( r .LE. tcradius ) then
                  do while (cond .LT. ncond)
                    plev(i,j,k) = (p0(i,j) - dp*exp(-(r/rp)**(1.5))*exp(-(ztc(i,j,k)/zp)**(2.0)))*((tv0(i,j) -  &
                           gam*ztc(i,j,k))/tv0(i,j))**(g/(rd*gam)) ! pressure coor + perturbation 
                    f = pmod(k) - plev(i,j,k)
                    dpdz = 2.0*dp*ztc(i,j,k)/zp**(2.0)*exp(-(r/rp)**(1.5))*exp(-(ztc(i,j,k)/zp)**(2.0))* &
                          ((tv0(i,j) - gam*ztc(i,j,k))/tv0(i,j))**(g/(rd*gam)) - g/(rd*tv0(i,j))* &
                          (p0(i,j) - dp*exp(-(r/rp)**(1.5))*exp(-(ztc(i,j,k)/zp)**(2.0)))* &
                          ((tv0(i,j) - gam*ztc(i,j,k))/tv0(i,j))**(g/(rd*gam) - 1.0)
                    ztc1 = ztc(i,j,k)
                    ztc(i,j,k)  = ztc(i,j,k) + f/dpdz
                    if (ztc(i,j,k) - ztc1 .LE. eps) then 
                      ztc(i,j,k) = ztc1
                      cond = ncond
                    endif 
                  enddo
              endif
              qtc(i,j,k) = q0(i,j)*exp((-ztc(i,j,k)/zq1))*exp(-(ztc(i,j,k)/zq2)**(2.0)) !q(z) dist for vortex
                                  !Refer to http://www.applet-magic.com/gradientwind.htm for Gradient wind Balance dep in Hemisphere
if (fc .LT. 0) then
r = -r
endif
              vttc = -abs(fc*r*0.5) + ( 0.25*fc*fc*r*r + abs( - & 
                     ( 1.5*(r/rp)**(1.5)*rd*(tv0(i,j)-gam*ztc(i,j,k)) ) / &
                     ( 1.0 + 2.0*rd*(tv0(i,j) - gam*ztc(i,j,k)) * ztc(i,j,k) / & 
                     ( g*zp*zp ) - p0(i,j)/dp*exp( (abs(r)/rp)**(1.5) ) * exp( ( ztc(i,j,k)/zp )**(2.0) ) ) ))**(0.5)
r = -r
              ttc(i,j,k) = (tv0(i,j)-gam*ztc(i,j,k))/(1.0+0.608*qtc(i,j,k))*(1.0+(2.0*rd*(tv0(i,j)-gam*ztc(i,j,k))*ztc(i,j,k))/  &
                          (g*zp**(2.0)*(1.0-p0(i,j)/dp*exp((r/rp)**(1.5))*exp((ztc(i,j,k)/zp)**(2.0)))))**(-1.0)
              !t(phi,lam,z) dist do vortex and bckgrd, t approach bckgrd temp as r becomes large  
            else
              qtc(i,j,k) = qt       !dry model doesn't need it
              ztc(i,j,k) = zt + rd*tvt(i,j)/g*log(pt(i,j)/pmod(k))
              plev(i,j,k) = pt(i,j)*exp(-((g*(ztc(i,j,k)-zt))/(rd*tvt(i,j))))
              ttc(i,j,k) = tvt(i,j)
              vttc = 0.0         !tangential wind
            endif
               utc(i,j,k) = vttc*d1/d                   !CHANGE THE CALC of d1, d2 in SH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               vtc(i,j,k) = vttc*d2/d
         enddo
            pstc(i,j,1)= pmod(1)/(1.0-ztc(i,j,1)*gam/tv0(i,j))**(g/(rd*gam)) ! from inversing 
      enddo
enddo
  allocate(divtc(size(lam),size(phi),nz),vortc(size(lam),size(phi),nz))
  allocate(ux(size(lam),size(phi),nz),uy(size(lam),size(phi),nz),vx(size(lam),size(phi),nz),vy(size(lam),size(phi),nz))
  call diffx1(ux,utc,dlam,size(lam),size(phi),nz)
  call diffx1(vx,vtc,dlam,size(lam),size(phi),nz)
  call diffy1(uy,vtc,dphi,size(lam),size(phi),nz)
  call diffy1(vy,utc,dphi,size(lam),size(phi),nz)
  divtc = ux + vy
  vortc = vx - uy 
open(unit=110,file='plev.dat',form='formatted',status='unknown',action='write')
open(unit=111,file='pstc.dat',form='formatted',status='unknown',action='write')
open(unit=112,file='utc.dat',form='formatted',status='unknown',action='write')
open(unit=113,file='vtc.dat',form='formatted',status='unknown',action='write')
open(unit=114,file='ztc.dat',form='formatted',status='unknown',action='write')
open(unit=115,file='ttc.dat',form='formatted',status='unknown',action='write')
open(unit=116,file='divtc.dat',form='formatted',status='unknown',action='write')
open(unit=117,file='vortc.dat',form='formatted',status='unknown',action='write')

      do k=1,nz
         write(110,'(11f15.5)') plev(:,:,k) !u
         write(111,'(11f15.5)') pstc(:,:,1) !u
         write(112,'(11f15.5)') utc(:,:,k) !u
         write(113,'(11f15.5)') vtc(:,:,k) !v
         write(114,'(11f15.5)') ztc(:,:,k) !div
         write(115,'(11f15.5)') ttc(:,:,k) !temp
         write(116,'(11f15.5)') divtc(:,:,k) !geopot
         write(117,'(11f15.5)') vortc(:,:,k) 
      enddo
  close(110)
  close(111)
  close(112)
  close(113)
  close(114)
  close(115)
  close(116)
  close(117)
!return
stop
!end subroutine initialize

contains

subroutine readdat(varout,nx,ny,nz,irec,noread)
implicit none
integer,intent(in) :: nx, ny, nz 
integer,intent(in) :: noread
integer :: k
integer, intent(inout) :: irec
real, intent(out) :: varout(nx,ny,nz) !!!!!!!!!!!!!there might be some problem with the ny, nx order 
!allocate(varout(nx,ny,nz))
do k=1,nz
   irec = irec+1
!   if (noread == 0) then
!      read(111,rec=irec)
!   else
      read(111,rec=irec)varout(:,:,k)  
!   endif
enddo
return
end subroutine readdat


!end module
end program
