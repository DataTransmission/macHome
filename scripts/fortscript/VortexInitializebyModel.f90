program vortexinit 
!module vortexinit
!use kinds
use diff

implicit none

real :: a, omg, rd, g, zt, qt, zq1, zq2
real :: gam, dp, rp, zp, phic, lamc, fc, e, eps, rad
real :: d1, d2, d
real :: r, ps, dpdz, f, tcradius, ztc1 !tangential velocity
real, allocatable, dimension(:,:,:) :: xxx, tb, logp0, vortb, divb, tv 
real, allocatable, dimension(:,:)  :: q0b,t0b, dtdp
real, allocatable, dimension(:,:)  :: p0, q0, t0, tv0, tvt, pt, psb
real, allocatable, dimension(:)    :: lam1, phi1, lam, phi
real, allocatable, dimension(:,:,:):: zb,ztc,ztcini, plev, qtc, ttc, utc, vtc, pstc, vttc 
real, allocatable, dimension(:,:,:):: divtc, vortc, ux, uy, vx, vy 
!real(r8),allocatable :: z(:,:,:), d(:,:,:), t(:,:,:), w(:,:,:)
!real(r8),allocatable :: ztcini(:,:,:), plevini(:,:,:)
integer :: i, j, k, ii, jj
integer :: tcxgrid, tcygrid
integer :: xc, yc, cond, ncond, irec
integer,parameter :: nx=48, ny=40, nz=5 
real :: dlam, tcdx!, dphi, tcdy
real, dimension(ny-1) :: dphi !for Gauss Lat
real, allocatable, dimension(:) :: tcdy
real, dimension(nx,ny,nz) :: pmod
real, parameter, dimension(nz) :: sigma = (/ 0.9, 0.7, 0.5, 0.3, 0.1 /) !sigma = pmod/ps 
real, parameter :: pi =  3.141592653589793
character*20 :: xsize, lamsize 
!In a p-coor model, we use the algorithm to calc the background surface pressure
!distribution for an initial TC. Then by using this ps(r), we could get the z distribution
!at each pmod(k) level. Since the analytical conversion of p to z is inaccurate in the TC,
!(due to warm core and pressure perturbation), so need to use Newton's method to make 
!the plev(k) converge to pmod(k).

!==========================================================================
!============PARAMETERS====================================================
!==========================================================================
allocate(q0b(nx,ny),t0b(nx,ny),dtdp(nx,ny))
tcxgrid=3                                !the x-dir/y-dir +-grids you have in your domain 
tcygrid=3
phic=10.0                           !phic(lamc)=center lat(lon) of initial vortex, 
lamc=150.0
q0b=0.021    !use 0.021 in paper                              !q0=sfc spec hum, 
!t0b=302.15               
a=6.37122E6                !a=radius of earth, 
omg=7.292115E-5           !omg=rot speed of earth, 
rd=287.04                       !rd=dry air gas ctc, 
g=9.80616                       !g=gravity, 
zt=1.5E4                       !zt=tropopaus ht,
e= 1E-30                       !e=small ctc to avoid division by zero, 
eps=2E-13                   !eps=convergence limit do fixed point iterations
qt=1E-8*1E-3                !qt=upper atm spec hum, 
zq1=3E3                        !zq1 zq2=ctc for spec hum profile, 
zq2=8E3
gam=0.0055                      !gam=tv lapse rate, 7 degrees per km 
!dp=11.15                       !dp=p0 - ps, sfc p diff btwn bckgrd sfc p0 and ps at ini vortex center, 
dp = 20   !artificially changed, NEED to CHANGE BACK!!!!!!!!!!!!!!!!!!!!!!!!!!!
rp=282E3                    !rp & zp=ctc for p fit, 
zp=7E3
tcradius = 1E9
rad = dble(pi/180)
dlam = 360.0/dble(nx)  !7.5 in dryGCM
allocate(phi1(ny))
phi1 = (/-89.841393873950310, -89.165361482951127, -87.953395498539663, -86.212513729241280, -83.953152745080843,& 
         -81.188892627198655, -77.936355289103304, -74.215100774997993, -70.047508628386709, -65.458642967093425,&
         -60.476101615276114, -55.129850070118195, -49.452041258561515, -43.476822151756018, -37.240128393444451,& 
         -30.779468174318218, -24.133696650652773, -17.342782263123389, -10.447566360772949,  -3.489517575544590,&
           3.489517575544546,  10.447566360772949,  17.342782263123389,  24.133696650652773,  30.779468174318218,&  
          37.240128393444451,  43.476822151756018,  49.452041258561515,  55.129850070118195,  60.476101615276114,& 
          65.458642967093425,  70.047508628386709,  74.215100774997993,  77.936355289103304,  81.188892627198655,&  
          83.953152745080843,  86.212513729241280,  87.953395498539663,  89.165361482951127,  89.841393873950310/)
do i = 1,ny-1
   dphi(i) = phi1(i+1) - phi1(i)  !only going to be ny-1 dphi's for Guass lat 
enddo
!dphi = 180.0/dble(ny-1)!Gauss-Lat

fc=2.0*omg*sin(phic*rad)            !fc=coriolis param, 
ncond=2
write(xsize,'(i10)') tcxgrid*2+1
irec=0
open(111,file='t_lastday.dat',form='unformatted',access='direct',recl=nx*ny)
allocate(xxx(nx,ny,nz),vortb(nx,ny,nz),divb(nx,ny,nz),tb(nx,ny,nz),logp0(nx,ny,1),psb(nx,ny)) !z(nx,ny,nz),d(nx,ny,nz),t(nx,ny,nz),w(nx,ny,nz),logp0(nx,ny,1))
!!!!!!!!!!!!INDEX SHOULD BE CONSISTENT WITH OTHER DIMS OF VAR

call readdat(vortb,nx,ny,nz,irec,1) !could be calc from grads and output binary
call readdat(divb ,nx,ny,nz,irec,1) !0 = noread; 1 = read
call readdat(tb   ,nx,ny,nz,irec,1)
!call readdat(xxx  ,nx,ny,nz,irec,0) ! skip diabatic heating
call readdat(logp0,nx,ny,1 ,irec,1) !this is called q in Ben's model which is logp0 = sqrt(2)*log(p0)
psb = 1000.0*exp(logp0(:,:,1)/2.0**(0.5))

do k = 1,nz
   pmod(:,:,k) = sigma(k)*psb(:,:)
enddo
     !these fields except tb and logp0 should be overlapped as anomolous field, tb and logp0 should be cropped and replaced
    !call readdat(q0b,nx,ny,1,irec) these two field might be used during CAM
    !call readdat(t0b,nx,ny,1,irec)
close(111)
    !open(unit=99,file='tb.dat',form='formatted',status='unknown',action='write')
    !do k=1,nz
    !   write(99,'(48f15.5)') t(:,:,k) !u
    !enddo
!=====suppose to be from CAM model=======
do i=1,nx
   do j=1,ny
      dtdp(i,j)=(log(tb(i,j,2)) - log(tb(i,j,1)))/(log(pmod(i,j,2))-log(pmod(i,j,1))) !
      t0b (i,j) = exp(log(tb(i,j,1)) - dtdp(i,j)*(log(pmod(i,j,1)) - (log(1000.0) + logp0(i,j,1)/2.0**0.5)))
                                     !log(1000) + logp0(i,j,1)/2**.t = model ln(ps)
   enddo
enddo
!==========================================================================
!==========================Genesis Zone====================================
!==========================================================================
allocate( lam1(nx) ) 
allocate( lam(2*tcxgrid+1), phi(2*tcygrid+1) )
allocate( p0(size(lam),size(phi)),q0(size(lam),size(phi)),t0(size(lam),size(phi)) )
do i = 1,nx                              !i grids do hi res
   lam1(i) = (i-1)*dlam
   if (abs(lam1(i) - lamc) < dlam) then
     xc = i
   endif
enddo
do i = 1,ny-1 !ny-1 b/c we only have ny-1 dphi's to compare 
   if ( abs(phi1(i) - phic) < dphi(i)) then
!   phi1(i) = -90.0 + (i-1)*dphi         !for non gauss latitude model
!   if ( abs(phi1(i) - phic) < dphi) then
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
phic=phic*rad                        !phic(lamc)=center lat(lon) of initial vortex, 
lamc=lamc*rad
!====================================================================
!=====================END OF GENESIS ZONE============================
!====================================================================

   
      !restrict sfc var to the area of TC initialize zone
      !p0(ii,jj)=1015                        !p0=bckgrd sfc p, 
allocate(tv0(size(lam),size(phi)),  tvt(size(lam),size(phi)),  pt(size(lam),size(phi)))
do i = 1,size(lam)
   do j = 1,size(phi)
      tv0(i,j)=t0(i,j)*(1.0+0.608*q0(i,j))          !tv0=sfc bckgrd tv(virtual temp), 
      tvt(i,j)=tv0(i,j) - gam*zt                 !tvt=upper atm t set to a ctc,
      pt(i,j) =p0(i,j)*((tvt(i,j)/tv0(i,j))**(g/(rd*gam)))
   enddo
enddo
!==========================================================================
!=============End of Parameter=============================================
!==========================================================================

allocate(utc(size(lam),size(phi),nz),   vtc(size(lam),size(phi),nz), &
         ztc(size(lam),size(phi),nz),   ttc(size(lam),size(phi),nz), &
         qtc(size(lam),size(phi),nz),  pstc(size(lam),size(phi),1),  &
         plev(size(lam),size(phi),nz), vttc(size(lam),size(phi),nz), &
         ztcini(size(lam),size(phi),nz), tv(size(lam),size(phi),nz))
allocate(zb(size(lam),size(phi),nz))
!==========================================================================
!================CALCULATE VAR=============================================
!==========================================================================
print*,gam
do i=1,size(lam) !1:15
      do j=1,size(phi) !1:8
         r  = a*acos(sin(phic)*sin(phi(j)) + cos(phic)*cos(phi(j))*cos(lam(i)-lamc))
         d1 = sin(phic)*cos(phi(j)) - cos(phic)*sin(phi(j))*cos(lam(i) - lamc)
         d2 = cos(phic)*sin(lam(i) - lamc)                                             !opposite sign in SH
         d  = max(e, (d1**(2.0) + d2**(2.0))**(0.5))
         ps = p0(i,j) - dp*exp(-(r/rp)**(1.5))  !ps is the sfc-p to initial TC, p0 is the surface p from the model output
         do k=1,nz
            zb(i,j,k) =tv0(i,j)/gam*(1.0 - (pmod(i,j,k)/p0(i,j))**(rd*gam/g)) !zb at each pmod level 
            if (pmod(i,j,k) .GE. pt(i,j) .AND. ps .GE. pmod(i,j,k)) then 
              cond = 1
              ztc(i,j,k) = tv0(i,j)/gam*(1.0 - (pmod(i,j,k)/ps)**(rd*gam/g)) !ztc at each pmod level 
!              ztcini(i,j,k)=ztc(i,j,k);
!!              plevini(i,j,k) = (p0(i,j) - dp*exp(-(r/rp)**1.5)*exp(-(ztc(i,j,k)/zp)**2))*((tv0 -  &
!!                      gam*ztc(i,j,k))/tv0)**(g/(rd*gam)); ! pressure coor + perturbation 
!              if ( r .LE. tcradius ) then
!                do while (cond .LT. ncond)
!                    plev(i,j,k) = (p0(i,j) - dp*exp(-(r/rp)**(1.5))*exp(-(ztc(i,j,k)/zp)**(2.0)))*((tv0(i,j) -  &
!                           gam*ztc(i,j,k))/tv0(i,j))**(g/(rd*gam)) ! pressure coor + perturbation 
!                    f = pmod(i,j,k) - plev(i,j,k)
!                    dpdz = 2.0*dp*ztc(i,j,k)/zp**(2.0)*exp(-(r/rp)**(1.5))*exp(-(ztc(i,j,k)/zp)**(2.0))* &
!                          ((tv0(i,j) - gam*ztc(i,j,k))/tv0(i,j))**(g/(rd*gam)) - g/(rd*tv0(i,j))* &
!                          (p0(i,j) - dp*exp(-(r/rp)**(1.5))*exp(-(ztc(i,j,k)/zp)**(2.0)))* &
!                          ((tv0(i,j) - gam*ztc(i,j,k))/tv0(i,j))**(g/(rd*gam) - 1.0)
!                    ztc1 = ztc(i,j,k)
!                    ztc(i,j,k)  = ztc(i,j,k) + f/dpdz
!                    if (ztc(i,j,k) - ztc1 .LE. eps) then 
!                      ztc(i,j,k) = ztc1
!                      cond = ncond
!                    endif 
!                 enddo
!              endif
               !this is the original 
              qtc(i,j,k) = q0(i,j)*exp((-ztc(i,j,k)/zq1))*exp(-(ztc(i,j,k)/zq2)**(2.0)) !q(z) dist for vortex
                                  !Refer to http://www.applet-magic.com/gradientwind.htm for Gradient wind Balance dep in Hemisphere
              vttc(i,j,k) = -abs(fc*r*0.5) + ( 0.25*fc*fc*r*r -  &            !THIS WILL CREATE PROBLEM IN SH 
                     ( 1.5*(r/rp)**(1.5)*rd*(tv0(i,j)-gam*ztc(i,j,k)) ) / &
                     ( 1.0 + 2.0*rd*(tv0(i,j) - gam*ztc(i,j,k)) * ztc(i,j,k) / & 
                     ( g*zp*zp ) - p0(i,j)/dp*exp( (r/rp)**(1.5) ) * exp( ( ztc(i,j,k)/zp )**(2.0) ) ) )**(0.5)
!              ttc(i,j,k) = (tv0(i,j)-gam*ztc(i,j,k))/(1.0+0.608*qtc(i,j,k))*(1.0+(2.0*rd*(tv0(i,j)-gam*ztc(i,j,k))*ztc(i,j,k))/  &
!                          (g*zp**(2.0)*(1.0-p0(i,j)/dp*exp((r/rp)**(1.5))*exp((ztc(i,j,k)/zp)**(2.0)))))**(-1.0)
               !this is using ztc instead zb 
!              qtc(i,j,k) = q0(i,j)*exp((-zb(i,j,k)/zq1))*exp(-(zb(i,j,k)/zq2)**(2.0)) !q(z) dist for vortex
!              vttc(i,j,k) = -abs(fc*r*0.5) + ( 0.25*fc*fc*r*r -  &            !THIS WILL CREATE PROBLEM IN SH 
!                     ( 1.5*(r/rp)**(1.5)*rd*(tv0(i,j)-gam*zb(i,j,k)) ) / &
!                     ( 1.0 + 2.0*rd*(tv0(i,j) - gam*zb(i,j,k)) * zb(i,j,k) / & 
!                     ( g*zp*zp ) - p0(i,j)/dp*exp( (r/rp)**(1.5) ) * exp( ( zb(i,j,k)/zp )**(2.0) ) ) )**(0.5)
!              ttc(i,j,k) = (tv0(i,j)-gam*zb(i,j,k))/(1.0+0.608*qtc(i,j,k))*(1.0+(2.0*rd*(tv0(i,j)-gam*zb(i,j,k))*zb(i,j,k))/  &
!                          (g*zp**(2.0)*(1.0-p0(i,j)/dp*exp((r/rp)**(1.5))*exp((zb(i,j,k)/zp)**(2.0)))))**(-1.0)
                             !t(phi,lam,z) dist for vortex and bckgrd, t approach bckgrd temp as r becomes large  
               !this is using t'
              ttc(i,j,k) = (tv0(i,j) - gam*zb(i,j,k))/(1.0+0.608*qtc(i,j,k))*(-2.0*rd*(tv0(i,j)-gam*zb(i,j,k))*zb(i,j,k))/&
                           (2.0*rd*(tv0(i,j)-gam*zb(i,j,k))*zb(i,j,k) + g*zp**(2.0)*(1.0-p0(i,j)/dp*exp((r/rp))**(1.5)*exp((zb(i,j,k)/zp)**(2.0))))
            else
              qtc(i,j,k) = qt       !dry model doesn't need it
              ztc(i,j,k) = zt + rd*tvt(i,j)/g*log(pt(i,j)/pmod(i,j,k))
              plev(i,j,k) = pt(i,j)*exp(-((g*(ztc(i,j,k)-zt))/(rd*tvt(i,j))))
              !this is the original
!              ttc(i,j,k) = tvt(i,j)
              ttc(i,j,k) = 0.0
              vttc(i,j,k) = 0.0         !tangential wind
            endif
               utc(i,j,k) = vttc(i,j,k)*d1/d                   !CHANGE THE CALC of d1, d2 in SH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               vtc(i,j,k) = vttc(i,j,k)*d2/d
         enddo
         pstc(i,j,1)= pmod(i,j,1)/(1.0-ztc(i,j,1)*gam/tv0(i,j))**(g/(rd*gam)) ! from inversing, i guess already calc in ps 
    enddo
enddo
  allocate(divtc(size(lam),size(phi),nz),vortc(size(lam),size(phi),nz))
  allocate(ux(size(lam),size(phi),nz),uy(size(lam),size(phi),nz),vx(size(lam),size(phi),nz),vy(size(lam),size(phi),nz))
  allocate(tcdy(size(phi)-1))
!  tcdx = 1.1E6*dlam
  call diffx(ux,utc,size(lam),size(phi),nz,phi,nx)
  call diffx(vx,vtc,size(lam),size(phi),nz,phi,nx)
!  tcdy = 1D6*dphi
!  call diffy(uy,vtc,tcdy,size(lam),size(phi),nz)
!  call diffy(vy,utc,tcdy,size(lam),size(phi),nz)
  tcdy = 1E6*dphi(yc-tcygrid:yc+tcygrid-1)
  call diffy_gauss(uy,vtc,tcdy,size(lam),size(phi),nz)
  call diffy_gauss(vy,utc,tcdy,size(lam),size(phi),nz)
!print*,tcdx,tcdy
!pause
  divtc = ux + vy
  vortc = vx - uy 
if (fc .LT. 0) then
   utc = -utc    !SH wind
   vtc = -vtc
endif
    !overlap vort, div; replace logps, temp
do k = 1,nz
  ii = 1
  do i = xc-tcxgrid,xc+tcxgrid
    jj=1
    do j = yc-tcygrid,yc+tcygrid
       vortb(i,j,k) = vortb(i,j,k) + vortc(ii,jj,k)*100    !!!!!!Artificially 10 mag, need to cancel this once tested
       divb(i,j,k)  = divb(i,j,k) + divtc(ii,jj,k)*100 
       tb(i,j,k) =  tb(i,j,k) + ttc(ii,jj,k)
       if (k==1) then
          logp0(i,j,1) = 2.0**0.5*log(pstc(ii,jj,1)/1000.0)
       endif
       jj=jj+1
    enddo
    ii=ii+1
  enddo
enddo
  !!!FORTRAN
  open(unit=90,file= 'vort_tot.dat',form='unformatted')
  open(unit=91,file=  'div_tot.dat',form='unformatted')
  open(unit=92,file=    't_tot.dat',form='unformatted')
  open(unit=93,file='logp0_tot.dat',form='unformatted')
  !!!MATLAB
!  open(unit=90,file= 'vort_tot.dat',form='formatted',status='unknown',action='write')!,recl=nx*ny)
!  open(unit=91,file=  'div_tot.dat',form='formatted',status='unknown',action='write')!,recl=nx*ny)
!  open(unit=92,file=    't_tot.dat',form='formatted',status='unknown',action='write')!,recl=nx*ny)
!  open(unit=93,file='logp0_tot.dat',form='formatted',status='unknown',action='write')!,recl=nx*ny)
            !if want to use matlab, change to formatted
  open(unit=110,file='plev.dat',form='formatted',status='unknown',action='write')
  open(unit=111,file='pstc.dat',form='formatted',status='unknown',action='write')
  open(unit=112,file='utc.dat',form='formatted',status='unknown',action='write')
  open(unit=113,file='vtc.dat',form='formatted',status='unknown',action='write')
  open(unit=114,file='ztc.dat',form='formatted',status='unknown',action='write')
  open(unit=115,file='ttc.dat',form='formatted',status='unknown',action='write')
  open(unit=116,file='divtc.dat',form='formatted',status='unknown',action='write')
  open(unit=117,file='vortc.dat',form='formatted',status='unknown',action='write')
  open(unit=118,file='vttc.dat',form='formatted',status='unknown',action='write')
  open(unit=119,file='ztcini.dat',form='formatted',status='unknown',action='write')
  open(unit=120,file='lam.dat',form='formatted',status='unknown',action='write')
  open(unit=121,file='phi.dat',form='formatted',status='unknown',action='write')
  write(lamsize,'(i10)') nx
     !!!!FORTRAN
      do k=1,nz
         write(90) vortb(:,:,k) !bg w/ tc
         write(91) divb(:,:,k) !bg w/ tc
         write(92) tb(:,:,k) !bg w/ tc
!         write(93) 1000.0*exp(logp0(:,:,k)/2.0**0.5) !bg w/ tc
         write(93) logp0(:,:,k) !bg w/ tc
      enddo
     !!!!MATLAB
      do k=1,nz
!         write(90,'('// lamsize //'f20.12)') vortb(:,:,k) !bg w/ tc
!         write(91,'('// lamsize //'f20.12)') divb(:,:,k) !bg w/ tc
!         write(92,'('// lamsize //'f20.12)') tb(:,:,k) !bg w/ tc
         write(110,'('// xsize //'f20.12)') plev(:,:,k) !u
         write(111,'('// xsize //'f20.12)') 2**0.5*log(pstc(:,:,1)/1000) !u
         write(112,'('// xsize //'f20.12)') utc(:,:,k) !u
         write(113,'('// xsize //'f20.12)') vtc(:,:,k) !v
         write(114,'('// xsize //'f20.12)') ztc(:,:,k) !div
         write(115,'('// xsize //'f20.12)') ttc(:,:,k) !temp
         write(116,'('// xsize //'f20.12)') divtc(:,:,k) !geopot
         write(117,'('// xsize //'f20.12)') vortc(:,:,k) 
         write(118,'('// xsize //'f20.12)') vttc(:,:,k) 
         write(119,'('// xsize //'f20.12)') ztcini(:,:,k) 
      enddo
!         write(93,'('// lamsize //'f20.12)') logp0(:,:,1) !bg w/ tc
         write(120,'('// xsize //'f20.12)') lam/rad
         write(121,'('// xsize //'f20.12)') phi/rad
  close(90)
  close(91)
  close(92)
  close(93)
  close(110)
  close(111)
  close(112)
  close(113)
  close(114)
  close(115)
  close(116)
  close(117)
  close(118)
  close(119)
  close(120)
  close(121)
!print*, lam,lam1*rad
stop

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
