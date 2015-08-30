program epflux
!===========================================================================
!=========To Solve EPFlux Div and Residual Circulation======================

implicit none
real,allocatable :: test(:,:),test2(:,:,:),temp(:,:,:,:),pottemp(:,:,:,:)
real,allocatable :: lam(:),phi(:),p(:,:,:,:),z(:),plevel(:)
real,allocatable :: u(:,:,:,:),v(:,:,:,:),w(:,:,:,:),vort(:,:,:,:),div(:,:,:,:)
real,allocatable :: rho0(:),ub(:,:,:,:),vb(:,:,:,:),wb(:,:,:,:),phttempb(:,:,:,:)
real,allocatable :: up(:,:,:,:),vp(:,:,:,:),wp(:,:,:,:),pottempp(:,:,:,:)
real,allocatable :: ubz(:,:,:,:),upvpb(:,:,:,:),vppottemppb(:,:,:,:),wpupb(:,:,:,:)
real,allocatable :: epfluxy(:,:,:,:),epfluxz(:,:,:,:)
real,allocatable :: epfluxdivy(:,:,:,:),epfluxdivz(:,:,:,:),epfluxdiv(:,:,:,:)
real,allocatable :: ps(:,:,:),q(:,:,:,:),dz(:)
real,allocatable :: pottempb(:,:,:,:),pottempbz(:,:,:,:),vbpottemppb(:,:,:,:)
real,allocatable :: upvp(:,:,:,:),vppottempp(:,:,:,:),wpup(:,:,:,:)
real,allocatable :: f(:),ubcosphiphi(:,:,:,:)
real,allocatable :: rho0ppb1(:,:,:,:),rho0ppb2(:,:,:,:),rho0ppbz(:,:,:,:) 
real,allocatable :: rho0ppbphi(:,:,:,:),vbstar(:,:,:,:),wbstar(:,:,:,:)
real,allocatable :: vbstaradv(:,:,:,:),wbstaradv(:,:,:,:)
real,allocatable :: epfluxdivyt(:,:,:,:),epfluxdivzt(:,:,:,:),epfluxdivt(:,:,:,:),
real,allocatable :: vbstart(:,:,:,:),wbstart(:,:,:,:),vbstaradvt(:,:,:,:)
real,allocatable :: wbstaradvt(:,:,:,:)
real, parameter  :: Rd = 286.9, cp=1004, omega=7.292*10e-5 ! J kg-1 K-1 
real, parameter  :: kp=Rd/cp, p0=1000 !ref pressure
real, parameter  :: g=9.8,h0=8.2e3 !inverse of scale height 
real, parameter  :: tr=h0*g/rd !scale height ref temp
real, parameter  :: rhor=p0/(Rd*tr) !ref rho at surface 
real, parameter  :: a=6384e3,pi=3.141592653589793238462 !earth radius
real             :: dx,dy
integer          :: nx,ny,nz,nt,irec
integer          :: i,j,k,t
 nx=48
 ny=40
 nz=5
 nt=1900
 irec=0
allocate(z(nz),dz(nz),plevel(nz),rho0(nz))
allocate(phi(ny),ubz(nx,ny,nz,nt),f(ny))
close(111)
 plevel=(/900,700,500,300,100/)
 z(:) = -h0*log(plevel(:)/p0)
 !dx=*2*pi/nx
 dy=2.0012e7/ny
 !500.3 km
 dz(:)=-h0*(-log((plevel(:)+100)/p0) + log(plevel(:)/p0))
 !863.9565       1094.958       1495.037       2358.993       5683.807 
 phi(1)=-0.5*pi 
   do j=2,ny
      phi(j)= pi/ny + phi(j-1)
   enddo
 do j=1,ny
    f(j) = 2*omega*sin(phi(j))
 enddo
!===================OPEN Model Ouput=====================
allocate(vort(nx,ny,nz,nt),div(nx,ny,nz,nt),temp(nx,ny,nz,nt))
allocate(q(nx,ny,nz,nt),u(nx,ny,nz,nt),v(nx,ny,nz,nt),ps(nx,ny,nt))

 open(111,file='t.dat',form='unformatted',access='direct',recl=nx*ny)
   do t=1,nt
    do k=1,nz
       irec = irec+1
       read(111,rec=irec)vort(:,:,k,t)
    enddo
    do k=1,nz
       irec=irec+1
       read(111,rec=irec)div(:,:,k,t) 
       !effected by prescribed forcing/heating and horizontal velocity 
       !reprsenting vertical velocity gradient
    enddo
    do k=1,nz
       irec=irec+1
       read(111,rec=irec)temp(:,:,k,t)
    enddo
!print*,temp(:,:,1,nt)
    do k=1,nz
       irec=irec+1
       read(111,rec=irec)q(:,:,k,t)
    enddo
       irec=irec+1
       read(111,rec=irec)ps(:,:,t)
   enddo
 open(54,file='ps.dat',action='write',form='formatted') 
   do t=1,3
         write(54,'(48f15.5)')ps(:,:,t)
   enddo
 close(54)

!print*,vort(:,:,1,nt)
 close(111)

33 irec=0
 open(112,file='uv.dat',form='unformatted',access='direct',recl=ny*nx)
    do t=1,nt
       do k=1,nz
          irec=irec+1
          read(112,rec=irec)u(:,:,k,t)
       enddo 
       do k=1,nz
          irec=irec+1
          read(112,rec=irec)v(:,:,k,t)
       enddo
    enddo
 close(112)

! open(56,file='u.dat',action='write',form='formatted') 
!   do t=1,3
!      do k=1,nz
!         write(56,'(48f15.5)')u(:,:,k,t)
!      enddo
!   enddo
 close(56)
deallocate(vort,q)
!========================================================
!=============Vertical Velocity==========================
allocate(w(nx,ny,nz,nt),p(nx,ny,nz,nt))
allocate(pottemp(nx,ny,nz,nt))
    do t=1,nt
       do k=1,nz
          p(:,:,k,t) = 1000*exp(ps(:,:,t))*exp(-z(k)/h0)
          !1000*exp(ps) is the exact ps
       enddo
    enddo
    w(:,:,1,:)=dz(1)*(div(:,:,1,:))
    do k=2,nz
       w(:,:,k,:)=w(:,:,k-1,:) + dz(k)*(div(:,:,k,:))
    enddo
!========================================================
!============Potential Temp==============================
 pottemp(:,:,:,:) = temp(:,:,:,:)*(p0/p(:,:,:,:))**kp
!========================================================
 rho0(:) = rhor*exp(-z(:)/h0)
deallocate(p,ps) 
allocate(ub(nx,ny,nz,nt),vb(nx,ny,nz,nt),wb(nx,ny,nz,nt))
allocate(up(nx,ny,nz,nt),vp(nx,ny,nz,nt),wp(nx,ny,nz,nt))
allocate(pottempb(nx,ny,nz,nt),pottempp(nx,ny,nz,nt))
 call zonalmean(u,ub,nx,ny,nz,nt)
 call zonalmean(v,vb,nx,ny,nz,nt)
 call zonalmean(w,wb,nx,ny,nz,nt)
 call zonalmean(pottemp,pottempb,nx,ny,nz,nt)
 call zonaldev(u,ub,up,nx,ny,nz,nt)
 call zonaldev(v,vb,vp,nx,ny,nz,nt)
 call zonaldev(w,wb,wp,nx,ny,nz,nt)
deallocate(u,v,w)
 call zonaldev(pottemp,pottempb,pottempp,nx,ny,nz,nt)
allocate(upvp(nx,ny,nz,nt),vppottempp(nx,ny,nz,nt),wpup(nx,ny,nz,nt))
allocate(vbpottemppb(nx,ny,nz,nt))
allocate(upvpb(nx,ny,nz,nt),vppottemppb(nx,ny,nz,nt),wpupb(nx,ny,nz,nt))
 upvp=up*vp
 vppottempp=vp*pottempp
 wpup=wp*up
 call zonalmean(upvp,upvpb,nx,ny,nz,nt)
 call zonalmean(vppottempp,vppottemppb,nx,ny,nz,nt)
 call zonalmean(wpup,wpupb,nx,ny,nz,nt)
deallocate(up,vp,wp,upvp,vppottempp,wpup)
allocate(pottempbz(nx,ny,nz,nt),ubcosphiphi(nx,ny,nz,nt))

 call zderiv(ub,ubz,nx,ny,nz,nt,dz,1)
 ! num=1 is ubz(:,:,1,:)= (ub(:,:,1,:)-0)/dz(1)*0.5
 ! asssume surface velocity is zero
 call zderiv(pottempb,pottempbz,nx,ny,nz,nt,dz,2)
 ! num=2 is pottempbz(:,:,1,:)= pottempbz(:,:,2,:)
 ! assume pottemp vert shear at surface 2 vert layers are very close 
 call phideriv(ub,ubcosphiphi,nx,ny,nz,nt,dy,phi)
open(unit=114,file='epfluxdivy1.dat',form='formatted',status='unknown',action='write')
open(unit=115,file='epfluxdivz1.dat',form='formatted',status='unknown',action='write')
open(unit=116,file='epfluxdiv1.dat',form='formatted',status='unknown',action='write')
 
allocate(epfluxy(nx,ny,nz,nt),epfluxz(nx,ny,nz,nt))
 do t=1,nt
    do k=1,nz
       do j=1,ny
        epfluxy(:,j,k,t) = rho0(k)*a*cos(phi(j))*(ubz(:,j,k,t)* &
                    vppottemppb(:,j,k,t)/pottempbz(:,j,k,t) - upvpb(:,j,k,t))
        epfluxz(:,j,k,t) = rho0(k)*a*cos(phi(j))*((f(j) - ubcosphiphi(:,j,k,t)/ &
           (a*cos(phi(j))))*vbpottemppb(:,j,k,t)/pottempbz(:,j,k,t) - wpupb(:,j,k,t))
       enddo
    enddo
 enddo

!=======Residual circulation(Hadley Circulation) Vstar Wstar=====================
open(unit=117,file='vbstar1.dat',form='formatted',status='unknown',action='write')
open(unit=118,file='wbstar1.dat',form='formatted',status='unknown',action='write')
open(unit=119,file='vbstaradv1.dat',form='formatted',status='unknown',action='write')
open(unit=120,file='wbstaradv1.dat',form='formatted',status='unknown',action='write')

allocate(rho0ppb1(nx,ny,nz,nt),rho0ppb2(nx,ny,nz,nt),rho0ppbz(nx,ny,nz,nt))
allocate(rho0ppbphi(nx,ny,nz,nt),vbstar(nx,ny,nz,nt),wbstar(nx,ny,nz,nt))
allocate(vbstaradv(nx,ny,nz,nt),wbstaradv(nx,ny,nz,nt))
 do k=1,nz
    rho0ppb1(:,:,k,:)=rho0(k)*vppottemppb(:,:,k,:)/pottempbz(:,:,k,:)
 enddo
 call zderiv(rho0ppb1,rho0ppbz,nx,ny,nz,nt,z,2)
 do k=1,nz
 vbstar(:,:,k,:) = vb(:,:,k,:) - rho0ppbz(:,:,k,:)/rho0(k)
 enddo

 do j=1,ny
    rho0ppb2(:,j,:,:)=cos(phi(j))*vppottemppb(:,j,:,:)/pottempbz(:,j,:,:)
 enddo
 call phideriv(rho0ppb2,rho0ppbphi,nx,ny,nz,nt,dy,phi)
 do j=1,ny
 wbstar(:,j,:,:) = wb(:,j,:,:) + rho0ppbphi(:,j,:,:)/(a*cos(phi(j)))
 enddo
do j=1,ny
   do k=1,nz
      vbstaradv(:,j,k,:)= vbstar(:,j,k,:)*(ubcosphiphi(:,j,k,:)/ &
                      (a*cos(phi(j)))-f(j))*rho0(k)*a*cos(phi(j))
      wbstaradv(:,j,k,:)= wbstar(:,j,k,:)*ubz(:,j,k,:)*rho0(k)*a*cos(phi(j))
   enddo
enddo
! vbstaradv wbstaradv have same unit as epfluxdiv
!================================================================================

deallocate(rho0,ub,ubz,vppottemppb,upvpb,ubcosphiphi,pottempbz,vbpottemppb,wpupb)
allocate(epfluxdivy(nx,ny,nz,nt),epfluxdivz(nx,ny,nz,nt),epfluxdiv(nx,ny,nz,nt))
 call phideriv(epfluxy,epfluxdivy,nx,ny,nz,nt,dy,phi)
 call zderiv(epfluxz,epfluxdivz,nx,ny,nz,nt,z,2)
 do j=1,ny
 epfluxdiv(:,j,:,:)=epfluxdivy(:,j,:,:)/(a*cos(phi(j))) + epfluxdivz(:,j,:,:)
 enddo
! do t=1,nt
!    do k=1,nz
!       write(114,'(48f15.5)') epfluxdivy(:,:,k,t)
!       write(115,'(48f15.5)') epfluxdivz(:,:,k,t)
!       write(116,'(48f15.5)') epfluxdiv(:,:,k,t)
!       write(117,'(48f15.5)') vbstar(:,:,k,t)
!       write(118,'(48f15.5)') wbstar(:,:,k,t)
!       write(119,'(48f15.5)') vbstaradv(:,:,k,t)
!       write(120,'(48f15.5)') wbstaradv(:,:,k,t)
!    enddo
! enddo
! do t=1,nt
!    do k=1,nz
!       write(114,'(1f15.5)') epfluxdivy(1,:,k,t)
!       write(115,'(1f15.5)') epfluxdivz(1,:,k,t)
!       write(116,'(1f15.5)') epfluxdiv(1,:,k,t)
!       write(117,'(1f15.5)') vbstar(1,:,k,t)
!       write(118,'(1f15.5)') wbstar(1,:,k,t)
!       write(119,'(1f15.5)') vbstaradv(1,:,k,t)
!       write(120,'(1f15.5)') wbstaradv(1,:,k,t)
!    enddo
! enddo
! close(114)
! close(115)
! close(116)
! close(117)
! close(118)
! close(119)
! close(120)

!=========================Monthly Mean=====================================
allocate(epfluxdivyt(1,ny,nz,int(nt/30)),epfluxdivzt(1,ny,nz,int(nt/30)))
allocate(epfluxdivt(1,ny,nz,int(nt/30)),vbstart(1,ny,nz,int(nt/30)))
allocate(wbstart(1,ny,nz,int(nt/30)),vbstaradvt(1,ny,nz,int(nt/30)))
allocate(wbstaradvt(1,ny,nz,int(nt/30)))

call timemean(epfluxdivy,epfluxdivyt,nx,ny,nz,nt,30)
call timemean(epfluxdivz,epfluxdivzt,nx,ny,nz,nt,30)
call timemean(epfluxdiv,epfluxdivt,nx,ny,nz,nt,30)
call timemean(vbstar,vbstart,nx,ny,nz,nt,30)
call timemean(wbstar,wbstart,nx,ny,nz,nt,30)
call timemean(vbstaradv,vbstaradvt,nx,ny,nz,nt,30)
call timemean(wbstaradv,wbstaradvt,nx,ny,nz,nt,30)
  do t=1,int(nt/30)
     do k=1,nz
        write(114,'(1f15.5)') epfluxdivyt(1,:,k,t)
        write(115,'(1f15.5)') epfluxdivzt(1,:,k,t)
        write(116,'(1f15.5)') epfluxdivt(1,:,k,t)
        write(117,'(1f15.5)') vbstart(1,:,k,t)
        write(118,'(1f15.5)') wbstart(1,:,k,t)
        write(119,'(1f15.5)') vbstaradvt(1,:,k,t)
        write(120,'(1f15.5)') wbstaradvt(1,:,k,t)
     enddo
  enddo
  close(114)
  close(115)
  close(116)
  close(117)
  close(118)
  close(119)
  close(120)
!===========================================================================
stop


!==============All SUBROUTINES==============================================
!========TimeAvg,ZonalAvg,ZonalDev,ZDeriv,YDeriv============================

 contains
 subroutine timemean(varin,varout,nx,ny,nz,nt,ndays)
 implicit none
 integer,intent(in)::nx,ny,nz,nt,ndays
 real,intent(in) :: varin(nx,ny,nz,nt)
 real,intent(out):: varout(1,ny,nz,int(nt/ndays))
 integer :: j,k, t, tt, tt1
 tt=1
 do t=1,int(nt/ndays)
    tt1=tt+ndays-1
 do j=1,ny
 do k=1,nz
 !nt/ndays=63 months
    varout(1,j,k,t) = sum(varin(1,j,k,tt:tt1))/(ndays)
 enddo 
 enddo
 enddo
    tt= tt + ndays 
 return
 end subroutine timemean

 subroutine zonalmean(varin,varbout,nx,ny,nz,nt)
 implicit none
 integer,intent(in) :: nx,ny,nz,nt
 real,intent(in) :: varin(nx,ny,nz,nt)
 real,intent(out):: varbout(nx,ny,nz,nt)
 integer  :: i,j,k,t
 do t=1,nt
 do k=1,nz
 do j=1,ny
 varbout(:,j,k,t)=sum(varin(:,j,k,t))/nx
 enddo
 enddo 
 enddo 
 return
 end subroutine zonalmean

 subroutine zonaldev(varin,varbin,varpout,nx,ny,nz,nt)
 implicit none
 integer,intent(in) :: nx,ny,nz,nt
 real,intent(in) :: varin(nx,ny,nz,nt),varbin(nx,ny,nz,nt)
 real,intent(out):: varpout(nx,ny,nz,nt)
 integer  :: i,j,k,t
 varpout(:,:,:,:) = varin(:,:,:,:)-varbin(:,:,:,:)
 return
 end subroutine zonaldev

 subroutine zderiv(varin,varout,nx,ny,nz,nt,dz,num)
 implicit none
 integer,intent(in) :: nx,ny,nz,nt
 real,intent(in) :: dz(nz)
 real,intent(in) :: varin(nx,ny,nz,nt)
 real,intent(out):: varout(nx,ny,nz,nt)
 integer :: k,num
 do k=2,nz
   varout(:,:,k,:)= (varin(:,:,k,:) - varin(:,:,k-1,:))/(0.5*dz(k)+0.5*dz(k-1))
 enddo
 if(num .eq. 1)then
 varout(:,:,1,:)= (varin(:,:,1,:)-0)/dz(1)*0.5
 else
 varout(:,:,1,:)= varout(:,:,2,:)
 endif
 return 
 end subroutine zderiv

 subroutine phideriv(varin,varout,nx,ny,nz,nt,dy,phi)
 implicit none
 integer,intent(in) :: nx,ny,nz,nt
 real,intent(in) :: phi(ny),dy
 real,intent(in) :: varin(nx,ny,nz,nt)
 real,intent(out):: varout(nx,ny,nz,nt)
 integer :: j
 do j=2,ny
   varout(:,j,:,:)= (varin(:,j,:,:)*cos(phi(j)) - varin(:,j-1,:,:)*cos(phi(j-1)))/dy
 enddo
 varout(:,1,:,:)=varout(:,ny,:,:)
 return
 end subroutine phideriv


 end program epflux
