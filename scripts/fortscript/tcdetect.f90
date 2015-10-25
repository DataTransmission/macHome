  program tcdetect 
  !use mini
  use netcdf !use the library when makefile
  use myncread
  
  implicit none
! Contact gchen@rsmas.miami.edu 
! thresholds:
!  pick summer season JJA for NH, DJF for SH, merge daily data into a singal matrix (Basinwise)
!  1) 850h-Pa vort > vort threshold (neg in SH)
!  2) 7x7 box max sfcspd > spd threshold 
!  3) min psl in 7x7 center
!  4) 7x7 tempa avg and three lev avg (300, 500, 700) > tempa threshold
!  5) 7x7 tempa avg over 300, 500, 700 > 0
!  6) 7x7 tempa avg over 300 > 850
!  7) 7x7 spd avg over 850 > 300
!  8) connect center of storm over time step (6hrly need 5.6 deg, daily need 8.5 deg)
!  9) defined as storm if criteria last > 2 day (1.5 day for 6hrly output) 
  !filename='case7.cam2.h1.0001-01-02-00000.nc'

! DESIGN experiments to test TC heat transport, 
!  1) detect disturbance and put in the TC 
!      PROBLEM: after bogus, don't loop over bogus
!      SOLUTION1: apply a moving circular coordinate tracer for this cyclone, 
!                 anything in this domain are not taken into detection domain.
!      SOLUTION2: see the TC as same one if the next time step distance (D) btw previous time step D= X(t1)- X(tn0) 
!                 is within maximal distance (MD) (use historical fastest moving TC)
!                 don't bogus if it's within this distance, recalc the distance (D) for every future D(tn+1)= X(tn+1)- X(tn)
!                 every time D(tn+1) < MD then don't bogus on that point
!  2) Bogus in only at the date of historical TC
!      PROBLEM: background not suitable for genesis
!      SOLUTION: stochastic generate TC, match with historical TC only in frequency
!  3) Build an hourly Markov model for TC
  integer:: nx, ny, nx1, nx2, ny1, ny2, nz
  integer:: zeroN, zeroS, fortyN, fortyS, twenty5N
  integer:: hfbox, boxarea, days
  integer,allocatable,dimension(:):: ni, wnp, enp, enp1, enp2, day
  integer,allocatable,dimension(:):: atl, atl1, atl2, si, aus, sp
  integer,allocatable,dimension(:):: tcRowInd, tcColInd, tclev
  real,allocatable,dimension(:)::    dx, dy
  real,allocatable,dimension(:)::    lat, lon, lev
  real,allocatable,dimension(:,:)::  sfspd, psl, vortthrmat, sfspdthrmat, tzintathrmat
  real,allocatable,dimension(:,:)::  vx, uy, vort850, tzint, tzintavg, tzinta, matrix, trackmap
  real,allocatable,dimension(:,:,:)::u, v, spd, t, ta, tavg, taavg
  real,allocatable,dimension(:)         :: zvmat !matrix to store tv to judge minval(zv)
  integer:: ii,jj, status, ncID
  integer:: i,j,k,test(2),pp(2)
  integer:: nyr,nmon,nday
  character(len=50):: filename
!  character:: string*80
  character(len=1)::x1,x2,x3
  character(len=4) :: fmt ! format descriptor
  fmt = '(I1)' ! an integer of width 5 with zeros at the left
  nyr=1; nmon=1; nday=1
   
   status= nf90_open('TCParam.nc',nf90_NoWrite, ncId)
   call readnc_int('nx'      ,nx,ncID)
   call readnc_int('ny'      ,ny,ncID)
   call readnc_int('nz'      ,nz,ncID)
   call readnc_int1d('tclev'   ,tclev,ncID)
   call readnc_int('nx1'     ,nx1,ncID)
   call readnc_int('nx2'     ,nx2,ncID)
   call readnc_int('ny1'     ,ny1,ncID)
   call readnc_int('ny2'     ,ny2,ncID)
   call readnc_real1d('dy'      ,dy,ncID)
   call readnc_real1d('dx'      ,dx,ncID)
   call readnc_int1d('day'     ,day,ncID)
   call readnc_int('days'    ,days,ncID)
   call readnc_int('zeroN'   ,zeroN,ncID)
   call readnc_int('zeroS'   ,zeroS,ncID)
   call readnc_int('fortyN'  ,fortyN,ncID)
   call readnc_int('fortyS'  ,fortyS,ncID)
   call readnc_int('twenty5N',twenty5N,ncID)
   call readnc_int('hfbox'   ,hfbox,ncID)
   call readnc_int('boxarea' ,boxarea,ncID)
   call readnc_int1d('ni'      ,ni,ncID)
   call readnc_int1d('wnp'     ,wnp,ncID)
   call readnc_int1d('enp'     ,enp,ncID)
   call readnc_int1d('enp1'    ,enp1,ncID)
   call readnc_int1d('enp2'    ,enp2,ncID)
   call readnc_int1d('atl'     ,atl,ncID)
   call readnc_int1d('atl1'    ,atl1,ncID)
   call readnc_int1d('atl2'    ,atl2,ncID)
   call readnc_int1d('si'      ,si,ncID)
   call readnc_int1d('aus'     ,aus,ncID)
   call readnc_int1d('sp'      ,sp,ncID)
   call readnc_int1d('tcRowInd',tcRowInd,ncID)
   call readnc_int1d('tcColInd',tcColInd,ncID)
   status= nf90_open('TCThrMat.nc',nf90_NoWrite, ncID)
   call readnc_real2d('vortthrmat'    ,vortthrmat,ncID,nx,ny)
   call readnc_real2d('sfspdthrmat'   ,sfspdthrmat,ncID,nx,ny)
   call readnc_real2d('tzintathrmat'  ,tzintathrmat,ncID,nx,ny)
   do nyr=nyr1,nyr2
     do nmon=nmon1,nmon2
       do nday=nday1,nday2
   write (x1,fmt) nyr ! converting integer to string using a 'internal file'
   write (x2,fmt) nmon
   write (x3,fmt) nday
   if (nyr < 10) then
     if (nday < 10) then
       filename=(['case7.cam2.h1.000'//x1//'-'//x2//'-0'//x3//'-00000.nc']);
     else  
       filename=(['case7.cam2.h1.000'//x1//'-'//x2//'-'//x3// '-00000.nc']);
     endif
   else
     if (nday < 10) then
       filename=(['case7.cam2.h1.00'//x1//'-'//x2//'-0'//x3//'-00000.nc']);
     else 
       filename=(['case7.cam2.h1.00'//x1//'-'//x2//'-'//x3// '-00000.nc']);
     endif
   endif
   filename = 'case7.cam2.h1.000'//x1//'-0'//x2//'-0'//x3//'-00000.nc'
   print*,filename
   allocate(vort850(nx,ny),u(nx,ny,nz),v(nx,ny,nz))
   allocate(sfspd(nx,ny),psl(nx,ny),t(nx,ny,nz))
   status= nf90_open(filename,nf90_NoWrite, ncId)
!   call readnc_real1d('lat'     ,lat,ncID     )
!   call readnc_real1d('lon'     ,lon,ncID     )
!   call readnc_real1d('lev'     ,lev,ncID     )
   call readnc_real2d('PSL'     ,psl,ncID,nx,ny)
   call readnc_real3d('U'       ,u,ncID,nx,ny,nz)
   call readnc_real3d('V'       ,v,ncID,nx,ny,nz)
   call readnc_real3d('T'       ,t,ncID,nx,ny,nz)
   sfspd = sqrt(u(:,:,nz)**2.0 +v(:,:,nz)**2.0)
   allocate(vx(nx,ny),uy(nx,ny))
   k= tclev(4) ! at 850
   vx(1,:) = (v(2,:,k) - v(1,:,k))/dx
   do i= 2,nx-1
     vx(i,:) = (v(i+1,:,k) - v(i-1,:,k))/(2.0*dx)
   enddo
   vx(nx,:) = (v(nx,:,k) - v(nx-1,:,k))/dx
   uy(:,1) = (u(:,2,k) - u(:,1,k))/dy
   do i = 2,ny-1
     uy(:,i) = (u(:,i+1,k) - u(:,i-1,k))/(2.0*dy)
   enddo
   uy(:,ny) = (u(:,ny,k) - u(:,ny-1,k))/dy
   
   vort850 = vx - uy
!select the center region (7*7) of tzint for aave
!v is speed
   allocate(spd(nx,ny,size(tclev)))
   do k=1,size(tclev)
     spd(:,:,k)=sqrt(u(:,:,tclev(k))**2 + v(:,:,tclev(k))**2)
   enddo
   !data tzintavg / ny*nx * 0.0 /!, tzinta/ ny*nx * 0.0/
   allocate(tavg(nx,ny,size(tclev)),ta(nx,ny,size(tclev)),taavg(nx,ny,size(tclev)))
   allocate(matrix(10,2) ! 10 Means we don't get more then 10 TC a day
   !tzint= sum(t(tclev(1:3),:,:),1); !300,500,700
   do j=fortyS-hfbox,fortyN+hfbox
     do i=si(1)-hfbox,atl2(-1)+hfbox
       do k=1,size(tclev)
         tavg(i,j,k)= sum(sum(t(i-hfbox:i+hfbox,j-hfbox:j+hfbox,tclev(k)),1),1)/boxarea
         ta(i,j,k)= t(i,j,tclev(k))- tavg(i,j,k)
       enddo
     enddo
   enddo
   do j=fortyS,fortyN
     do i=si(1),atl2(-1)
       do k=1,size(tclev)
         taavg(i,j,k)= sum(sum(ta(i-hfbox:i+hfbox,j-hfbox:j+hfbox,k),1))/boxarea       
       enddo
     enddo
   enddo
   ii=0
   do i=1,size(tcRowInd)
           pp= minloc(psl(tcColInd(i)-hfbox:tcColInd(i)+hfbox,tcRowInd(i)-hfbox:tcRowInd(i)+hfbox)) 
       if (vort850(tcColInd(i),tcRowInd(i)) .GT. vortthrmat(tcColInd(i),tcRowInd(i)) .AND. &
           maxval(sfspd(tcColInd(i)-hfbox:tcColInd(i)+hfbox,tcRowInd(i)-hfbox:tcRowInd(i)+hfbox))> sfspdthrmat(tcColInd(i),tcRowInd(i)) .AND. &
           pp(1)==hfbox+1 .AND. pp(2)==hfbox+1 .AND. &
!           abs(sum(taavg(tcColInd(i),tcRowInd(i),1:3))/3.0) .GT. tzintathrmat(tcColInd(i),tcRowInd(i)) .AND. &
!           taavg(tcColInd(i),tcRowInd(j),1)> 0 .AND. & !300,500,700 > 0
!           taavg(tcColInd(i),tcRowInd(j),2)> 0 .AND. & !300,500,700 > 0
!           taavg(tcColInd(i),tcRowInd(j),3)> 0&! .AND. & !300,500,700 > 0
!           taavg(tcColInd(i),tcRowInd(j),2)> taavg(tcColInd(i),tcRowInd(j),4) & !.AND. & !300> 850 change to 500>850
           sum(sum(spd(tcColInd(i)-hfbox:tcColInd(i)+hfbox,tcRowInd(j)-hfbox:tcRowInd(j)+hfbox,4),1))/boxarea> &
           sum(sum(spd(tcColInd(i)-hfbox:tcColInd(i)+hfbox,tcRowInd(j)-hfbox:tcRowInd(j)+hfbox,1),1))/boxarea & !too hard to satisfy 
           )then
         ii=ii+ 1
print*,vort850(tcColInd(i),tcRowInd(i))
print*,tcRowInd(i),tcColInd(i)
         matrix(ii,1)=tcRowInd(i)
         matrix(ii,2)=tcColInd(i)
       endif    
   enddo
!write(*,'(4F13.8)'),matrix(:,1)
open(unit=111,file='trackmap.dat',form='formatted',status='unknown',action='write')
write(111,'(2f15.8)') matrix(:,1),matrix(:,2)
close(111)
stop

end program
