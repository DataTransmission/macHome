  program tcdetect 
  !use mini
  use netcdf !use the library when makefile
  use myncread
  
  implicit none
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
  !filenameame='case7.cam2.h1.0001-01-02-00000.nc'

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
  integer,allocatable,dimension(:):: tcrow, tccol, tclev
  real,allocatable,dimension(:)::    dx
  real,allocatable,dimension(:)::    lat, lon, lev
  real,allocatable,dimension(:,:)::  sfspd, psl, vortthrmat, sfspdthrmat, tzintathrmat
  real,allocatable,dimension(:,:)::  vx, uy, vort850, tzint, tzintavg, tzinta, matrix 
  real,allocatable,dimension(:,:,:)::u, v, us, vs, spd, t, ta, tavg, taavg, trackmap, p
  real,allocatable,dimension(:)         :: zvmat !matrix to store tv to judge minval(zv)
  real:: signs, dy
  integer:: nmonth(12)
  integer:: ii,jj,tt, status, ncID
  integer:: i,j,k,test(2),pp(2)
  integer:: nyr,nmon,nday
  character(len=50):: filename, filenamebsf
!  character:: string*80
  character(len=2)::x1,x2,x3,yrten,mnten,dyten
  character(len=4) :: fmt ! format descriptor
  real,parameter:: pref= 1000 ! for calculating pressure varying on fixed sigma coord
  fmt = '(I1)' ! an integer of width 5 with zeros at the left
   status= nf90_open('TCParam.nc',nf90_NoWrite, ncID)
   call readnc_int('nx'      ,nx,ncID)
   call readnc_int('ny'      ,ny,ncID)
   call readnc_int('nz'      ,nz,ncID)
   call readnc_int1d('tclev'   ,tclev,ncID) !850 700 500 300
   call readnc_int('nx1'     ,nx1,ncID)
   call readnc_int('nx2'     ,nx2,ncID)
   call readnc_int('ny1'     ,ny1,ncID)
   call readnc_int('ny2'     ,ny2,ncID)
   call readnc_real('dy'      ,dy,ncID)
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
   call readnc_int1d('tcrow',tcrow,ncID)
   call readnc_int1d('tccol',tccol,ncID)
   status=nf90_close(ncID)
   status= nf90_open('TCThrMat.nc',nf90_NoWrite, ncID)
   call readnc_real2d('vortthrmat'    ,vortthrmat,ncID,nx,ny)
   call readnc_real2d('sfspdthrmat'   ,sfspdthrmat,ncID,nx,ny)
   call readnc_real2d('tzintathrmat'  ,tzintathrmat,ncID,nx,ny)
   status=nf90_close(ncID)
   nmonth=(/1,2,3,4,5,6,7,8,9,10,11,12/)
   !day(1)='31';day(2)='31';day(3)='30';day(4)='30';day(5)='31';day(6)='28'
   allocate(vort850(nx,ny),u(nx,ny,size(tclev)),v(nx,ny,size(tclev)),us(nx,ny,nz),vs(nx,ny,nz))
   allocate(sfspd(nx,ny),psl(nx,ny),t(nx,ny,size(tclev)))
   allocate(spd(nx,ny,size(tclev)))
   allocate(tavg(nx,ny,size(tclev)),ta(nx,ny,size(tclev)),taavg(nx,ny,size(tclev)))
   allocate(matrix(5,800)) ! 10 Means we don't get more then 10 TC a day
   allocate(vx(nx,ny),uy(nx,ny))
   tt=0
   ii=0
   do nyr= 9,9
     do nmon= 1,size(nmonth)
       do nday= 1,day(nmon)
         yrten=''; mnten=''; dyten=''
         if (nyr < 10) yrten='0'; if (nmon < 10) mnten='0'; if (nday < 10) dyten='0'
         write(x1,'(I2)')nyr; write(x2,'(I2)')nmon; write(x3,'(I2)')nday! converting integer to string using a 'internal file'
         filenamebsf='bsf.case7.cam2.h1.00'//trim((yrten))//trim(adjustl(x1))//'-'//trim((mnten))//trim(adjustl(x2))//'-'//trim((dyten))//trim(adjustl(x3))//'-00000.nc'
         print*,filenamebsf
         status= nf90_open(filenamebsf,nf90_NoWrite, ncID)
         if (status==2) goto 11  ! if no files exist(status==2) then goto 11
         call readnc_real3d('U'       ,u,ncID,nx,ny,size(tclev))
         call readnc_real3d('V'       ,v,ncID,nx,ny,size(tclev))
         call readnc_real3d('T'       ,t,ncID,nx,ny,size(tclev))
         status= nf90_close(ncID)
         filename='case7.cam2.h1.00'//trim((yrten))//trim(adjustl(x1))//'-'//trim((mnten))//trim(adjustl(x2))//'-'//trim((dyten))//trim(adjustl(x3))//'-00000.nc'
         status= nf90_open(filename,nf90_NoWrite, ncID)
         if (status==2) goto 11  ! if no files exist(status==2) then goto 11
         status= nf90_open(filename,nf90_NoWrite, ncID)
         call readnc_real2d('PSL'     ,psl,ncID,nx,ny)
         call readnc_real3d('U'     ,us,ncID,nx,ny,nz)
         call readnc_real3d('V'     ,vs,ncID,nx,ny,nz)
         status= nf90_close(ncID)
         sfspd = sqrt(us(:,:,nz)**2.0 +vs(:,:,nz)**2.0)
         !find 850mb that's on each sigma lev k
         k=1
            vx(1,2:ny-1) = (v(2,2:ny-1,k) - v(1,2:ny-1,k))/dx(2:ny-1)
            do i= 2,nx-1
               vx(i,2:ny-1) = (v(i+1,2:ny-1,k) - v(i-1,2:ny-1,k))/(2.0*dx(2:ny-1))
            enddo
            vx(nx,2:ny-1) = (v(nx,2:ny-1,k) - v(nx-1,2:ny-1,k))/dx(2:ny-1)
            uy(1:nx,1) = (u(1:nx,2,k) - u(1:nx,1,k))/dy
            do j = 2,ny-1
              uy(1:nx,j) = (u(1:nx,j+1,k) - u(1:nx,j-1,k))/(2.0*dy)
            enddo
            uy(1:nx,ny) = (u(1:nx,ny,k) - u(1:nx,ny-1,k))/dy
         vort850(:,2:ny-1) = vx(:,2:ny-1) - uy(:,2:ny-1)
     !select the center region (7*7) of tzint for aave
      !v is speed
         !data tzintavg / ny*nx * 0.0 /!, tzinta/ ny*nx * 0.0/
         !tzint= sum(t(tclev(1:3),:,:),1); !300,500,700
         do i=si(1)-hfbox,atl2(size(atl2))+hfbox
           do j=fortyS-hfbox,fortyN+hfbox
             do k=1,size(tclev)
               spd(i,j,k)= sqrt(u(i,j,k)**2 + v(i,j,k)**2)
               tavg(i,j,k)= sum(sum(t(i-hfbox:i+hfbox,j-hfbox:j+hfbox,k),1))/boxarea
               ta(i,j,k)= t(i,j,k)- tavg(i,j,k)
               if (j>=fortyS .AND. j<=fortyN .AND. i>=si(1) .AND. i<= atl2(size(atl2)))  taavg(i,j,k)= sum(sum(ta(i-hfbox:i+hfbox,j-hfbox:j+hfbox,k),1))/boxarea       
             enddo
           enddo
         enddo
         do i=1,size(tcrow)
             pp= minloc(psl(tccol(i)-hfbox:tccol(i)+hfbox,tcrow(i)-hfbox:tcrow(i)+hfbox)) 
             signs=1.0  !the default is NH
!             print*,tzintathrmat(tccol(i),tcrow(i)) 
             if (tcrow(i)<=zeroS) signs=-1.0 !this is to find negative vorticity for TC in SH
             if (signs*vort850(tccol(i),tcrow(i)) .GT. vortthrmat(tccol(i),tcrow(i)) .AND. &
                 maxval(sfspd(tccol(i)-hfbox:tccol(i)+hfbox,tcrow(i)-hfbox:tcrow(i)+hfbox)) .GT. sfspdthrmat(tccol(i),tcrow(i)) .AND. &
                 pp(1)==hfbox+1 .AND. pp(2)==hfbox+1 .AND. &
!                 abs(sum(taavg(tccol(i),tcrow(i),1:3)))*2 .GT. tzintathrmat(tccol(i),tcrow(i)) .AND. &
                 taavg(tccol(i),tcrow(i),2) .GT. 0 .AND. & !300,500,700 > 0
                 taavg(tccol(i),tcrow(i),3) .GT. 0 .AND. & !300,500,700 > 0
                 taavg(tccol(i),tcrow(i),4) .GT. 0 .AND. &! .AND. & !300,500,700 > 0
                 taavg(tccol(i),tcrow(i),4) .GT. taavg(tccol(i),tcrow(i),1) .AND. & !.AND. & !300> 850 change to 500>850
                 sum(sum(spd(tccol(i)-hfbox:tccol(i)+hfbox,tcrow(i)-hfbox:tcrow(i)+hfbox,1),1))/boxarea .GT.  &
                 sum(sum(spd(tccol(i)-hfbox:tccol(i)+hfbox,tcrow(i)-hfbox:tcrow(i)+hfbox,4),1))/boxarea & !too hard to satisfy 
                 )then
               ii=ii+ 1
print*,ii
               matrix(1,ii)= nyr
               matrix(2,ii)= nmonth(nmon)
               matrix(3,ii)= nday
               matrix(4,ii)= tcrow(i)
               matrix(5,ii)= tccol(i)
             endif    
         enddo
       enddo
11     continue
     enddo
    write(x1,'(I2)')nyr 
    open(unit=153,file='tcmp'//trim(adjustl(x1))//'.dat',form='formatted',status='unknown',action='write')
    write(153,'(5f15.1)') matrix(:,1:ii)
    close(153)
   enddo
   print*,"Done Detecting"
   stop

end program
