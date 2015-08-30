           program vortexinit 
           !module vortexinit
           !use kinds
           !use diff
           use netcdf !use the library when makefile
           use myncread
           implicit none
           
           real*8 :: a, omg, rd, g, zt, qt, zq1, zq2
           real*8 :: gam, dp, rp, zp, phic, lamc, fc, e, eps, rad
           real*8 :: d1, d2, d
           real*8 :: r, ps, dpdz, f, tcradius, ztc1 !tangential velocity
           real*8 :: dlam, dphi
           real*8, allocatable, dimension(:,:,:) :: xxx, temp !logp0 
           real*8, allocatable, dimension(:,:)  :: ts, dtdp, psl
           real*8, allocatable, dimension(:,:)  :: p0, q0, t0, tv0, tvt, pt
           real*8, allocatable, dimension(:,:)  :: trckmp 
           real*8, allocatable, dimension(:)    :: lon, lat, lam, phi
           real*8, allocatable, dimension(:)    :: hyam, hybm
           real*8, allocatable, dimension(:,:,:):: ztc,ztcini, plev, qtc, ttc, utc, vtc, pstc, vttc 
           real*8, allocatable, dimension(:,:,:):: divtc, vortc, ux, uy, vx, vy 
           real*8, allocatable, dimension(:,:,:):: q0b, pmod, u, v,  q
           real*8, allocatable, dimension(:,:,:):: u_tcincl, v_tcincl, t_tcincl
           real*8, allocatable, dimension(:,:)  :: psl_tcincl
           !real(r8),allocatable :: z(:,:,:), d(:,:,:), t(:,:,:), w(:,:,:)
           !real(r8),allocatable :: ztcini(:,:,:), plevini(:,:,:)
           integer :: nx, ny, nz, rowsize, colsize
           integer :: i, j, k, ii, jj, itc, inx, nyr, nmon, nday
           integer :: tcxgrid, tcygrid
           integer :: xc, yc, cond, ncond, irec, status, ncID
           !real,parameter, dimension(nz) :: pmod = (/ 900.0, 700.0, 500.0, 300.0, 100.0 /)
           real*8, parameter :: pi= 3.141592653589793, pref= 1000.0
           integer,allocatable,dimension(:):: tclev, day
           character(len=2)::x1,x2,x3,yrten,mnten,dyten      
           character*20 :: xsize 
           character(len=50):: filename
           status= nf90_open('TCParam.nc',nf90_NoWrite, ncID)
           call readnc_int('nx'      ,nx,ncID)
           call readnc_int('ny'      ,ny,ncID)
           call readnc_int('nz'      ,nz,ncID)
           call readnc_int1d('day'     ,day,ncID)
           filename= 'case7.cam2.h1.0009-10-26-00000.nc'
           status= nf90_open(filename,nf90_NoWrite, ncID)
           if (status==2)  print*,'error' 
           call readnc_real1d('hyam',hyam,ncID)
           call readnc_real1d('hybm',hybm,ncID)
           status= nf90_close(ncID)
           !status= nf90_open('TCtrckmp.nc',nf90_NoWrite,ncID)
           !call readnc_int('rowsize',rowsize,ncID)
           !call readnc_int('colsize',colsize,ncID)
!print*,rowsize,colsize
           !call readnc_real2d('trckmp2',trckmp,ncID,colsize,rowsize)
           !status= nf90_close(ncID)
           allocate(trckmp(127,8)) !try to figure out what happened to matlab when saving trckmp2 into nc file
           open(54,file='TCtrckmp.dat')
           read(54,*)trckmp(:,:)
           close(54)
            itc=124 ! this is trackmum 91
!           do nyr= 8,8
!              do nmon= 1,12
!                 do nday= 1,day(nmon)
!                    if (trckmp(1,ii)==nyr .AND. trckmp(2,ii)==nmon .AND. trckmp(3,ii)==nday) then
             nyr=trckmp(itc,1)
             nmon=trckmp(itc,2)
             nday=trckmp(itc,3)
             yrten=''; mnten=''; dyten=''
             if (nyr < 10) yrten='0'; if (nmon < 10) mnten='0'; if (nday < 10) dyten='0'
             write(x1,'(I2)')nyr; write(x2,'(I2)')nmon; write(x3,'(I2)')nday! converting integer to string using a 'internal file'
             filename='case7.cam2.h1.00'//trim((yrten))//trim(adjustl(x1))//'-'//trim((mnten))//trim(adjustl(x2))//'-'//trim((dyten))//trim(adjustl(x3))//'-00000.nc'
             status= nf90_open(filename,nf90_Write, ncID)
             call readnc_real2d('PSL'     ,psl,ncID,nx,ny)
             psl= psl*1D-2
print*,filename
             allocate(pmod(nx,ny,nz))
             do i=1,nx
                do j=1,ny
                   do k=1,nz
                   pmod(i,j,k)= hyam(k)*pref+ hybm(k)* psl(i,j)
                   enddo
                enddo
             enddo
             call readnc_real3d('U'       ,u,ncID,nx,ny,nz)
             call readnc_real3d('V'       ,v,ncID,nx,ny,nz)
             call readnc_real3d('T'       ,temp,ncID,nx,ny,nz)
             call readnc_real2d('TS'      ,ts,ncID,nx,ny)
             call readnc_real3d('Q'       ,q0b,ncID,nx,ny,nz)
             call readnc_real1d('lon'       ,lon,ncID)
             call readnc_real1d('lat'       ,lat,ncID)
             !close(ncID)
!                   endif
!                   ii=ii+1
!                 enddo
!              enddo
!           enddo
           !In a p-coor model, we use the algorithm to calc the background surface pressure
           !distribution for an initial TC. Then by using this ps(r), we could get the z distribution
           !at each pmod(k) level. Since the analytical conversion of p to z is inaccurate in the TC,
           !(due to warm core and pressure perturbation), so need to use Newton's method to make 
           !the plev(k) converge to pmod(k).
                
           !==========================================================================
           !============PARAMETERS====================================================
           !==========================================================================
           rad = dble(pi/180)
           tcxgrid=3                                !the x-dir/y-dir +-grids you have in your domain 
           tcygrid=3
           phic=trckmp(itc,7)*rad                           !phic(lamc)=center lat(lon) of initial vortex, 
           lamc=trckmp(itc,8)*rad
           !q0b=0.025                               !q0=sfc spec hum, 
           !ts=302.15                !background sfc temp
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
           dp=11.15                       !dp=p0 - ps, sfc p diff btwn bckgrd sfc p0 and ps at ini vortex center, 
           rp=282D3                    !rp & zp=ctc for p fit, 
           zp=7D3
           fc=2.0*omg*sin(phic)            !fc=coriolis param, 
           ncond=2
           tcradius = 1D8
           write(xsize,'(i10)') tcxgrid*2+1
           irec=0
           !open(111,file='t_lastday.dat',form='unformatted',access='direct',recl=nx*ny)
           !ii=1
           !allocate(t(nx,ny,nz),logp0(nx,ny,1)) !z(nx,ny,nz),d(nx,ny,nz),t(nx,ny,nz),w(nx,ny,nz),logp0(nx,ny,1))
           !!!!!!!!!!!!INDEX SHOULD BE CONSISTENT WITH OTHER DIMS OF VAR
           !call readdat(xxx   ,nx,ny,nz,irec,0) !could be calc from grads and output binary
           !call readdat(xxx   ,nx,ny,nz,irec,0) !0 = noread; 1 = read
           !call readdat(t,nx,ny,nz,irec,1)
           !call readdat(logp0,nx,ny,1 ,irec,1) !this is called q in Ben's model which is logp0 = sqrt(2)*log(p0)
           !call readdat(q0b,nx,ny,1,irec)
           !call readdat(ts,nx,ny,1,irec)
           !close(111)
           !open(unit=99,file='tb.dat',form='formatted',status='unknown',action='write')
           !do k=1,nz
           !   write(99,'(48f15.5)') t(:,:,k) !u
           !enddo
           !=====suppose to be from model=======
           !do i=1,nx
           !   do j=1,ny
                 !dtdp(i,j)=(log(t(i,j,2)) - log(t(i,j,1)))/(log(pmod(2))-log(pmod(1))) !
                 !ts (i,j) = exp(/log(t(i,j,1)) - dtdp(i,j)*(log(pmod(1)) - (log(1000.0) + logp0(i,j,1)/2**.5)))
                                                !log(1000) + logp0(i,j,1)/2**.t = model ln(ps)
           !   enddo
           !enddo
           !==========================================================================
           !==========================Genesis Zone====================================
           !==========================================================================
           !allocate( lam1(nx), phi1(ny) ) 
           allocate( lam(2*tcxgrid+1), phi(2*tcygrid+1) )
           allocate( p0(size(lam),size(phi)),q0(size(lam),size(phi)),t0(size(lam),size(phi)) )
           dlam = abs(lon(1)-lon(2))
           dphi = abs(lat(1)-lat(2))
           xc= trckmp(itc,5)
           yc= trckmp(itc,4)
           !dlam = 360.0/dble(nx-1)
           !dphi = 180.0/dble(ny-1)
           !do i = 1,nx                              !i grids for hi res
           !   lam1(i) = (i-1)*dlam
           !   if (abs(lam1(i) - lamc) < dlam) then
           !     xc = i
           !   endif
           !enddo
           !do i = 1,ny 
           !   phi1(i) = -90.0 + (i-1)*dphi
           !   if ( abs(phi1(i) - phic) < dphi) then
           !      yc = i
           !   endif
           !enddo
           ii=1
           do i = xc-tcxgrid,xc+tcxgrid    ! 142 143 144 1 2 3 4   if center is at 1
              if (i <= 0) then             !  -2  -1   0 1 2 3 4   <== i will be
                 inx=nx ! go to the other side of 360
              else
                 inx=0
              endif
              p0(ii,1:2*tcygrid+1)=  psl(i+inx,yc-tcygrid:yc+tcygrid)
              q0(ii,1:2*tcygrid+1)=  q0b(i+inx,yc-tcygrid:yc+tcygrid,nz)
              t0(ii,1:2*tcygrid+1)=  ts(i+inx,yc-tcygrid:yc+tcygrid)
              lam(ii) = lon(i+inx)*rad  ! which radiance longitude the TC is at, save it in lam 
              ii = ii+1
           enddo
           phi = lat(yc-tcygrid:yc+tcygrid)*rad
print*,phic,lamc!,lam!,phi/rad,lam/rad
           !phic=phic*rad                        !phic(lamc)=center lat(lon) of initial vortex, 
           !lamc=lamc*rad
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
                    ztcini(size(lam),size(phi),nz))
           !==========================================================================
           !================CALCULATE VAR=============================================
           !==========================================================================
!print*,pmod(1,1,1)
!print*,nz
!stop
           do i=1,size(lam) !1:15
                 do j=1,size(phi) !1:8
                    r  = a*acos(sin(phic)*sin(phi(j)) + cos(phic)*cos(phi(j))*cos(lam(i)-lamc))
                    d1 = sin(phic)*cos(phi(j)) - cos(phic)*sin(phi(j))*cos(lam(i) - lamc)
                    d2 = cos(phic)*sin(lam(i) - lamc)                                             !opposite sign in SH
                    d  = max(e, (d1**(2.0) + d2**(2.0))**(0.5))
                    ps = p0(i,j) - dp*exp(-(r/rp)**(1.5))  !ps is the sfc-p to initial TC, p0 is the surface p from the model output
                  do k=1,nz
                       if (pmod(i,j,k) .GE. pt(i,j) .AND. ps .GE. pmod(i,j,k)) then 
                           cond = 1
                           ztc(i,j,k) = tv0(i,j)/gam*(1.0 - (pmod(i,j,k)/ps)**(rd*gam/g)) !ztc at each pmod level 
                           ztcini(i,j,k)=ztc(i,j,k);
             !              plevini(i,j,k) = (p0(i,j) - dp*exp(-(r/rp)**1.5)*exp(-(ztc(i,j,k)/zp)**2))*((tv0 -  &
             !                      gam*ztc(i,j,k))/tv0)**(g/(rd*gam)); ! pressure coor + perturbation 
!print*,tcradius
                           if ( r .LE. tcradius ) then
                              do while (cond .LT. ncond)
                                 plev(i,j,k) = (p0(i,j) - dp*exp(-(r/rp)**(1.5))*exp(-(ztc(i,j,k)/zp)**(2.0)))*((tv0(i,j) -  &
                                        gam*ztc(i,j,k))/tv0(i,j))**(g/(rd*gam)) ! pressure coor + perturbation 
                                 f = pmod(i,j,k) - plev(i,j,k)
                                 dpdz = 2.0*dp*ztc(i,j,k)/zp**(2.0)*exp(-(r/rp)**(1.5))*exp(-(ztc(i,j,k)/zp)**(2.0))* &
                                       ((tv0(i,j) - gam*ztc(i,j,k))/tv0(i,j))**(g/(rd*gam)) - g/(rd*tv0(i,j))* &
                                       (p0(i,j) - dp*exp(-(r/rp)**(1.5))*exp(-(ztc(i,j,k)/zp)**(2.0)))* &
                                       ((tv0(i,j) - gam*ztc(i,j,k))/tv0(i,j))**(g/(rd*gam) - 1.0)
                                 ztc1 = ztc(i,j,k)
!print*,ztc1
                                 ztc(i,j,k)  = ztc(i,j,k) + f/dpdz
                                 if (ztc(i,j,k) - ztc1 .LE. eps) then 
                                   ztc(i,j,k) = ztc1
                                   cond = ncond
                                 endif 
                              enddo
                           endif
                           qtc(i,j,k) = q0(i,j)*exp((-ztc(i,j,k)/zq1))*exp(-(ztc(i,j,k)/zq2)**(2.0)) !q(z) dist for vortex
                                               !Refer to http://www.applet-magic.com/gradientwind.htm for Gradient wind Balance dep in Hemisphere
                           vttc(i,j,k) = -abs(fc*r*0.5) + ( 0.25*fc*fc*r*r -  &            !THIS WILL CREATE PROBLEM IN SH 
!print*,r,rp,rd,tv0(i,j),gam,ztc(i,j,k),zp,p0(i,j),dp
!print*,tv0(i,j)-gam*ztc(i,j,k)
!print*,(2*rd*(tv0(i,j)-gam*ztc(i,j,k))*ztc(i,j,k))/(g*zp**2)
                                 ( 1.5*(r/rp)**(1.5)*rd*(tv0(i,j)-gam*ztc(i,j,k)) ) / &
                               ( 1.0 + 2.0*rd*(tv0(i,j) - gam*ztc(i,j,k)) * ztc(i,j,k) / & 
                                  ( g*zp*zp ) - p0(i,j)/dp*exp( (r/rp)**(1.5) ) * exp( ( ztc(i,j,k)/zp )**(2.0) ) ) )**(0.5)
                           ttc(i,j,k) = (tv0(i,j)-gam*ztc(i,j,k))/(1.0+0.608*qtc(i,j,k))*(1.0+(2.0*rd*(tv0(i,j)-gam*ztc(i,j,k))*ztc(i,j,k))/  &
                                       (g*zp**(2.0)*(1.0-p0(i,j)/dp*exp((r/rp)**(1.5))*exp((ztc(i,j,k)/zp)**(2.0)))))**(-1.0)
                                          !t(phi,lam,z) dist for vortex and bckgrd, t approach bckgrd temp as r becomes large  
                      else
                         qtc(i,j,k) = qt       !dry model doesn't need it
                         ztc(i,j,k) = zt + rd*tvt(i,j)/g*log(pt(i,j)/pmod(i,j,k))
                         plev(i,j,k) = pt(i,j)*exp(-((g*(ztc(i,j,k)-zt))/(rd*tvt(i,j))))
                         ttc(i,j,k) = tvt(i,j)
                         vttc(i,j,k) = 0.0         !tangential wind
                      endif
                          utc(i,j,k) = vttc(i,j,k)*d1/d               !CHANGE THE CALC of d1, d2 in SH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          vtc(i,j,k) = vttc(i,j,k)*d2/d
                  enddo
                       pstc(i,j,1)= pmod(i,j,nz)/(1.0-ztc(i,j,nz)*gam/tv0(i,j))**(g/(rd*gam)) ! from inversing 
                enddo
           enddo
           filename='case7.cam2.h1.00'//trim((yrten))//trim(adjustl(x1))//'-'//trim((mnten))//trim(adjustl(x2))//'-'//trim((dyten))//trim(adjustl(x3))//'-00000.nc'
             status= nf90_open(filename,nf90_Write, ncID)
             call readnc_real2d('PSL'     ,psl,ncID,nx,ny)
             call readnc_real3d('U'       ,u,ncID,nx,ny,nz)
             call readnc_real3d('V'       ,v,ncID,nx,ny,nz)
             call readnc_real3d('T'       ,temp,ncID,nx,ny,nz)
           allocate(psl_tcincl(nx,ny),u_tcincl(nx,ny,nz),v_tcincl(nx,ny,nz),t_tcincl(nx,ny,nz))
           do i=1,nx
              do j=1,ny
                 psl_tcincl(i,j)=psl(i,j)
                 do k=1,nz
                    u_tcincl(i,j,k)=u(i,j,k)
                    v_tcincl(i,j,k)=v(i,j,k)
                    t_tcincl(i,j,k)=temp(i,j,k)
                 enddo
              enddo
           enddo
           ii=1
print*,xc-tcxgrid
print*,yc-tcygrid
           do i=xc-tcxgrid,xc+tcxgrid
              jj=1
              do j=yc-tcygrid,yc+tcygrid
                 psl_tcincl(i,j)=pstc(i,j,1)*1D2
                 do k=1,nz
                    u_tcincl(i,j,k)=utc(ii,jj,k)+ u_tcincl(i,j,k)
                    v_tcincl(i,j,k)=vtc(ii,jj,k)+ v_tcincl(i,j,k)
                    t_tcincl(i,j,k)=ttc(ii,jj,k)+ t_tcincl(i,j,k)
                 enddo
                 jj=jj+1
              enddo
              ii=ii+1
           enddo
           call writenc_real3d('U',u_tcincl,ncID,nx,ny,nz)
           call writenc_real3d('V',v_tcincl,ncID,nx,ny,nz)
           call writenc_real3d('T',t_tcincl,ncID,nx,ny,nz)
           status= nf90_close(ncID)
stop
             !allocate(divtc(size(lam),size(phi),nz),vortc(size(lam),size(phi),nz))
             !allocate(ux(size(lam),size(phi),nz),uy(size(lam),size(phi),nz),vx(size(lam),size(phi),nz),vy(size(lam),size(phi),nz))
             !call diffx1(ux,utc,dlam,size(lam),size(phi),nz)
             !call diffx1(vx,vtc,dlam,size(lam),size(phi),nz)
             !call diffy1(uy,vtc,dphi,size(lam),size(phi),nz)
             !call diffy1(vy,utc,dphi,size(lam),size(phi),nz)
             !divtc = ux + vy
             !vortc = vx - uy 
           if (fc .LT. 0) then  !should be careful when the TC is close to the EQ
             utc = -utc    !SH wind
             vtc = -vtc
           endif
             open(unit=110,file='plev.dat',form='formatted',status='unknown',action='write')
             open(unit=111,file='pstc.dat',form='formatted',status='unknown',action='write')
             open(unit=112,file='utc.dat',form='formatted',status='unknown',action='write')
             open(unit=113,file='vtc.dat',form='formatted',status='unknown',action='write')
             !open(unit=114,file='ztc.dat',form='formatted',status='unknown',action='write')
             open(unit=115,file='ttc.dat',form='formatted',status='unknown',action='write')
             !open(unit=116,file='divtc.dat',form='formatted',status='unknown',action='write')
             !open(unit=117,file='vortc.dat',form='formatted',status='unknown',action='write')
             !open(unit=118,file='vttc.dat',form='formatted',status='unknown',action='write')
             !open(unit=119,file='ztcini.dat',form='formatted',status='unknown',action='write')
             open(unit=120,file='lam.dat',form='formatted',status='unknown',action='write')
             open(unit=121,file='phi.dat',form='formatted',status='unknown',action='write')
                 do k=1,nz
                    write(110,'('// xsize //'f15.5)') plev(:,:,k) !u
                    !write(111,'('// xsize //'f15.5)') 2**0.5*log(pstc(:,:,1)/1000) !u
                    write(111,'('// xsize //'f15.5)') pstc(:,:,1) !u
                    write(112,'('// xsize //'f15.5)') utc(:,:,k) !u
                    write(113,'('// xsize //'f15.5)') vtc(:,:,k) !v
                    !write(114,'('// xsize //'f15.5)') ztc(:,:,k) !div
                    write(115,'('// xsize //'f15.5)') ttc(:,:,k) !temp
                    !write(116,'('// xsize //'f15.5)') divtc(:,:,k) !geopot
                    !write(117,'('// xsize //'f15.5)') vortc(:,:,k) 
                    !write(118,'('// xsize //'f15.5)') vttc(:,:,k) 
                    !write(119,'('// xsize //'f15.5)') ztcini(:,:,k) 
                 enddo
                    !write(120,'('// xsize //'f15.5)') lam/rad
                    !write(121,'('// xsize //'f15.5)') phi/rad
             close(110)
             close(111)
             close(112)
             close(113)
             !close(114)
             close(115)
             !close(116)
             !close(117)
             !close(118)
             !close(119)
             close(120)
             close(121)
           stop
           
           end program
