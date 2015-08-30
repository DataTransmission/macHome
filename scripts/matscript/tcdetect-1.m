%changing vgl criteria doesn't really make a differece
%changing min pressure criteria does make a difference
%changing stdev criteria for vort850 makes a sig difference, the
%the bigger criteria gives less connected track, so one stdev
%is better than the suggested 2stdev
%models of coarse resolution should have a lower maximum wind threshold
%because wind is harder to form  
%avg ... stdev for vort850, sfcspd, tzint
%vort binary data pattern
%e=1, t=1,  z=850, xy
%e=1, t=2,  z=850, xy
%...
%e=1, t=60, z=850, xy
%e=2, t=1,  z=850, xy
%e=2, t=2,  z=850, xy
%...
%e=41,t=60, z=850, xy
%==============parameter======================
% detect the grid exceed all variable necessary for TC definition, 
% then save the grid that exceeds this criteria
% ATM master fields:
% 'T','AEROD_v','CLDHGH','CLDICE','CLDLIQ','CLDLOW','CLDMED','CLDTOT','CLOUD','CMFDQ','CMFDQR','CMFDT','CMFMC','CMFMCDZM','CONCLD','DCQ','DTCOND','DTV','FICE','FLDS','FLDSC','FLNS','FLNSC','FLNT','FLNTC','FLUT','FLUTC','FREQSH','FREQZM','FSDS','FSDSC','FSNS','FSNSC','FSNT','FSNTC','FSNTOA','FSNTOAC','FSUTOA','GCLDLWP','ICEFRAC','ICLDIWP','ICLDTWP','LANDFRAC','LCLOUD','LHFLX','LWCF','NDROPCOL','NDROPMIX','NDROPSNK','NDROPSRC','OCNFRAC','ODV_SSLTA','ODV_SSLTC','ODV_bcar1','ODV_bcar2','ODV_dust1','ODV_dust2','ODV_dust3','ODV_dust4','ODV_ocar1','ODV_ocar2','ODV_sulf','OMEGA','OMEGAT','PBLH','PCONVB','PCONVT','PHIS','PRECC','PRECCDZM','PRECL','PRECSC','PRECSH','PRECSL','PRECT','PS','PSL','Q','QC','QFLX','QREFHT','QRL','QRS','RELHUM','RHREFHT','SFCLDICE','SFCLDLIQ','SHFLX','SNOWHICE','SNOWHLND','SOLIN','SRFRAD','SWCF','T','TAUX','TAUY','TGCLDIWP','TGCLDLWP','TMQ','TREFHT','TREFMNAV','TREFMXAV','TROP_P','TROP_T','TROP_Z','TS','TSMN','TSMX','U','US','UU','V','VD01','VQ','VS','VT','VU','VV','WTKE','Z3'

% thresholds:
%  pick summer season JJA for NH, DJF for SH, merge daily data into a singal matrix (Basinwise)
%  1) 850h-Pa vort > vort threshold (neg in SH)
%  2) 7x7 box max sfcspd > spd threshold 
%  3) min psl in 7x7
%  4) 7x7 tempa avg and three lev avg (300, 500, 700) > tempa threshold
%  5) 7x7 tempa avg over 300, 500, 700 > 0
%  6) 7x7 tempa avg over 300 > 850
%  7) 7x7 spd avg over 850 > 300
%  8) connect center of storm over time step (6hrly need 5.6 deg, daily need 8.5 deg)
%  9) defined as storm if criteria last > 2 day (1.5 day for 6hrly output) 
  %filename='case7.cam2.h1.0001-01-02-00000.nc'

% temp threshold:
% 1) calc daily vertically summed temp (DVST) (700 500 300 p-lev)
% 2) 7x7 box avg of DVST for each grid = avgDVST
% 3) DVSTa = DVST - avgDVST
% 4)  
run ~/matscript/startup.m
%var(1)='vort850';var(2)='spd';var(3)='z';var(4)='t';var(5)='tzint'
  filename= 'case7.cam2.h1.0001-01-02-00000.nc'
  temp = nc_varget(filename,'T');
  phi  = nc_varget(filename,'lat'); 
  thet  = nc_varget(filename,'lon'); 
  nx=size(temp,3); ny=size(temp,2); nz=size(temp,1); 
  nx1=ceil(nx*100/360); nx2=ceil(nx*160/360); ny1=ceil(ny*90/180); ny2=ceil(ny*130/180); 
  nt=60; nt1=1; nt2=60;  
  hfbox= 10; boxarea=(hfbox*2+1)^2; nmin=1; %halfbox grid #
  vx1=1; vx2=320; vy1=41; vy2=120;
  r=6370*10^3; %earth radius
  dy = (phi(2)- phi(1))*110*10^3;
  dx = (thet(2)- thet(1))*110*10^3*cosd(phi);
  var = {'vort850', 'sfspd'};
  varthr = {'vortthr','sfspdthr'};
  bas = {'ni' 'wnp' 'enp' 'atl' 'si' 'aus' 'sp'};
  vgl.nh= 0; vgl.sh= 0;
  mon = {'07' '08' '09' '12' '01' '02'}
  days = 31+31+30+31+31+30;
% NI(0:40N,40E:100E); WNP(0:40N,100E:160W); ENP(0:40N,160W:60W); ATL(0:40N,100W:20W)
% SI(0:40S,30E:110E); AUS(0:40S,110E:170E); SP(0:40S,170E:120W) 
  zeroN  = find(abs(0.01 -phi)== min(abs(0.01 -phi))); %to avoid having two index of eqtr(NH & SH), us 0.01 instead of 0
  zeroS  = find(abs(-0.01-phi)== min(abs(-0.01-phi))); 
  fortyN = find(abs(40   -phi)== min(abs(40   -phi)));
  fortyS = find(abs(-40  -phi)== min(abs(-40  -phi)));
  twenty5N= find(abs(25   -phi)== min(abs(25   -phi)));
%==NH==
  ni  = find(abs(40 -thet)==min(abs(40 -thet))):   find(abs(100-thet)==min(abs(100-thet)));
  wnp = find(abs(100-thet)==min(abs(100-thet)))+1: find(abs(200-thet)==min(abs(200-thet)));
  enp1= find(abs(200-thet)==min(abs(200-thet)))+1: find(abs(252-thet)==min(abs(252-thet)));
  enp2= find(abs(252-thet)==min(abs(252-thet)))+1: find(abs(291.5-thet)==min(abs(291.5-thet)));
  atl1= enp2; %there is an overlap upper trig is atl lower trig is enp 
  atl2= find(abs(291.5-thet)==min(abs(291.5-thet)))+1: find(abs(340-thet)==min(abs(340-thet)));
%==SH==
  si  = find(abs(30-thet)==min(abs(30 -thet))):    find(abs(110-thet)==min(abs(110-thet)));
  aus = find(abs(110-thet)==min(abs(110-thet)))+1: find(abs(170-thet)==min(abs(170-thet)));
  sp  = find(abs(170-thet)==min(abs(170-thet)))+1: find(abs(240-thet)==min(abs(240-thet)));
  zeromat = zeros(length(zeroN:twenty5N),length(enp2)); %create a matrix containing NaN at upper triangle and 0 at lower triangle.
  zeromat(:,:)=NaN;testl=zeromat;testu=zeromat;
  testl(find(isnan(tril(zeromat))))=0; 
  testu(find(isnan(triu(zeromat))))=0;
  for nv= 1:2
    for nb= 1:length(bas)
       eval([varthr{nv} '.' bas{nb} '= []']);
    end
  end

  for nmon = 1:6
    for nday = 1:31
      if nmon== 3 | nmon== 6 & nday== 31
      else
        if nday < 10
          filename=(['case7.cam2.h1.0001-' mon{nmon} '-0' num2str(nday) '-00000.nc']);
        else
          filename=(['case7.cam2.h1.0001-' mon{nmon} '-' num2str(nday) '-00000.nc']);
        end
      end
      temp = nc_varget(filename,'T');
      u    = nc_varget(filename,'U');
      v    = nc_varget(filename,'V');
      psl  = nc_varget(filename,'PSL'); %sea level pressure
      spd  = sqrt(u.^2 + v.^2);
      sfspd = squeeze(spd(26,:,:));
%  tempa = 
      tzint= temp; 
      for k= 1:nz
       vx(k,:,1) = (v(k,:,2)' - v(k,:,1)')./dx;
         for i = 2:nx-1
            vx(k,:,i) = (v(k,:,i+1)' - v(k,:,i-1)')./(2.*dx);
         end
       vx(k,:,nx) = (v(k,:,nx)' - v(k,:,nx-1)')./dx;
      end
      for k= 1:nz
       uy(k,1,:) = (squeeze(u(k,2,:))' - squeeze(u(k,1,:))')./dy;
         for i = 2:ny-1
            uy(k,i,:) = (squeeze(u(k,i+1,:))' - squeeze(u(k,i-1,:))')./(2.*dy);
         end
       uy(k,ny,:) = (squeeze(u(k,ny,:))' - squeeze(u(k,ny-1,:))')./dy;
      end
      vort = vx - uy;
%  contourf(thet,phi(2:length(phi)-1),squeeze(vort(26,2:length(phi)-1,:)),30,'linestyle','none');colorbar %avoid the pole Inf value
      vort850 = squeeze(vort(find(min(abs(850-nz))),:,:)); 
%==variable threshold===
      for nv= 1:2
        if nmon <= 3 %NH basin
          for nb= 1:3
            eval(['tempvar.' varthr{nv} '.' bas{nb} '=' varthr{nv} '.' bas{nb}]);
          end
          eval([varthr{nv} '.ni = (reshape(' var{nv} '(zeroN:fortyN,ni ),length(zeroN:fortyN)*length(ni),1))']);
          eval([varthr{nv} '.wnp= (reshape(' var{nv} '(zeroN:fortyN,wnp),length(zeroN:fortyN)*length(wnp),1))']);
        
          eval([varthr{nv} '.enp= reshape(' var{nv} '(zeroN:fortyN,enp1),length(zeroN:fortyN)*length(enp1),1)']); %reshape into a vector
          test.enp = eval([var{nv} '(zeroN:twenty5N,enp2)']) + testl; %finish creating a matix with lower triagle containing vort850 and NaN upper triangle.
          eval([varthr{nv} '.enp = cat(1,reshape(test.enp,size(test.enp,1)*size(test.enp,2),1),' varthr{nv} '.enp)']);

          eval([varthr{nv} '.atl= reshape(' var{nv} '(zeroN:fortyN,atl2),length(zeroN:fortyN)*length(atl2),1)']);
          eval([varthr{nv} '.atl= cat(1,reshape(' var{nv} '(twenty5N:fortyN,atl1),length(twenty5N:fortyN)*length(atl1),1),' varthr{nv} '.atl)']);
          test.atl = eval([var{nv} '(zeroN:twenty5N,atl1)']) + testu;
          eval([varthr{nv} '.atl = cat(1,reshape(test.atl,size(test.atl,1)*size(test.atl,2),1),' varthr{nv} '.atl)']);
          for nb= 1:3
            eval([varthr{nv} '.' bas{nb} '= cat(1,tempvar.' varthr{nv} '.' bas{nb} ',' varthr{nv} '.' bas{nb} ')']);
          end     
        else %SH basin
          for nb= 4:7
            eval(['tempvar.' varthr{nv} '.' bas{nb} '=' varthr{nv} '.' bas{nb}]);
          end
          eval([varthr{nv} '.si = (reshape(' var{nv} '(fortyS:zeroS,si ),length(fortyS:zeroS)*length(si),1))']);
          eval([varthr{nv} '.aus= (reshape(' var{nv} '(fortyS:zeroS,aus),length(fortyS:zeroS)*length(aus),1))']);
          eval([varthr{nv} '.sp = (reshape(' var{nv} '(fortyS:zeroS,sp ),length(fortyS:zeroS)*length(sp),1))']);
          for nb= 4:7
            eval([varthr{nv} '.' bas{nb} '= cat(1,tempvar.' varthr{nv} '.' bas{nb} ',' varthr{nv} '.' bas{nb} ')']);
          end     
        end
      end
      tempvgl.nh = vgl.nh; tempvgl.sh = vgl.sh;
      if nmon <= 3
        vgl.nh= (sum(sum(sfspd(zeroN:fortyN,ni )))+  sum(sum(sfspd(zeroN:fortyN,wnp)))+  sum(sum(sfspd(zeroN:fortyN,enp1)))+ ...
            sum(sum(sfspd(zeroN:fortyN,enp2)))+ sum(sum(sfspd(zeroN:fortyN,atl2)))) + tempvgl.nh;
      else
        vgl.sh= (sum(sum(sfspd(fortyS:zeroS,si )))+ sum(sum(sfspd(fortyS:zeroS,aus)))+ sum(sum(sfspd(fortyS:zeroS,sp)))) + tempvgl.sh;
      end
    end
  end
  vglsum= (vgl.nh+vgl.sh)/...
       ((length(ni)+length(wnp)+length(enp1)+length(enp2)+length(atl2))*length(zeroN:fortyN)+ (length(si)+length(aus)+length(sp))*length(fortyS:zeroS))/days;
  for nb=1:length(bas)
    eval(['sfspdthr.' bas{nb} '=vgl+sfspdthr.' bas{nb}]);
  end
  for nv=1:2
    for nb=1:length(bas)
      if nv==1
        eval([varthr{nv} '.' bas{nb} '= 2*nanstd(' varthr{nv} '.' bas{nb} ')']);
      else
        eval([varthr{nv} '.' bas{nb} '= nanstd(' varthr{nv} '.' bas{nb} ')']);
      end
    end
  end
pause
%select the center region (7*7) of tzint for aave
%v is speed
  ii=0;
%==form a map of threshold===
  thrmtrx={'vortthrmtrx','sfspdmtrx'};
  for nv=1
    eval([thrmtrx{nv} '= zeros(size(temp,2),size(temp,3))']);
    eval([thrmtrx{nv} '(zeroN:fortyN,ni) =' varthr{nv} '.ni']);
    eval([thrmtrx{nv} '(zeroN:fortyN,wnp)=' varthr{nv} '.wnp']);
    eval([thrmtrx{nv} '(fortyS:zeroS,si)='  varthr{nv} '.si']);
    eval([thrmtrx{nv} '(fortyS:zeroS,aus)=' varthr{nv} '.aus']);
    eval([thrmtrx{nv} '(fortyS:zeroS,sp)='  varthr{nv} '.sp']);
    
    eval([thrmtrx{nv} '(zeroN:fortyN,enp1)=' varthr{nv} '.enp']);
    zeromat = zeros(length(zeroN:twenty5N),length(enp2)); %create a matrix containing NaN at upper triangle and 0 at lower triangle.
    zeromat(:,:)=NaN;test.atl=zeromat;
    eval(['test.atl(find(isnan(triu(zeromat,1))))=' varthr{nv} '.atl']);
    eval(['test.atl(find(isnan(test.atl))) =' varthr{nv} '.enp']); %for enp2
    test.atl = flipud(test.atl);
    eval([thrmtrx{nv} '(zeroN:twenty5N,atl1)= test.atl']);
    eval([thrmtrx{nv} '(twenty5N:fortyN,atl1)=' varthr{nv} '.atl']);
    eval([thrmtrx{nv} '(zeroN:fortyN,atl2)=' varthr{nv} '.atl']);
  end 
   

  for nbas=1:7
%    nx1=eval([bas{nbas}]);ny1=eval([bas{
    for i=nx1
      for j=ny1
          condt = 0;
%          do k=1,2
%           if (temp(i,j) - sum(temp(j-hfbox:j+hfbox,i-hfbox:i+hfbox))/boxarea > 0) then
%               condt = condt + 1
%            endif
%          enddo

%          condz = 0
%          do k=1,nz
%          nii=0
%            do ni=i-hfbox,i+hfbox
%             do nj=j-hfbox,j+hfbox
%                 nii=nii+1
%                 zvmat(nii) = zv(ni,nj,t,e)
%              enddo
%            enddo
%            do r=1,nmin
%             %write(*,'(2f15.8)'),zv(i,j,1,t,e),minval(zvmat)
%              if (minval(zvmat) == zv(i,j,t,e)) then
%              else
%                call minimum(zvmat,boxarea,condz)%find the nmin-th minimum value for minimum zv value
%              endif
%            enddo
%          enddo
           if zv(i,j,1)==min(zv(i-hfbox:i+hfbox,j-hfbox:j+hfbox,1))  & ...  
               sum(spd(j-hfbox:j+hfbox,i-hfbox:i+hfbox,1))/boxarea > sum(spd(i-hfbox:i+hfbox,j-hfbox:j+hfbox,3))/boxarea 
               
%               condt == 2 .AND. ...
                
%               condz < nz*nmin .AND. ...
%               tv(i,j,2,t,e) - sum(tv(i-hfbox:i+hfbox,j-hfbox:j+hfbox,2,t,e))/boxarea > ...
%               tv(i,j,3,t,e) - sum(tv(i-hfbox:i+hfbox,j-hfbox:j+hfbox,3,t,e))/boxarea ...% .AND. ...
              
             ii=ii+1; 
             matrix(1,ii)= vort850(i,j,1);
             matrix(2,ii)= max(spd(i-hfbox:i+hfbox,j-hfbox:j+hfbox,1,));
             matrix(3,ii)= tzint(i,j,1) - sum(tzint(i-hfbox:i+hfbox,j-hfbox:j+hfbox,1))/boxarea; %vert int'd temp anomaly
             matrix(4,ii)= i;
             matrix(5,ii)= j;
             matrix(6,ii)= t;
             matrix(7,ii)= e;
         end
      end
    end

for e=1,ne
  for t=1,nt
    for k=1,nz
      vgl(1,1,t,e) = sum(spd(vx1:vx2,vy1:vy2,t,e))/((vx2-vx1+1)*(vy2-vy1+1));
%      write(110,'(320f15.8)') zv(:,:,t,e)
    end
  end
end
stdev(1) = 
call stdevsub(matrix(1,1:ii),stdev(1),ii)
call stdevsub(matrix(2,1:ii),stdev(2),ii)
call stdevsub(matrix(3,1:ii),stdev(3),ii)
call stdevsub(matrix(4,1:ii),stdev(4),ii)
write(*,'(4F13.8)'),stdev(1),stdev(2),stdev(3),stdev(4)
jj=0
k=1
open(unit=111,file='trackmap.dat',form='formatted',status='unknown',action='write')
do i=1,ii
  if ( ...
      matrix(1,i) > stdev(1) ... % .AND. ...                               %vort850
%      matrix(2,i) > sum(vgl(1,1,:,:))/(nt*ne) + stdev(2) .AND. ... %max(spd850)
%      matrix(3,i) > stdev(3) ...                                      %tzinta
     ) then
%write(*,'(2f15.8)'),matrix(3,i),stdev(3)
     write(111,'(4f15.8)') matrix(4,i),matrix(5,i),matrix(6,i),matrix(7,i)
  endif  
enddo
%open(unit=110,file='zdist.dat',form='formatted',status='unknown',action='write')
%close(110)
close(111)
stop

%====================================================================================
%====================================================================================
%====================================================================================
contains
subroutine readdat(varname,varout,nx,ny,nz,nt,ne)
implicit none
integer,intent(in) :: nx,ny,nz,nt,ne
integer :: e,t,irec
character(10), intent(in) :: varname
real, intent(out) :: varout(nx,ny,nz,nt,ne)
irec = 0 
open(111,file=trim(varname)//'.dat',form='unformatted',access='direct',recl=nx*ny)
do e=1,ne
  do t=1,nt
   do k=1,nz
      irec = irec+1
      read(111,rec=irec)varout(:,:,t,e)
   enddo
  enddo
enddo
close(111)
return
end subroutine readdat

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

subroutine stdevsub(varin,varout,ii)
implicit none
real,intent(in)::varin(1,ii)
real,intent(out)::varout
integer,intent(in)::ii
real::sumdiff2,diff2
sumdiff2 = 0
do i=1,ii
diff2 = (varin(1,i) - sum(varin(1,:))/ii)**2
sumdiff2 = diff2 + sumdiff2
enddo
varout = sqrt(sumdiff2/(ii-1))
return
end subroutine stdevsub


subroutine avg(varin,varout,nx,nx1,nx2,ny,ny1,ny2,nz,nt,ne,ind)
implicit none
integer,intent(in)::ind
real,intent(in) :: varin(nx,ny,nz,nt,ne)
integer,intent(in):: nx1,nx2,ny1,ny2
real,intent(out):: varout(nx,ny,nz,nt,ne)
integer,intent(inout)::nx,ny,nz,nt,ne
real ::rnx, rny, rnz, rnt, rne
integer :: i,j,t,e
rnx=1/real(nx); rny=1/real(ny); rnz=1/real(nz); rnt=1/real(nt); rne=1/real(ne)
select case (ind) 
  case (1) %area avg
nx=1;ny=1
  do e=1,ne
    do t=1,nt
      do k=1,nz
        varout(1,1,t,e) = sum(varin(nx1:nx2,ny1:ny2,t,e))*rnx*rny
      enddo
    enddo
  enddo
  case (2) %z avg
nz=1
  do e=1,ne
    do t=1,nt
      do i=1,nx
        do j=1,ny
          varout(i,j,1,t,e) = sum(varin(i,j,:,t,e))*rnz
        enddo
      enddo
    enddo
  enddo
  case (3) %t avg
nt=1
  do e=1,ne
    do i=1,nx
      do j=1,ny
        do k=1,nz
          varout(i,j,1,e) = sum(varin(i,j,:,e))*rnt
        enddo
      enddo
    enddo
  enddo
  case (4) %e avg
ne=1
  do t=1,nt
    do i=1,nx
      do j=1,ny
        do k=1,nz
          varout(i,j,t,1) = sum(varin(i,j,t,:))*rne
        enddo
      enddo
    enddo
  enddo
end select
return
end subroutine avg 
end program
