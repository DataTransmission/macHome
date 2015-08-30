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
% 1) calc daily vertically summed temp (DVST) (700 500 300 plev)
% 2) 7x7 box avg of DVST for each grid = avgDVST
% 3) DVSTa = DVST - avgDVST
% 4)  
run ~/matscript/startup.m
%var(1)='vort850';var(2)='spd';var(3)='z';var(4)='t';var(5)='tzint'
  filename= 'case7.cam2.h1.0002-01-02-00000.nc'
  temp = nc_varget(filename,'T');
  phi  = nc_varget(filename,'lat'); 
  thet = nc_varget(filename,'lon'); 
  plev = nc_varget(filename,'lev');
  nx= size(temp,3); ny= size(temp,2); nz= size(temp,1); 
  nx1= ceil(nx*100/360); nx2= ceil(nx*160/360); 
  ny1= ceil(ny*90/180);  ny2= ceil(ny*130/180); 
  tclev= [find((abs(300-plev)==(min(abs(300-plev))))),...  
         find((abs(500-plev)==(min(abs(500-plev))))),...
         find((abs(700-plev)==(min(abs(700-plev))))),... 
         find((abs(850-plev)==(min(abs(850-plev)))))];
  hfbox= 3; boxarea= (hfbox*2+1)^2; nmin=1; %halfbox grid #
  r=      6370*10^3; %earth radius
  dy=     (phi(2)- phi(1))* 110*10^3;
  dx=     (thet(2)- thet(1))* 110*10^3* cosd(phi);
  var=    {'vort850', 'sfspd','spd3_85','psl','ta3_5_7_85',};
  varvec= {'vort850vec','sfspdvec','tzintavec'};
  varthr= {'vort850thr','sfspdthr','tzintathr','psl'};
  bas=    {'ni' 'wnp' 'enp' 'atl' 'si' 'aus' 'sp'};
  mon=    {'07' '08' '09' '11' '01' '02'};
  day=    [31,31,30,30,31,28];
  inv3=  1.0/3.0;
  nbb= 1:7;
%  dayl= ;lunar year
  days = sum(day);
% NI(0:40N,40E:100E); WNP(0:40N,100E:160W); ENP(0:40N,160W:60W); ATL(0:40N,100W:20W)
% SI(0:40S,30E:110E); AUS(0:40S,110E:170E); SP(0:40S,170E:120W) 
  zeroN   = find(abs(0.01 -phi)== min(abs(0.01 -phi))); %to avoid having two index of eqtr(NH & SH), us 0.01 instead of 0
  zeroS   = find(abs(-0.01-phi)== min(abs(-0.01-phi))); 
  fortyN  = find(abs(40   -phi)== min(abs(40   -phi)));
  fortyS  = find(abs(-40  -phi)== min(abs(-40  -phi)));
  twenty5N= find(abs(25   -phi)== min(abs(25   -phi)));
%==NH==
  ni  = find(abs(40 -thet)==min(abs(40 -thet))):...
        find(abs(100-thet)==min(abs(100-thet)));
  wnp = find(abs(100-thet)==min(abs(100-thet)))+1:...
        find(abs(200-thet)==min(abs(200-thet)));
  enp = find(abs(200-thet)==min(abs(200-thet)))+1:...
        find(abs(291.5-thet)==min(abs(291.5-thet)));
  enp1= find(abs(200-thet)==min(abs(200-thet)))+1:...
        find(abs(252-thet)==min(abs(252-thet)));
  enp2= find(abs(252-thet)==min(abs(252-thet)))+1:...
        find(abs(291.5-thet)==min(abs(291.5-thet)));
  atl = find(abs(252-thet)==min(abs(252-thet)))+1:...
        find(abs(340-thet)==min(abs(340-thet))); 
  atl1= enp2; %there is an overlap upper trig is atl lower trig is enp 
  atl2= find(abs(291.5-thet)==min(abs(291.5-thet)))+1:...
        find(abs(340-thet)==min(abs(340-thet)));
%==SH==
  si  = find(abs(30-thet)==min(abs(30 -thet))):...
        find(abs(110-thet)==min(abs(110-thet)));
  aus = find(abs(110-thet)==min(abs(110-thet)))+1:...
        find(abs(170-thet)==min(abs(170-thet)));
  sp  = find(abs(170-thet)==min(abs(170-thet)))+1:...
        find(abs(240-thet)==min(abs(240-thet)));


  for i=nbb
    if i<=4
      hem(i,:)= zeroN:fortyN;
    else
      hem(i,:)= fortyS:zeroS;
    end
  end
  zeromat = zeros(length(zeroN:twenty5N),length(enp2)); %create a matrix containing NaN at upper triangle and 0 at lower triangle.
  zeromat(:,:)= NaN; trimat= zeromat; 
  trimat(find(isnan(tril(zeromat,1))))= 3;
  trimat(find(isnan(triu(trimat,1))))= 4;
  trimatflip=flipud(trimat); 
    % 3 3 3 4    
    % 3 3 4 4    this is trimatflip
    % 3 4 4 4
  basmat = zeros(ny,nx);
  basmat(zeroN:twenty5N,enp2) = trimatflip; % plug in trimatflip
  for nb= 1:length(bas)
    if nb==3
      evalc([ 'basmat(hem(nb,:),enp1) = nb' ]); 
      ind{nb} = find(basmat==nb) ;
      [indrow{nb}  indcol{nb}] = find(basmat==nb);
    elseif nb==4
      evalc([ 'basmat(twenty5N:fortyN,atl1) = nb' ]); 
      evalc([ 'basmat(hem(nb,:),atl2) = nb' ]); 
      ind{nb} = find(basmat==nb) ;
      [indrow{nb}  indcol{nb} ]= find(basmat==nb);
    else
      evalc([ 'basmat(hem(nb,:),' bas{nb} ') = nb' ]);
      ind{nb} = find(basmat==nb) ;
      [indrow{nb}  indcol{nb} ]= find(basmat==nb);
    end
  end
  
  basmat(find(basmat==0))=NaN;
  for nv= 1:3
    for nb= 1:length(bas)
       evalc([varvec{nv} '.' bas{nb} '= []']);
    end
  end
  vgl.nh= 0; vgl.sh= 0;
  tcRowInd=cat(1,indrow{:}) ;
  tcColInd=cat(1,indcol{:}) ;
  % Open netCDF file.
  TCParam={'nx','ny','nz','tclev','nx1','nx2','ny1','ny2','dy','dx','day','days',...
         'zeroN','zeroS','fortyN','fortyS','twenty5N','hfbox','boxarea'...
         'ni','wnp','enp','enp1','enp2','atl','atl1','atl2','si','aus','sp','tcRowInd','tcColInd'};
  for n=1:length(TCParam)
    dim(n)=length(eval([TCParam{n}]));
  end
  nc = netcdf('TCParam.nc','clobber');% 'clobber'); 
  nc.description ='Variables contain basin indices and TC restricted region';
  nc.source='CCSM4 daily output';
  nc.resolution =([num2str(thet(2)-thet(1)) ' lon ' num2str(phi(2)-phi(1))  ' lat']);
 %nc.date ='Apr, 2010';
  for n=1:length(TCParam)
    nc(TCParam{n})=dim(n)
    nc{TCParam{n}}=TCParam{n}
    nc{TCParam{n}}(:)=eval([TCParam{n}])
  end
  close(nc)
pp
pause
  %could check variable by test = nc{'nx'}(:)
tic
  for nmon = 1:6
nmon

    for nday = 1:day(nmon)
tic
nday
      if nday < 10
        filename=(['case7.cam2.h1.0002-' mon{nmon} '-0' num2str(nday) '-00000.nc']);
      else
        filename=(['case7.cam2.h1.0002-' mon{nmon} '-' num2str(nday) '-00000.nc']);
      end
      temp = nc_varget(filename,'T');
      u    = nc_varget(filename,'U');
      v    = nc_varget(filename,'V');
      psl  = nc_varget(filename,'PSL'); %sea level pressure
%      spd  = sqrt(u.^2 + v.^2);
      sfspd = sqrt(squeeze(u(nz,:,:)).^2+squeeze(v(nz,:,:)).^2);
      for k= tclev(4) % at 850
       vx(:,1) = (v(k,:,2)' - v(k,:,1)')./dx;
         for i = 2:nx-1
            vx(:,i) = (v(k,:,i+1)' - v(k,:,i-1)')./(2.*dx);
         end
       vx(:,nx) = (v(k,:,nx)' - v(k,:,nx-1)')./dx;
      end
      for k= tclev(4) 
       uy(1,:) = (squeeze(u(k,2,:))' - squeeze(u(k,1,:))')./dy;
         for i = 2:ny-1
            uy(i,:) = (squeeze(u(k,i+1,:))' - squeeze(u(k,i-1,:))')./(2.*dy);
         end
       uy(ny,:) = (squeeze(u(k,ny,:))' - squeeze(u(k,ny-1,:))')./dy;
      end
      vort850 = vx - uy;
%  contourf(thet,phi(2:length(phi)-1),squeeze(vort(26,2:length(phi)-1,:)),30,'linestyle','none');colorbar %avoid the pole Inf value
%==variable threshold===
      tzint= squeeze(sum(temp(tclev(1:3),:,:),1)); %300,500,700
      tzinta = zeros(ny,nx);tzintavg = zeros(ny,nx); 
      for nv= 1:2
        for nb= nbb
          evalc([ varvec{nv} '.' bas{nb} '= cat(1,' varvec{nv} '.' bas{nb} ',' var{nv} '(ind{nb}))' ]);
        end
      end
tic
      for nb=nbb
        for j=1:length(indrow{nb})
          for i=1:length(indcol{nb})
            tzintavg(indrow{nb}(j),indcol{nb}(i))= sum(sum(tzint(indrow{nb}(j)-hfbox:indrow{nb}(j)+hfbox,indcol{nb}(i)-hfbox:indcol{nb}(i)+hfbox)))/boxarea;
          end
        end
        tzinta= tzint- tzintavg;
        evalc([ 'tzintavec.' bas{nb} '= cat(1,tzinta(find(basmat==nb)),tzintavec.' bas{nb} ')' ]);
      end
toc
      if nmon <= 3 %NH has 1st three summer months
        vgl.nh= sum(sum(sfspd(hem(4,:),ni(1):atl2(end)))) + vgl.nh;
      else
        vgl.sh= sum(sum(sfspd(hem(5,:),si(1):sp(end) ))) + vgl.sh;
      end
toc
    end
  end
toc

  vglob= (vgl.nh+vgl.sh)/...
      (( length(ni(1):atl2(end))*length(hem(4,:))+...
         length(si(1):sp(end))*length(hem(5,:)) )* days);
  for nb=nbb
    evalc(['vort850thr.' bas{nb} '= 2* std(vort850vec.' bas{nb} ')' ]);
    evalc(['sfspdthr.' bas{nb} '= vglob+ std(sfspdvec.' bas{nb} ')']);
    evalc(['tzintathr.' bas{nb} '= std(tzintavec.' bas{nb} ')' ]);
  end
%select the center region (7*7) of tzint for aave
%v is speed
%==form a map of threshold===
  thrmat={'vortthrmat','sfspdthrmat','tzintathrmat'};
  vortthrmat = zeros(ny,nx); sfspdthrmat = zeros(ny,nx); tzintathrmat = zeros(ny,nx);
  for nv=1:length(thrmat)
    for nb=nbb
      evalc([ thrmat{nv} '(find(basmat == nb))=' varthr{nv} '.' bas{nb} ]);
    end
  end 
%  for nv=1:3
%    for nb=nbb
%      fid = fopen(eval([ thrmat{nv} '.' bas{nb} '.dat' ]); 
%      fwrite(fid,eval([ thrmat{nv} '.' bas{nb} ]);
%    end
%  end

%save('var.mat','basmat','indrow','indcol','ind','vortthr','sfspdthr','tzintathr',...
%     'vortthrmat','sfspdthrmat','tzintathrmat','-append') 

  nc = netcdf('TCThrMat.nc','clobber');% 'clobber'); 
  nc.description ='Three Threshold Matrix with values in TC defined region';
  nc.source='CCSM4 daily output';
  nc.resolution =([num2str(thet(2)-thet(1)) ' lon ' num2str(phi(2)-phi(1))  ' lat']);
  nc('nx')= nx; %define dim
  nc('ny')= ny;
  for n=1:length(thrmat)
    nc{thrmat{n}}={'ny','nx'} %define dim for matrix
    nc{thrmat{n}}=thrmat{n}
    nc{thrmat{n}}(:)=eval([thrmat{n}])
  end


%============================================================================
%===========================DETECT TC========================================
%============================================================================

%  for nmon = 1:6
%    for nday = 1:day(nmon)
%      if nday < 10
%        filename=(['case7.cam2.h1.0001-' mon{nmon} '-0' num2str(nday) '-00000.nc']);
%      else
%        filename=(['case7.cam2.h1.0001-' mon{nmon} '-' num2str(nday) '-00000.nc']);
%      end
%      temp = nc_varget(filename,'T');
%      u    = nc_varget(filename,'U');
%      v    = nc_varget(filename,'V');
%      psl  = nc_varget(filename,'PSL'); %sea level pressure
%%      spd  = sqrt(u.^2 + v.^2);
%      sfspd = sqrt(squeeze(u(nz,:,:)).^2+squeeze(v(nz,:,:)).^2);
%      for k= plev(4) % at 850
%       vx(:,1) = (v(k,:,2)' - v(k,:,1)')./dx;
%         for i = 2:nx-1
%            vx(:,i) = (v(k,:,i+1)' - v(k,:,i-1)')./(2.*dx);
%         end
%       vx(:,nx) = (v(k,:,nx)' - v(k,:,nx-1)')./dx;
%      end
%      for k= plev(4) 
%       uy(1,:) = (squeeze(u(k,2,:))' - squeeze(u(k,1,:))')./dy;
%         for i = 2:ny-1
%            uy(i,:) = (squeeze(u(k,i+1,:))' - squeeze(u(k,i-1,:))')./(2.*dy);
%         end
%       uy(ny,:) = (squeeze(u(k,ny,:))' - squeeze(u(k,ny-1,:))')./dy;
%      end
%      vort850 = vx - uy;
%      for j=fortyS-hfbox: fortyN+hfbox
%        for i=si(1)-hfbox: atl2(end)+hfbox
%          for k=1:length(plev)
%            tavg(k,j,i)= mean2(temp(plev(k),j-hfbox:j+hfbox,i-hfbox:i+hfbox));
%            ta(k,j,i)= temp(plev(k),j,i)- tavg(k,j,i);
%          end
%        end
%      end
%      for j=fortyS:fortyN
%        for i=si(1):atl2(end)
%          for k=1:length(plev) 
%            taavg(k,j,i)= mean2(ta(k,j-hfbox:j+hfbox,i-hfbox:i+hfbox));
%          end
%            if ...
%            vort 
%            (mean2(spd(plev(4),j-hfbox:j+hfbox,i-hfbox:i+hfbox))- ...
%             mean2(spd(plev(1),j-hfbox:j+hfbox,i-hfbox:i+hfbox)))> 0);
%        end
%      end
%   end
% end  
%
%    dd=dd+1;
%for nv=1:3
%  for nb=nbb
%    fid = fopen(eval([ varthrmat{nv} '.' bas{nb} '.dat' ]); 
%    fwrite(fid,eval([ varthrmat{nv} '.' bas{nb} ]);
%  end
%end
%  for nhem=1:2
%  for nbas=1:7
%    nx1=evalc([bas{nbas}]);ny1=hem(nhem);
%    for i=nx1
%      for j=ny1
%          condt = 0;
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
%           if zv(i,j,1)==min(zv(i-hfbox:i+hfbox,j-hfbox:j+hfbox,1))  & ...  
%               sum(spd(j-hfbox:j+hfbox,i-hfbox:i+hfbox,1))/boxarea > sum(spd(i-hfbox:i+hfbox,j-hfbox:j+hfbox,3))/boxarea 
               
%               condt == 2 .AND. ...
                
%               condz < nz*nmin .AND. ...
%               tv(i,j,2,t,e) - sum(tv(i-hfbox:i+hfbox,j-hfbox:j+hfbox,2,t,e))/boxarea > ...
%               tv(i,j,3,t,e) - sum(tv(i-hfbox:i+hfbox,j-hfbox:j+hfbox,3,t,e))/boxarea ...% .AND. ...
              
%             ii=ii+1; 
%             matrix(1,ii)= vort850(i,j,1);
%             matrix(2,ii)= max(spd(i-hfbox:i+hfbox,j-hfbox:j+hfbox,1,));
%             matrix(3,ii)= tzint(i,j,1) - sum(tzint(i-hfbox:i+hfbox,j-hfbox:j+hfbox,1))/boxarea; %vert int'd temp anomaly
%             matrix(4,ii)= i;
%             matrix(5,ii)= j;
%             matrix(6,ii)= t;
%             matrix(7,ii)= e;
%         end
%      end
%    end
%
%for e=1,ne
%  for t=1,nt
%    for k=1,nz
%      vgl(1,1,t,e) = sum(spd(vx1:vx2,vy1:vy2,t,e))/((vx2-vx1+1)*(vy2-vy1+1));
%      write(110,'(320f15.8)') zv(:,:,t,e)
%    end
%  end
%end
%stdev(1) = 
%call stdevsub(matrix(1,1:ii),stdev(1),ii)
%call stdevsub(matrix(2,1:ii),stdev(2),ii)
%call stdevsub(matrix(3,1:ii),stdev(3),ii)
%call stdevsub(matrix(4,1:ii),stdev(4),ii)
%write(*,'(4F13.8)'),stdev(1),stdev(2),stdev(3),stdev(4)
%jj=0
%k=1
%open(unit=111,file='trackmap.dat',form='formatted',status='unknown',action='write')
%do i=1,ii
%  if ( ...
%      matrix(1,i) > stdev(1) ... % .AND. ...                               %vort850
%%      matrix(2,i) > sum(vgl(1,1,:,:))/(nt*ne) + stdev(2) .AND. ... %max(spd850)
%%      matrix(3,i) > stdev(3) ...                                      %tzinta
%     ) then
%write(*,'(2f15.8)'),matrix(3,i),stdev(3)
%     write(111,'(4f15.8)') matrix(4,i),matrix(5,i),matrix(6,i),matrix(7,i)
%  endif  
%enddo
%open(unit=110,file='zdist.dat',form='formatted',status='unknown',action='write')
%close(110)
%close(111)
%stop

%====================================================================================
%====================================================================================
%====================================================================================

%subroutine minimum(zvmat,boxarea,condz)
%implicit none          
%integer,intent(inout) :: condz, boxarea
%real,intent(inout) :: zvmat(boxarea)
%integer :: ni
%real :: minzvmat
%minzvmat=minval(zvmat)
%do ni=1,boxarea
%  if (zvmat(ni) /= minzvmat) then
%  else
%    zvmat(ni) = 10e5 
%  endif
%enddo
%condz = condz + 1
%return
%end subroutine minimum
%
%subroutine stdevsub(varin,varout,ii)
%implicit none
%real,intent(in)::varin(1,ii)
%real,intent(out)::varout
%integer,intent(in)::ii
%real::sumdiff2,diff2
%sumdiff2 = 0
%do i=1,ii
%diff2 = (varin(1,i) - sum(varin(1,:))/ii)**2
%sumdiff2 = diff2 + sumdiff2
%enddo
%varout = sqrt(sumdiff2/(ii-1))
%return
%end subroutine stdevsub
%

