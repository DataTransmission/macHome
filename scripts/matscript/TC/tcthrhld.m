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
addpath('/bkirtman/gchen/data/ccsm4_0/archive/case7/atm')
addpath('/bkirtman/gchen/data/ccsm4_0/archive/case7/atm/TC_h2b_ncdata')
%var(1)='vort850';var(2)='spd';var(3)='z';var(4)='t';var(5)='tzint'
  filenamebsf= 'bsf.case7.cam2.h1.0008-01-02-00000.nc'
  filename= 'case7.cam2.h1.0008-01-02-00000.nc'
  fabsf= netcdf(filenamebsf,'nowrite');
  fa= netcdf(filename,'nowrite');
  %tem = nc_varget(filename,'T');
  phi  = fa{'lat'}(:); 
  thet = fa{'lon'}(:); 
  lev  = fa{'lev'}(:);
  tclev= fabsf{'lev_p'}(:);
  %hyai = nc_varget(filename,'hyai');
  %hybi = nc_varget(filename,'hybi');
  %pref = 1000;
  %plev = nc_varget(filename,'lev');
  close(fa);
  close(fabsf);
  nx= size(thet,1); ny= size(phi,1); nz= size(lev,1); nzbsf= size(tclev,1);
  nx1= ceil(nx*100/360); nx2= ceil(nx*160/360); 
  ny1= ceil(ny*90/180);  ny2= ceil(ny*130/180); 
  %tclev= [find((abs(300-plev)==(min(abs(300-plev))))),...  
  %       find((abs(500-plev)==(min(abs(500-plev))))),...
  %       find((abs(700-plev)==(min(abs(700-plev))))),... 
  %       find((abs(850-plev)==(min(abs(850-plev)))))];
  hfbox= 1; boxarea= (hfbox*2+1)^2; nmin=1; %halfbox grid #
  r=      6370*10^3; %earth radius
  dy=     (phi(2)- phi(1))* 110*10^3;
  dx=     (thet(2)- thet(1))* 110*10^3* cosd(phi);
  var=    {'vort850', 'sfspd','spd3_85','psl','ta3_5_7_85',};
  varvec= {'vort850vec','sfspdvec','tzintavec'};
  varthr= {'vort850thr','sfspdthr','tzintathr','psl'};
  bas=    {'ni', 'wnp','enp','atl',   'si','aus','sp'};
  basind= [  1,2,3,  -1,-2,-3,  0];
  day=    [31,28,31,30,31,30,31,31,30,31,30,31];
  nmon=   1:12;
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
  ii=1;
  for nb= 1:length(bas)
    if nb==3
      basmat(hem(nb,:),enp1) = nb; 
    elseif nb==4
      basmat(twenty5N:fortyN,atl1) = nb; 
      basmat(hem(nb,:),atl2) = nb; 
    else
      evalc([ 'basmat(hem(nb,:),' bas{nb} ') = nb' ]);
    end
     ind{nb} = find(basmat==nb) ;
     [indrow{nb}  indcol{nb} ]= find(basmat==nb);
     evalc([ 'basin=' bas{nb} ]);
     index(1,ii:ii+length(find(basmat==nb))-1)= nb;
     [index(2,ii:ii+length(find(basmat==nb))-1)  index(3,ii:ii+length(find(basmat==nb))-1)] = find(basmat==nb);
     ii=ii+ length(find(basmat==nb));
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
tic
  %could check variable by test = nc{'nx'}(:)
  for nmon=nmon 
    for nday = 1:day(nmon)
      if nmon< 10; mnten='0'; else mnten=''; end
      if nday< 10; dyten='0'; else dyten=''; end
      filename=(['case7.cam2.h1.0008-' mnten num2str(nmon) '-' dyten num2str(nday) '-00000.nc']);
      filenamebsf=(['bsf.case7.cam2.h1.0008-' mnten num2str(nmon) '-' dyten num2str(nday) '-00000.nc']);
      fa= netcdf(filename,'nowrite');
      fabsf= netcdf(filenamebsf,'nowrite');
      psl  = fa{'PSL'}(:,:); psl(psl>10^10)=NaN; %sea level pressure
      temp = fabsf{'T'}(2:4,:,:); temp(temp>10^10)=NaN; % temp at 300,500,700
      u    = fabsf{'U'}(1,:,:); u(u>10^10)=NaN; % u at 850
      v    = fabsf{'V'}(1,:,:); v(v>10^10)=NaN;
      us=fa{'U'}(1,nz,:,:); us(us>10^10)=NaN;
      vs=fa{'V'}(1,nz,:,:); vs(vs>10^10)=NaN;
%      spd  = sqrt(u.^2 + v.^2);
      close(fa)
      sfspd = sqrt(us.^2+vs.^2);
      vx(:,1) = (v(:,2) - v(:,1))./dx;
        for i = 2:nx-1
           vx(:,i) = (v(:,i+1) - v(:,i-1))./(2.*dx);
        end
      vx(:,nx) = (v(:,nx) - v(:,nx-1))./dx;
      uy(1,:) = (u(2,:) - u(1,:))./dy;
        for i = 2:ny-1
           uy(i,:) = (u(i+1,:) - u(i-1,:))./(2.*dy);
        end
      uy(ny,:) = (u(ny,:) - u(ny-1,:))./dy;
      vort850 = vx - uy;
%  contourf(thet,phi(2:length(phi)-1),squeeze(vort(26,2:length(phi)-1,:)),30,'linestyle','none');colorbar %avoid the pole Inf value
%==variable threshold===
      tzint= squeeze(sum(temp,1)); %300,500,700
      tzinta = zeros(ny,nx); tzintavg = zeros(ny,nx);
tic
      for i=1:length(index)
         tzintavg(index(2,i),index(3,i))= sum(sum(tzint(index(2,i)-hfbox:index(2,i)+hfbox,index(3,i)-hfbox:index(3,i)+hfbox)))/boxarea;
      end
      tzinta= tzint- tzintavg;
      if (nmon==1 | nmon==2 | nmon==3 | nmon==4 | nmon==12)   nbb=[5,6,7];
      elseif( nmon== 6 | nmon== 7 | nmon== 8 | nmon== 9 | nmon==10)    nbb=[1,2,3,4];
      else   nbb=1;
      end
      for nb=nbb
        evalc([ 'tzintavec.' bas{nb} '= cat(1,tzinta(find(basmat==nb)),tzintavec.' bas{nb} ')' ]);
        evalc([ 'vort850vec.' bas{nb} '= cat(1,vort850(find(basmat==nb)),vort850vec.' bas{nb} ')' ]);
        evalc([ 'sfspdvec.' bas{nb} '= cat(1,sfspd(find(basmat==nb)),sfspdvec.' bas{nb} ')' ]);
      end
toc
      if nmon <=11 & nmon >=6 %NH summer months
        vgl.nh= sum(sum(sfspd(hem(4,:),ni(1):atl2(end)))) + vgl.nh;
      else
        vgl.sh= sum(sum(sfspd(hem(5,:),si(1):sp(end) ))) + vgl.sh;
      end
    end
  end
toc

  ii=1;
  for i=1:length(tcRowInd)
      if temp(1,tcRowInd(i),tcColInd(i)) < 10e10
         tcrow(ii) = tcRowInd(i) ;
         tccol(ii) = tcColInd(i);
         ii=ii+1;
      end
  end

  TCParam={'nx','ny','nz','tclev','nx1','nx2','ny1','ny2','dy','dx','day','days',...
         'zeroN','zeroS','fortyN','fortyS','twenty5N','hfbox','boxarea'...
         'ni','wnp','enp','enp1','enp2','atl','atl1','atl2','si','aus','sp','tcrow','tccol'};
  for n=1:length(TCParam) % 32 variables saved
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
'Done Saving TCParam.nc'
  vglob= (vgl.nh+vgl.sh)/...
      (( length(ni(1):atl2(end))*length(hem(4,:))+...
         length(si(1):sp(end))*length(hem(5,:)) )* days);
  for nb=1:length(bas)
    evalc(['vort850thr.' bas{nb} '= 2* nanstd(vort850vec.' bas{nb} ')' ]);
    evalc(['sfspdthr.' bas{nb} '= vglob+ nanstd(sfspdvec.' bas{nb} ')']);
    evalc(['tzintathr.' bas{nb} '= nanstd(tzintavec.' bas{nb} ')' ]);
  end
%select the center region (7*7) of tzint for aave
%v is speed
%==form a map of threshold===
  thrmat={'vortthrmat','sfspdthrmat','tzintathrmat'};
  vortthrmat = zeros(ny,nx); sfspdthrmat = zeros(ny,nx); tzintathrmat = zeros(ny,nx);
  for nv=1:length(thrmat)
    for nb=1:length(bas)
      evalc([ thrmat{nv} '(find(basmat == nb))=' varthr{nv} '.' bas{nb} ]);
    end
  end 

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

%%%%USE FORTRAN SCRIPT TCDETECT.F90 TO FIND TC BY USING THE TCTHRHOLD MAPS!!!!!





