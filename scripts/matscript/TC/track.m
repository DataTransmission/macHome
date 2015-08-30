
condition= input('Detect tracks? 1;  Track Field Anomalies? 2')
switch condition
case 1
clear
run ~/matscript/startup.m
addpath('/bkirtman/gchen/data/ccsm4_0/archive/case7/ocn')
%addpath('/bkirtman/gchen/data/ccsm4_0/archive/case7/ocn/old')
addpath('/bkirtman/gchen/data/ccsm4_0/archive/case7/atm')
%addpath('/bkirtman/gchen/data/ccsm4_0/archive/case7/atm/old')
  zeroS= nc_varget('TCParam.nc','zeroS')
  zeroN= nc_varget('TCParam.nc','zeroN')
  afname= 'case7.cam2.h1.0008-01-02-00000.nc'; % this file just to get phi thet plev
  af= netcdf(afname,'nowrite');
  temp = af{'T'}(:,:,:);
  alat  = af{'lat'}(:);
  alon = af{'lon'}(:);
  ofname= 'case7.pop.h.nday1.0005-07.nc';
  of= netcdf(ofname,'nowrite');
  olat = of{'ULAT'}(:,:);
  olon= of{'ULONG'}(1,:);
  ng=2; % plotting grid numbers
  olatng=10;
  olonng=20;
  alatng=round(olatng*(olat(2)-olat(1))/(alat(2)-alat(1)))
  alonng=round(olonng*(olon(2)-olon(1))/(alon(2)-alon(1)))
  tsng=10; % ts avg'ing grid numbers
  ocol_ind_ext(1:olonng)=length(olon)-olonng+1:length(olon); 
  ocol_ind_ext(olonng+1:olonng+length(olon))=1:length(olon); 
  ocol_ind_ext(end+1:end+olonng)=1:olonng; % extended longitude
  acol_ind_ext(1:alonng)=length(alon)-alonng+1:length(alon); 
  acol_ind_ext(alonng+1:alonng+length(alon))=1:length(alon); 
  acol_ind_ext(end+1:end+alonng)=1:alonng; % extended longitude
  plev = of{'lev'}(:);
  z_t  = of{'z_t'}(:);
  close(af);close(of);
  %tclev= [find((abs(300-plev)==(min(abs(300-plev))))),...
         %find((abs(500-plev)==(min(abs(500-plev))))),...
         %find((abs(700-plev)==(min(abs(700-plev))))),...
         %find((abs(850-plev)==(min(abs(850-plev)))))];
  nx= size(temp,3); ny= size(temp,2); %nz= size(temp,1); %not at pressure level
  dy= (alat(2)- alat(1))* 110*10^3;
  dx= (alon(2)- alon(1))* 110*10^3* cosd(alat);
  nz1=1;
  nz2=10;
  tcmp=load('tcmp.dat')'; %tcmp(i,:) 1:yr 2:mon 3:day 4:latind 5:lonind 
  nlon=0;
  F=1; % fraction of heat transported below mixed layer
  rho0=1020; % kgm-3, ocn density
  cp=3900; % Jkg-1C-1, ocn heat capacity 
  dW=400*10^3; % average TC width
  dL=100*10^3; % average TC travels 3 grids, I'm guessing from the output
%  linestyle={'.','<','+','o','x','s','d','^','v','>','h'};
%  colors={'r','g','b','m','c'};
  mon=    {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
  mday=    [31,28,31,30,31,30,31,31,30,31,30,31];  
  ii=1; %ii day of tracks
  jj=1;
  tdy=1;
  trckng=3; %TC within n grid points afr are counted as a track
  nd=0;
  yr=unique(tcmp(1,:));
  mn=unique(tcmp(2,:));
  day=unique(tcmp(3,:));
  saveday=-10:10;
  ii=1;  dayvec=[]; md1=saveday(end)+1;  md2=md1+mday(mn(1))
  if mn(ii)==1; mnone=12; else mnone=mn(ii)-1; end
  dayvec(1,1:saveday(end))=mnone; 
  dayvec(2,1:saveday(end))=mday(mnone)-saveday(end)+1:mday(mnone); % add Dec 29 30 31 to start of dayvec
  while (ii <=length(mn)-1);
    dayvec(1,md1:md1+mday(mn(ii))-1)=mn(ii);  % iith months
    dayvec(2,md1:md1+mday(mn(ii))-1)=1:mday(mn(ii));  % 1:last day of iith month
    md1=md1+mday(mn(ii));
    ii=ii+1;  
  end
  dayvec(1,md1:md1+mday(mn(ii))-1)=mn(ii);
  dayvec(2,md1:md1+mday(mn(ii))-1)=1:mday(mn(ii));  
  if mn(ii)==12; mnone=1; else mnone=mn(ii)+1; end
  dayvec(1,md1+mday(mn(ii)):md1+mday(mn(ii))+saveday(end)-1)=mnone;  % add Jan 1 2 3 to the end of dayvec
  dayvec(2,md1+mday(mn(ii)):md1+mday(mn(ii))+saveday(end)-1)=1:saveday(end);
%================================================================================================================
%
%
%                        Detect the consecutive days and points and save the points as a tracks 
%
%
%================================================================================================================
  ii=1; %count track
  for yy=1:length(yr)
     mm=1; %count month
     dd=saveday(end)+1; %count day
     while dd <= length(dayvec)-1-2*saveday(end)
        if (dayvec(2,dd)==mday(mm));  onemore=1;  else onemore=0;  
        end  %if last day of a month, then go to first day of next month to compare criteria
        itmp1= find(tcmp(3,:)==dayvec(2,dd) & tcmp(2,:)==mn(mm) & tcmp(1,:)==yr(yy)); %  index of same date TC 1st day
        itmp2= find(tcmp(3,:)==dayvec(2,dd+1) & tcmp(2,:)==mn(mm)+onemore & tcmp(1,:)==yr(yy)); %  index of same date TC 2nd day
        if (onemore==1);  mm=mm+1;  
        end
        if (isempty(itmp1)==0   &   isempty(itmp2)==0)
           dayct2=1;   % TC points at conti two days  
           while (dayct2<=length(itmp2));  dayct1=1;  % march one point at a time to compare with all value in 1st day 
              while (dayct1<=length(itmp1))  % fix 2nd day point and compare with all point in 1st day once track found break the loop and go to 2nd day 2nd point
                 if (abs(tcmp(4,itmp2(dayct2))- tcmp(4,itmp1(dayct1)))<= trckng)   &   (abs(tcmp(5,itmp2(dayct2))- tcmp(5,itmp1(dayct1)))<= trckng)
                    trckmp(:,ii)= tcmp(:,itmp1(dayct1));
                    trckmp(:,ii+1)= tcmp(:,itmp2(dayct2));
                    ii=ii+2;
                    if ii>3 & trckmp(:,ii-2)==trckmp(:,ii-3) %remove repeated points
                      trckmp(:,ii-2)=trckmp(:,ii-1);
                      ii=ii-1;
                    end
                    dayct1=inf;dayct2=inf;
                 end;  dayct1= dayct1+1;
              end;  dayct2= dayct2+1;
           end
        end;  dd=dd+1;
     end
  end
  pp=1;ii=1;
%=======================================================================
%
%           label number of tracks
%
%=======================================================================
  for i=1:size(trckmp,2)-1
     if trckmp(3,i+1)-trckmp(3,i) == 1   &  abs(trckmp(4,i+1)-trckmp(4,i)) <=ng  & abs(trckmp(5,i+1)-trckmp(5,i)) <=ng...
          |  trckmp(3,i+1) - trckmp(3,i) + mday(trckmp(2,i)) + trckmp(1,i+1) - trckmp(1,i)== 1 
        trckmp(6,i)=pp;
     else 
        trckmp(6,i)=pp;
        pp=pp+1;ii=1;
     end
  end
  trckmp(6,end)=pp;


  disp('Done detecting tracks')

case 2
%================================================================================================================
%
%
%                    Detect the tracks that crosses ocn and save the fields (TS,PSL,HMXL)
%
%
%================================================================================================================
  yrstr={'000','00','0',''}; ii=1;
  sshcorr=[];pslcorr=[];
  for i=1:size(trckmp,2)
      if trckmp(1,i)< 10; yrten=yrstr{1}; elseif trckmp(1,i)<100; yrten=yrstr{2}; elseif trckmp(1,i)<1000; yrten=yrstr{3}; else yrten=yrstr{4}; end
      if trckmp(2,i)< 10; mnten='0'; else mnten=''; end
      if trckmp(3,i)< 10; dyten='0'; else dyten=''; end
      afname=(['case7.cam2.h1.' yrten num2str(trckmp(1,i)) '-' mnten num2str(trckmp(2,i)) '-' dyten num2str(trckmp(3,i)) '-00000.nc']);
      ofname=(['case7.pop.h.nday1.' yrten num2str(trckmp(1,i)) '-' mnten num2str(trckmp(2,i)) '.nc']);
      ofnamexx=(['case7.pop.h.nday1.' yrten num2str(trckmp(1,i)-1) '-' mnten num2str(trckmp(2,i)) '.nc']);
      afnamexx=(['case7.cam2.h0.' yrten num2str(trckmp(1,i)-1) '-' mnten num2str(trckmp(2,i)) '.nc']);
      af= netcdf(afname,'nowrite');
      of= netcdf(ofname,'nowrite');
%      lhflx=  af{'LHFLX'}(:,:);
%      psl=    af{'PSL'}(:,:);
%      sst=    of{'SST'}(trckmp(3,i),:,:);
%      octemp= of{'TEMP'}(trckmp(3,i),nz1:nz2,:,:);
      hmxl=   of{'HMXL_2'}(trckmp(3,i),:,:);  %find the grid it's on and take a fix lat zonal plot of depth
%      ssh=    of{'SSH'}(:,:,:);
      [r r1]= find(abs(olat  -alat(trckmp(4,i)))==min(min(abs(olat-  alat(trckmp(4,i)))))); r=r(1);
       c    = find(abs(olon-alon(trckmp(5,i)))== ...
              min(min(abs(olon-alon(trckmp(5,i))))));   %the longitude goes to zero at the pole, so just pick a fixed lat
%=============================================
%
%     Save dTS, dPSL, dHMXL, Q in trckmp(6:8,:)
%    (Save the 1st day of track field and subtract the following day with the 1st day)
%=============================================
      if hmxl(r,c)< 9e+10  %not NaN
%         hmxlx2= of{'HMXL_2'}(:,:,:);
%         hmxlx2= squeeze(mean(hmxlx2,1)); %monthly mean 
%         hmxlx2(hmxlx2>10^10)=NaN; 
         ocols=ocol_ind_ext(c:c+2*olonng);
         orows=r-olatng:r+olatng;
         arows=trckmp(4,i)-alatng:trckmp(4,i)+alatng; 
         acols=acol_ind_ext(trckmp(5,i):trckmp(5,i)+2*alonng);
         trckmp2(1:6,ii)= trckmp(1:6,i); % 6 track #,  7 saveday
         trckmp2(7,ii)=   alat(trckmp(4,i));
         trckmp2(8,ii)=   alon(trckmp(5,i));
         ii=ii+1;
%         for sd=1:length(saveday);
            disp(saveday(sd))
            itmp=find(dayvec(1,saveday(end)+1:end-saveday(end))==trckmp(2,i)  &  dayvec(2,saveday(end)+1:end-saveday(end))==trckmp(3,i));
            itmp=itmp+saveday(end); %avoid 1st few index that are end of Dec
            if (itmp+saveday(sd)> saveday(end)); one=0; elseif (itmp+saveday(sd)> length(dayvec)-saveday(end)); one=1; else one=-1;  end 
          %end % less than Jan 3
            yyy= trckmp(1,i)+one; % less than Jan 1 if one=-1, more than Dec 31 if one=1
            mmm= dayvec(1,itmp+saveday(sd)); % yyy,mmm,ddd are for ocn file use. find 1 less index, could be same month could be one less 
            ddd= dayvec(2,itmp+saveday(sd));
            % last year same month avg
            if yyy< 10; yrten=yrstr{1}; elseif yyy<100; yrten=yrstr{2}; elseif yyy<1000; yrten=yrstr{3}; else yrten=yrstr{4};  end
            if mmm< 10; mnten='0'; else mnten='';  end
            if ddd< 10; dyten='0'; else dyten='';  end
            if yyy==trckmp(1,i) & mmm==trckmp(2,i) % if 1 less index in the same year and month, just use the same ocn.nc file
%               hmxlx=   of{'HMXL_2'}(ddd,:,:); 
%               sshx=    of{'SSH'}(ddd,:,:);
%               sstx=    of{'SST'}(ddd,:,:); 
%               octempx= of{'TEMP'}(ddd:ddd+2,nz1:nz2,:,:);
            else
               ofnamex=(['case7.pop.h.nday1.' yrten num2str(yyy) '-' mnten num2str(mmm) '.nc']);
               ofx= netcdf(ofnamex,'nowrite');
%               hmxlx=    ofx{'HMXL_2'}(ddd,:,:);  %find the grid it's on and take a fix lat zonal plot of depth
%               sshx=    ofx{'SSH'}(ddd,:,:);
%               sstx=    ofx{'SST'}(ddd,:,:);
%               octempx= of{'TEMP'}(ddd,ddd+2,nz1:nz2,:,:);
               close(ofx)
            end
            %afnamex=(['case7.cam2.h1.' yrten num2str(yyy) '-' mnten num2str(mmm) '-' dyten num2str(ddd) '-00000.nc']); %back one day
            %afx=netcdf(afnamex,'nowrite');
%            ux= afx{'U'}(:,:,:);
%            vx= afx{'V'}(:,:,:);
            %close(afx);
%            yyy1=trckmp(1,i); 
%            mmm1=trckmp(2,i);
%            if yyy1< 10; yrten=yrstr{1}; elseif yyy1<100; yrten=yrstr{2}; elseif yyy1<1000; yrten=yrstr{3}; else yrten=yrstr{4}; end
%            if mmm1-1< 10; mnten='0'; else mnten='';            end
%            if mmm1==1; one=-1; mmm1=12; mnten=''; else one=0; mmm1=mmm1-1; end  %mmm1=12 if mmm1=Jan, go one less year to mmm1=Dec
%            ofnamex=(['case7.pop.h.nday1.' yrten num2str(yyy1+one) '-' mnten num2str(mmm1) '.nc']);
%            afnamex=(['case7.cam2.h0.' yrten num2str(yyy1+one) '-' mnten num2str(mmm1) '.nc']);
%            ofx=netcdf(ofnamex,'nowrite');
%            afx=netcdf(afnamex,'nowrite');
%            sstx2= ofx{'SST'}(:,:,:);
%            sstx2= squeeze(mean(sstx2,1));
%            close(ofx);close(afx);
%            ofx=netcdf(ofnamexx,'nowrite');
%            afx=netcdf(afnamexx,'nowrite');
%            sshxx=ofx{'SSH'}(:,:,:);
%            sshxx=squeeze(mean(sshxx,1));
%            close(ofx);close(afx);
%            sshx2= of{'SSH'}(:,:,:);
%            sshx2= squeeze(mean(sshx2,1));
%            ssh(ssh>10^10)=NaN;
%            sshx2(sshx2>10^10)=NaN; 
            hmxlx(hmxlx>10^10)=NaN; 
%           sstx2(sstx2>10^10)=NaN; 
%           sstx(sstx>10^10)=NaN; 
%           sst(sst>10^10)=NaN;
%           ssh(ssh>10^10)=NaN;
%           sshx(sshx>10^10)=NaN;
%           sshxx(sshxx>10^10)=NaN;
%            pslx2= af{'PSL'}(:,:);  
%            octempx2= of{'TEMP'}(:,nz1:nz2,:,:);           
%            octempx2= squeeze(mean(mean(octempx2(:,:,r-tsng:r+tsng,:),1),3));
%              tsx=  nc_varget(afnamex,'TS'); 
%              pslx= af{'PSL'}(ddd,:,:);
%            octempx(octempx>10^10)=NaN;octemp(octemp>10^10)=NaN;
%           trckmp2(8,ii)= mean2(ts(trckmp(4,i)-tsng:trckmp(4,i)+tsng,trckmp(5,i)-tsng:trckmp(5,i)+tsng))-... 
%                               mean2(tsx(trckmp(4,i)-tsng:trckmp(4,i)+tsng,trckmp(5,i)-tsng:trckmp(5,i)+tsng)); % save the ts area-avg'd around the track to calc corr with psl 
%           trckmp2(7,ii)= (nanmean(nanmean(hmxl(orows,ocols)- hmxlx2(orows,ocols))));
%           sshcorr= cat(1,reshape((ssh(orows,ocols)- sshx2(orows,ocols)),length(orows)*length(ocols),1),sshcorr);
%           pslcorr= cat(1,reshape((psl(arows,acols)- pslx2(arows,acols)),length(arows)*length(acols),1),pslcorr);
%           trckmp2(8,ii)= mean2(sst(r-tsng:r+tsng,c-tsng:c+tsng))-... 
%                                mean2(sstx(r-tsng:r+tsng,c-tsng:c+tsng)); % save the ts area-avg'd around the track to calc corr with psl 
%           trckmp2(10,ii)= psl(trckmp(4,i),trckmp(5,i))- pslx(trckmp(4,i),trckmp(5,i));
%           trckmp2(9,ii)= sign(saveday(sd))*nanmean(nanmean(hmxlx(rows,cols)- hmxl(rows,cols))); 
%           trckmp2(12,ii)= F* cp* trckmp2(7,ii)* trckmp2(9,ii)* dW* dL;  % F=fraction of heat transported from mixed layer to deeper layer, cp=ocean heat capacity, 
%====================================================================
%
%      save the MXL depth anom at TC points to calc difference 
%      btw two days on the same jth track
%
%====================================================================
%           plot(olon(1,c),olat(r,5),'k*');hold on
%           plot(alon(trckmp(5,(trckmp(6,:)==trckmp(6,i)))),alat(trckmp(4,(trckmp(6,:)==trckmp(6,i)))),'*-')
%hold on    
%           mycoastshift(nlon);       
%pause;% clf
%disp('mix layer anom deepend + ')
%           contourf(sign(saveday(sd))*(hmxl(orows,ocols)- hmxlx2(orows,ocols)),50); shading('flat') % hmxl - last month mean
%colorbar; pause; clf
%disp('psl -')
%           contourf(sign(saveday(sd))*psl(arows,acols),50); shading('flat') % hmxl - last month mean
%colorbar; pause; clf
%disp('ssh anom lower -')
%           contourf(sign(saveday(sd))*(ssh(orows,ocols)- sshx2(orows,ocols)),50); shading('flat')
%colorbar; pause; clf
%           contourf(sst(orows,ocols),40);
%           contourf(sign(saveday(sd))*(sstx(orows,ocols)- nanmean(nanmean(sst(orows,ocols)))),40);
%colormap(jet); caxis([-1 1]); colorbar; pause; clf      
%           contourf((sst(orows,ocols)- nanmean(nanmean(sstx(orows,ocols)))));
%colormap(jet); caxis([-1 1]); colorbar; pause; clf      
%           contourf((sst(orows,ocols)- ((sstx2(orows,ocols)))));
%colormap(jet); caxis([-1 1]); colorbar; pause; clf      
%           contourf(lhflx(arows,acols),40)
%hold on;
%           quiver(squeeze(u(end,arows,acols)),squeeze(v(end,arows,acols)))
%colormap(mycmap(0)); shading('flat'); colorbar; pause; clf      
%           plot(mean(ssh(orows,ocols),1) - mean(sshxx(orows,ocols),1))
%           contourf(sshx(orows,ocols))%-sshxx(orows,ocols))


%           contourf(olon(ocols),olat(orows),hmxlx(orows,ocols),30) 
%t=colorbar('peer',gca); set(get(t,'ylabel'),'String', 'cm'); colormap(mycmap(1)); caxis([1000 5000]); shading('flat'); hold on;
%           ux(end,arows(3),acols(end-1))=-15;vx(end,arows(3),acols(end-1))=0;
%           quiver(alon(acols),alat(arows),squeeze(ux(end,arows,acols)),squeeze(vx(end,arows,acols)))
%xlim([olon(ocols(1)) olon(ocols(end))]); ylim([olat(orows(1)) olat(orows(end))]);
%pause; clf

%colormap(jet); colorbar; pause; clf      
%           contourf(olon(1,c-tsng:c+tsng),z_t(nz1:nz2),...
%                    squeeze(octemp(:,r,c-tsng:c+tsng))...
%                   - squeeze(nanmean(nanmean(octempx(:,:,r-ng:r+ng,c-tsng:c+tsng),1),3)),40); shading('flat'); % anomoly against days ago
%                    - octempx2(:,c-tsng:c+tsng),40); shading('flat');  % anomaly against last month avg
% colorbar; pause; clf
%           plot(squeeze(hmxl(trckmp(3,i),r,c(1)-5:c(1)+5))); pause      
%           plot(olon(1,c(1)),olat(r,5),'k*');  pause

%           contourf(olon(r,c(1))-ng:olon(r,c(1))+ng,olat(r,c(1))-ng:olat(r,c(1))+ng,  ...
%                squeeze(hmxl(trckmp(3,i),r-ng:r+ng,c(1)-ng:c(1)+ng)),35,'linestyle','none');

%           contourf(alon,alat,squeeze(ts(:,:)),20,'linestyle','none');
%colorbar;%caxis([205 325])
%===============================================================
%           Save the 1st lagging day and subtract each day by it
%===============================================================
            if sd==1; hmxlx1=hmxlx(r,ocols); end 
%            hmxl_lon_timeseries1(sd,:)=hmxlx(r,ocols)-hmxlx1;
%            hmxl_lon_timeseries2(sd,:)=hmxlx(r,ocols)-hmxlx2(r,ocols);
%            hmxl_lon_timeseries3(ii,:)=sign(saveday(sd))* (hmxlx(r,ocols)- hmxl(r,ocols));
%         end
%====================================
%   Lead Day - Lag Day (1st saveday)
%====================================
    tracknum=trckmp(6,i)
%    contourf(ocols,saveday,hmxl_lon_timeseries1,30);
%caxis([-2500 2500]);shading('flat');colorbar;colormap(mycmap(2));pause;clf
%================================================
%   Lead Day - Lag Day (3rd saveday, for a track)
%================================================
%if trckmp(6,i+1)~=trckmp(6,i)
%    contourf(hmxl_lon_timeseries3,30);
%caxis([-2500 2500]);shading('flat');colorbar;colormap(mycmap(2));pause;clf
%===============================
%   Lead Day - Lag Monthly Avg
%===============================
%    contourf(ocols,saveday,hmxl_lon_timeseries2,30);
%caxis([-1000 1000]);shading('flat');colorbar;pause;clf
    close(af);
    close(of);
    end
 end 
'Done Saving fields in track, (y,m,d,lat,lon,trck#,dhmxl,dssh)'
%corr(trckmp2(9,find(isnan(trckmp2(9,:))==0))',trckmp2(8,find(isnan(trckmp2(8,:))==0))')

%==================================================
%
%          Save trckmp2 (the tracks that cross ocn
%
%==================================================
  save('TCtrckmp.dat','trckmp2','-ascii')
  trckrowsize=size(trckmp2,1);
  trckcolsize=size(trckmp2,2);
  trckmpStr={'trckmp2'}
  nc = netcdf('TCtrckmp.nc','clobber');% 'clobber'); 
  nc.description ='The TC tracks that cross ocn only, no land';
  nc.source='CCSM4 daily output';
  nc.resolution =([num2str(alon(2)-alon(1)) ' lon ' num2str(alat(2)-alat(1))  ' lat']);
  nc('rowsize')= rowsize; %define dim
  nc('colsize')= colsize;
  for n=1:length(trckmpStr)
    nc{trckmpStr{n}}={'rowsize','colsize'}  
    nc{trckmpStr{n}}=trckmpStr{n} %TC that cross ocn
    nc{trckmpStr{n}}(:,:)=trckmp2%eval([trckmpStr{n}]) %TC that cross ocn
  end


end
%================================================================================================================
%
%
%                   Calc Corr between Variables 
%
%
%================================================================================================================
%  scl=mean(nonzeros(trckmp(6,:)))/mean(nonzeros(trckmp(7,:)));  clf %scale to plot 
%  corr(nonzeros(trckmp(6,:)),nonzeros(trckmp(7,:)))
%  plot(nonzeros(trckmp(6,:)));  hold on
%  plot(nonzeros(trckmp(7,:))*scl)
 
  


