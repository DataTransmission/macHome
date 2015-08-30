  % data vector structure:
  %   row: 1~10 = 10 yrs of Jan; 11~20 = 10 yrs of Feb and so on.
  %   col: spatial row vector of 432 grid points for the selected Eq Pacific zone
  % Eigenvector structure:
  %   each data vector has m=432 dimensions, so each eigenvector should be m=432 long
  %   and a total of m=432 # of eigenvector.
  % mth principle component: the projection of the data vetor x' onto the mth eigenvector e_m
  %   a scalar quantity with 10 data for each mode, representing 10 yrs of projection
%================Parameter========================
  lat = nc_varget('gT42NoOMP.cam.1333-10.nc','lat');
  lon = nc_varget('gT42NoOMP.cam.1333-10.nc','lon');
  nyr = 200;
  inyr=1/nyr;
  dt = 1; %1yr
  nmon = 12;
  mm=1;
  nmode=25, inmode=1/nmode;
  nlat= 5;
%  var = {'ssta'};
  var = {'ssta' 'tauxa' 'tauya'};
  iv = 1; clim= 1.5;%1 is sst, 2 is taux, 3 is tauy
  sm = 1; 
  testing= input('1. Build M (transition matrix) 2.Test M and plot  (1 or 2)?   ')
  if testing==1
    yy=800, build= 1, test= 0 %build M but don't test
  else
    yy=400, build= 0, test= 1 %test M and plot it 
  end
  isave= 0; %1 if eof, 2 if meof
  if nyr <= nmode
    input('nmode cannot exceed nyr')
    break
  end
%=================================================
%====================Read in variables=============
  for iy = 1:nyr
    for im=1:nmon
      if im<10
        zero=0;
      else
        zero=[];
      end
      a = nc_varget(['gT42NoOMP.cam.0' num2str(yy) '-' num2str(zero) num2str(mm) '.nc'],'TS');
      a1= nc_varget(['gT42NoOMP.cam.0' num2str(yy) '-' num2str(zero) num2str(mm) '.nc'],'TAUX');
      a2= nc_varget(['gT42NoOMP.cam.0' num2str(yy) '-' num2str(zero) num2str(mm) '.nc'],'TAUY');
        % restrict to Eq Pacific
      a = a(find(lat>-nlat & lat<nlat),find(lon>=120.0 & lon<=270.0));
      a1= a1(find(lat>-nlat & lat<nlat),find(lon>=120.0 & lon<=270.0));
      a2= a2(find(lat>-nlat & lat<nlat),find(lon>=120.0 & lon<=270.0));
      ii=1;
      for i = 1:size(a,1)
         sst(iy+(im-1)*nyr,ii:(ii-1)+size(a,2))  = a(i,:); % make a row vector
        taux(iy+(im-1)*nyr,ii:(ii-1)+size(a1,2)) = -a1(i,:); % make a row vector
        tauy(iy+(im-1)*nyr,ii:(ii-1)+size(a2,2)) = -a2(i,:); % make a row vector
          % row(nyr*nmon); col(nlat*nlon)
          % put all Jan in 1st nyr rows and all Feb in next nyr rows and so on.
        ii=ii+size(a,2);
      end
    mm=mm+1;
    end
  yy=yy+1;
  mm=1;
  end

%===========================Do EOF twice========================
  ii=1;
  for im = 1:nmon
    for iy = 1:nyr
      ssta(iy+(im-1)*nyr,:)  = sst(iy+(im-1)*nyr,:)  - sum( sst(ii:ii-1+nyr,:))*inyr;
      tauxa(iy+(im-1)*nyr,:) = taux(iy+(im-1)*nyr,:) - sum(taux(ii:ii-1+nyr,:))*inyr;
      tauya(iy+(im-1)*nyr,:) = tauy(iy+(im-1)*nyr,:) - sum(tauy(ii:ii-1+nyr,:))*inyr;
      varmat(:,:,iy,im) = reshape(eval([var{iv} '(iy+(im-1)*nyr,:)']),size(a,2),size(a,1))';% anomaly matrix with 0 mean (for plotting)
    end 
    jj=1;
    for ivar = 1:length(var)
%=======Calc EOF once using state variable independently=================
%      [EV, EOFs, PC, error, norms] = EOF(double(eval([var{ivar} '(ii:ii-1+nyr,:)'])),nmode,[]);
       [ C, lambda, EOFs ] = svds(double(eval([var{ivar} '(ii:ii-1+nyr,:)'])), nmode);
       PC= C* lambda;
       loads{1}=PC;
       loads{2}=EOFs;
       [sgns, loads]=sign_flip(loads,double(eval([var{ivar} '(ii:ii-1+nyr,:)'])));
       PC=loads{1};
       EOFs=loads{2};
       clear loads
        %PC*EOFs'-sst(ii:ii-1+nyr,:); %confirm if PC*Eigenvectors consist the variable in space time
        %L(1)/sum(L) %total of L eigenvalue, how much variance the 1 eigen mode explains
       stdev(im,ivar) = sqrt(sum(sum(PC'*PC*inyr)));
       EOFs1_monthly(:,:,im,ivar) = EOFs;
       PC1(:,jj:jj-1+nmode) = PC/stdev(im,ivar); %Principle component: projection of variable onto eigenvector directions.   
       jj=jj+nmode;
    end
%=======Calc EOF once again using combined PC from each state variable===
%    [EV, EOFs, PC, error, norms] = EOF(double(PC1),nmode,[]);
    [ C, lambda, EOFs ] = svds(double(PC1), nmode);
    PC= C* lambda; % Field = C * Lambda * EOFs', 
                   % Field * EOFs = C * Lambda = PC
                   % PC coeff: Project the fields onto Eigenvectors and see how 
                   %           big the response were, each Eigenvector has
                   %           it's own PC time series, to reflect how the 
                   %           field behaves in that eigen direction with time
    loads{1}=PC;
    loads{2}=EOFs;
    [sgns, loads]=sign_flip(loads,double(PC1));
    PC=loads{1};
    EOFs=loads{2};
    clear loads
    EOFs_monthly(:,:,im)= EOFs;
    PC_monthly(:,:,im) = PC; 
    if im>1 & build == 1  
      dd  = PC_monthly(1,:,im)'* PC_monthly(1,:,im); % cov matrix of PC= vector*vectortranspose
      dd1 = PC_monthly(1,:,im)'* PC_monthly(1,:,im-1);
      for iy = 2:nyr
        dd = PC_monthly(iy,:,im)'* PC_monthly(iy,:,im) + dd; %autocov matrix, add all yr PC and avg
        dd1 = PC_monthly(iy,:,im)'* PC_monthly(iy,:,im-1) + dd1; %lag-1 cov matrix
      end
%=======Make Seasonal Transition Matrix/Markov Model from autocov and lag-1 cov matrix===
      M(:,:,im-1)= dd1* (inyr* inmode)* inv(dd* (inyr* inmode)); 
    end
    ii=ii+nyr;
  end
  if build ==1
    dd  = PC_monthly(1,:,1)'* PC_monthly(1,:,1);
    dd1 = PC_monthly(1,:,1)'* PC_monthly(1,:,12);
    for iy=2:nyr
      dd  = PC_monthly(iy,:,1)'* PC_monthly(iy,:,1) + dd; %autocov matrix, add all yr PC and avg
      dd1 = PC_monthly(iy,:,1)'* PC_monthly(iy,:,12) + dd1; %lag-1 cov matrix
    end
    M(:,:,12) =  dd1* (inyr* inmode)* inv(dd* (inyr* inmode)); 
  end

%===================Test Markov Model==================================================
%======Test Markov Model, if nm=6, means use 6 seasonal M to estimate the 7th month state variable===
  if test==1
    for sm=1:12 %sm starting month
      testPC= M(:,:,sm)* PC_monthly(1,:,sm)'; %first test month predicting 2nd test month for all yrs
      ii=1+sm*nyr; %ssta start from predicting 3rd month
      jj=1;
      bb=sm; %restart after passing the 12th month of 1st yr of prediction
      yrcount=1;
      for nm=1:46 %nm month lead
            % test the (nm+1)th month ssta using M (transition matrix) 6 times
        if nm>1
          testPC= M(:,:,bb+jj-1)* testPC;
        end
        if rem(bb+jj,13)==0
          jj=0;bb=1; 
        end
        testvar= testPC'* EOFs_monthly(1:nmode,:,bb+jj)'* EOFs1_monthly(:,:,bb+jj,iv)'* stdev(bb+jj,iv);
        if ii==nmon*nyr+yrcount
          ii=rem(ii,nmon*nyr)+1;
          yrcount= yrcount+1;
        end
        corr_testvar_var(nm,sm)=corr(testvar',eval([var{iv} '(ii,:)'])'); %spacial avg correlation
%        testvarmat= reshape(testvar,size(a,2),size(a,1))';
%        figure;subplot(2,1,1);contourf(lon(find(lon>=120.0 & lon<=270.0)),lat(find(lat>-nlat & lat<nlat)),varmat(:,:,yrcount,bb+jj),20,'linestyle','none');
%         colorbar;colormap(b2r);caxis([-clim clim]);
%        subplot(2,1,2);contourf(lon(find(lon>=120.0 & lon<=270.0)),lat(find(lat>-nlat & lat<nlat)),testvarmat,20,'linestyle','none');colorbar;caxis([-clim clim]);
        jj=jj+1;
        ii=ii+nyr;
        if isave==1
          print(gcf,['eof_tauxa_' num2str(nm) 'monthlead'],'-djpeg','-r250')
        elseif isave==2
          print(gcf,['meof_tauxa_' num2str(nm) 'monthlead'],'-djpeg','-r250')
        end
      end
    end
  end
%contourf(abs(corr_testvar_var'),20,'linestyle','none');colorbar;colormap(flipud(bone))
contourf(1:46,1:12,(corr_testvar_var'),20,'linestyle','none');colorbar;colormap(flipud(b2r))
caxis([-0.9 0.9])
%title('Correlation between CCSM3 and Predicted SSTa - SSTa as predictor')
title('Correlation between CCSM3 and Predicted SSTa - multiple predictor')
xlabel('leading month');ylabel('starting month')
%print(gcf,'corr_ccsm3_ssta_sstaonly','-djpeg','-r250')
print(gcf,'corr_ccsm3_ssta_multiplepredictor_sgnfixed','-djpeg','-r250')
