  % 'LE' 'PC'
  % ODE45: explicit Runge-Kutta (4,5) formula, Dormand-PRince pair  
  run ~/matscript/startup
   
  global PR Ny % define global variable for LeoODE to use
  
%======================================== 
%       Specify the Perturbed PR here 
%======================================== 
  N       = 5;  % N-th order Polynomial  
  alpha   = 0.8;  
  PR      = alpha;
  PRrange = 0.1; % range of the PR
  ndPR    = 1e4; % # of Perturbations for the given PR
  pts     = [-1:2/(ndPR-1):1]';
  [p1 p2 PRend PRpts] = myPRPert(PR,PRrange,ndPR);
  [Poly qpts w]       = myPoly('Legendre',N);
  qpt_shift           = p1 + p2*qpts;
    % Legendre for Projection onto Param Space
    % Leg(qdtr,k) 0th order Legendre_qdtr for all six state var
  Poly_qdtr(1:N+1,1)=1; 
  Poly_qdtr(1:N+1,2:N+1) = subs(Poly(2:N+1),qpts); 
    % basis function at the shifted quadrature pts for all six State var
    % Legendre for non-quadrature pts that one wants to know
  Poly_sample(1:ndPR,1) = 1; % Leg(pts,k)
  Poly_sample(1:ndPR,2:N+1) = subs(Poly(2:N+1),pts); 

%========================================
%       Parameters
%======================================== 
  Ny               = 6;   % # of State Variables
  y0(1:Ny)         = 0;   % IC for all State Var
  y0(Ny+1:Ny+Ny^2) = 1; % IC for all 36 Jacobian element
  InitialTime      = 0;
  TStep            = 0.01;
  TStepNum         = 10;
  IterationTime    = TStep*TStepNum; % Time to update and normalize Jacobian Matrix
  nt_eq            = 200  % initial time for PC to start, this purpose is to run
                        % until equilibrium and start the integration at this stage
  dt               = 100;  % time interval size to sample 
  ndt              = 20*dt; % # of interpolation, for fixed time interval 
%  nt             = dt; % total integrate time, no sampling 
  nt               = 250; % integrate time, with sampling
  numBins          = 9*10e3; % numBins+1 = numbers of bin to determine 
                           % which solution is most likely for different parameters
             %y(1)=xa,y(2)=ya,y(3)=za,y(4)=xo,y(5)=yo,y(6)=zo

%========================================
% ODE integration function setup  
%========================================
  RelTol=1e-5;
  AbsTol=1e-6; %Tolerence Level
  Tolvec=ones(1,Ny+Ny^2)*AbsTol;
  options = odeset('RelTol',RelTol,'AbsTol',Tolvec);
%========================================

type = input('LE or PC');

switch type

case 'LE'
%=====================================================================
%
%            Lyapunov Exponent
%   (the time average of Lyapunov number summed)
%
%=====================================================================
%  Recover the Jacobian matrix and find the eigenvalues 
%  of the J^T*J matrix at each time step to find the maximum growth rate
%  InitialTime = timestep;
  for nqpts = 1:N+1
      Sum = 0;
      itr=0;
      T1 = InitialTime;
      T2 = T1 + IterationTime;
      TSpan = [T1:TStep:T2];
      while (T2 < nt)
          PR = qpt_shift(nqpts); 
          itr=itr+1; 
              % itrth iteration updating the time averaged Lyapunov Exponent
          T=itr*IterationTime; % the same as T2 if T1=0, else 
          [t,Y] = ode45(@LeoODE,TSpan,y0,options); 
              % time interval [0 1], initial state yi, 
              % error: Y size changes bc of numerical precision
          ii=Ny+1;
          for j = 1:Ny
              for i = 1:Ny
                  J(i,j) = Y(end,ii);
                  ii=ii+1;
              end
          end
          [Q R] = qr(J);
          permission=1;
          for i=1:Ny
             if R(i,i)==0
                permission=0;
                break;
             end
          end
          if (permission)
             Sum = Sum + log(abs(diag(R)));
             lam = Sum / T;  
                 % Lyapunov exponent = time avg of summed Lyapunov number 
             Lambda(:,itr,nqpts) = flipud(sort(lam));  
                 % sort lam to get the maximum Lambda at each iteration
          end
          y0 = [Y(end,1:Ny)';Q(:)];  
                 % replace Jacobian with the orthonormal vectors Q(:) 
                 % to normalize Jacobian matrix at t2
          T1 = T1 + IterationTime;
          T2 = T2 + IterationTime;
          TSpan = [T1:TStep:T2];
      end
      figure
      plot(lam(nqpts,:));  hold on;
  end



%====================================================================
%
%       INTERACTIVE ENSEMBLE CASE 
%    (Perturbed only at the initial time to get ensemble runs, 
%      and average at given time intervals for the atmosphere)
%     
%    1. Solve the system of equations for all state variable 
%       by plugging in few Initial Perturbations and run it 
%       for deterministic solutions for comparison 
%
%====================================================================  

case 'PC'
%====================================================================
%
%       CONTROL CASE 
%     (Perturbed only at the initial time to get ensemble runs)
%     
%    1. Solve the system of equations for all state variable by 
%       plugging in few Initial Perturbations and run it for 
%       deterministic solutions for comparison 
%
%====================================================================
  [T_ctrl,Y_ctrl] = ode45(@LeoODE,[InitialTime InitialTime+nt_eq],y0,options); 
  figure;
  for i=1:Ny  
      plot(0:nt/(length(T_ctrl)-1):nt,Y_ctrl(:,i));hold on; 
  end

%====================================================================
%
%       PERTURBED CASE (Perturbed at each given time interval)
%    1. Put the perturbed coefficient (on the quadrature points) 
%       into the equations and run it for a given time interval
%    2. Project final state variables onto Parameter Space by 
%       using orthonormal basis getting the coefficients
%    3. Approximate the solution by series expansion of coeff and basis
%    4. Sample the solution 10e4 times with coeff between [-1 1], 
%       find the maximum frequency solution to restart all the 
%       state variables, repeat 1 2 3 
%    5. Calculate the mean by 0th order Coeff, variance by summing up 
%       the square of (1st to nth order coeff) and (basis norm) 
%
%====================================================================  


  % projection matrix
  for k = 1:N+1
      Proj(:,k) = (Poly_qdtr(:,k) .* w(:)) / (myNormSqr(k,'Legendre'));  
              %int(Leg(k)^2,-1,1)); %Proj(qpts,k), Leg(qpts,k), w(qpts)
  end
  ii= 1;
  t = InitialTime+nt_eq;
  if isempty(Y_ctrl)==0 % if not empty
     y0 = Y_ctrl(end,:);
  end
  while t< nt
        T_intrp=[t:dt/(ndt-1):t+dt]'; 
              % Since T is always different under Runge-Kutta Scheme,so
              % fix the time interval to T_intrp and intrp the Soln onto these time 
        for nqpts=1:N+1 % # of param perturb
            PR = qpt_shift(nqpts);
            [T,Y] = ode45(@LeoODE,[t t+dt],y0,options); 
              % Y(T,Var) time interval to pert [t t+dt], initial state y, 
              % error: Y size changes bc of numerical precision
            Y_intrp(:,nqpts,ii:ii+length(T_intrp)-1) = interp1(T,Y,T_intrp)'; 
              % Y(Var,qpts,T) Intrp onto T_intrp time intervals, 
              % Y(T,Var), Y_intrp(T,Var,PR)
            nqpts
        end
        %===============================================================
        %
        %  Take one parm sample and plot the time series
        %   ( if the estimated soln not exact then the 
        %     time series should have jumps )
        %
        %===============================================================
        for tt = 1:length(T_intrp)
            yhat = Y_intrp(1,:,tt) * Proj; 
              % take the atm Var for example, the ocn var are slowly varing
              % yhat(1,k)
            y_1pt(tt) = yhat * Poly_sample(99,:)'; % take any sample pt
        end
        figure; plot(y_1pt)
        %================================================================
        %
        %   Calc |y_realization - y_estimate| = dy on the qpts 
        %    ( this shows how much time needed for the next restart
        %      in order to avoid projected solution not exact )
        %
        %================================================================
        ttt=0; % # of tt count 
        for tt = 2:length(T_intrp)/10:length(T_intrp)
            ttt=ttt+1;
            yhat = Y_intrp(1:Ny,:,tt) * Proj(:,:);
              % yhat(Var,k)), Y(Var,qpts,T), Proj(qpts,k) 
              % the qpts will give exact integration (quadrature),
              % the k-th coeff needs to multiply the kth basis function 
              % to get the solution in a polynomial form;
            y_qpt = yhat * Poly_qdtr'; % estimated solution on the qpts
              % y(Var,qpt)
            dy = y_qpt - Y_intrp(1:Ny,:,tt); 
              % dy(Var,qpt,ttt)
            dy_norm(:,ttt) = sqrt(sum(dy.^2,2)) ./ sqrt(sum(y_qpt.^2,2)) 
              % dy_norm(Var,ttt)
              % diff btw int'd and estimated solution on all qpts
            logyhat(:,:,ttt) = log(abs(yhat'));
            figure; plot(logyhat(:,:,ttt));
              % plot logyhat, expect yhat decrease monotonically with k
        end
        figure; plot(T_intrp(2:length(T_intrp)/10:length(T_intrp)),log10(dy_norm)); 
        xlabel('Nondimensional time'); 
        ylabel('$ \log_{10}( \sqrt{\sum_q{dy^2}} / \sqrt{\sum_q{y^2}} ) $','interpreter','latex')     
        title('Time Series of Relative Difference between Realization and Estimated Solution summed over qpts 5th order')
              % difference between integrated solution and projected solution
              % ydiff becomes smaller when decrease integration time dt, 
              % or when increase order of polynomial
'done'
pause
	      % expect yhat decreases with k ( yhat approaches a monotonically 
              % decreasing line if integration time dt decreases, if dt is 50
              % then yhat doesn't drop monotonically even when taken 10th order,
              % it remains a flat line )
        y_pt = yhat(1:Ny,:) * Poly_sample'; 
              % y(Var,npt), yhat(Var,k), Leg(npt,k), IC for perturbed solution
        figure; plot(pts,y_pt(6,:)); 
        % expect the expanded solution converges as k-order and time increases
        xlabel('Parameters');  
        ylabel('Single State Variable Magnitude');  
        title('Single State Variable Expanded in Parameter Space at 1st Sampling time')
        for nv = 1:Ny
            botEdge = min(y_pt(nv,:)); 
            topEdge = max(y_pt(nv,:)); 
              % sets the bin edges according to max 
              % and min values of Y_intrp solutions
            bins = linspace(botEdge, topEdge, numBins);
            [numValInEachBin whichBin] = histc(y_pt(nv,:)',bins');
            y0(nv) = mean(bins(find(numValInEachBin==max(numValInEachBin))));
        end
        ii= ii+length(T_intrp);
        t=t+dt+1
  end


%===================================================
%
% Plot the mean of PC expansion with time
%
%===================================================
  figure
%  aa = Proj * Legendre_sample'
  for t = 1:size(Y_intrp,3) 
%      y_sample(:,:,t) = Y_intrp(:,:,t) * aa; 
         %Proj * Legendre_sample'; 
         %y_sample(Var,npt,T),Proj(qpts,k),Legendre_sample(npt,k)
      y_mean(:,t) =   Y_intrp(:,:,t) * Proj(:,1); 
         % yhat(Var,k)), Y_intrp(Var,qpts,T), Proj(qpts,k) 
         % the qpts will give exact integration (quadrature),
      y_variance(:,t) = (Y_intrp(:,:,t) * Proj(:,2:N+1)).^2 * myNormSqr(2:N+1,1)'.^2 ;
         % y_variance(Var,T), Y_intrp(Var,qpts,T), Proj(qpts,k), norm(1,k) 
  end
  for i=6%1:Ny % Var
      plot(0:nt/(size(Y_intrp,3)-1):nt,y_mean(i,:));hold on; 
      plot(0:nt/(size(Y_intrp,3)-1):nt,(y_mean(i,:) +...
           2*sqrt(y_variance(i,:)/(N+1-1+1-N))),'r'); 
         % numbers of quadrature points detemines the sample size (N+1) 
      plot(0:nt/(size(Y_intrp,3)-1):nt,(y_mean(i,:) -... 
           2*sqrt(y_variance(i,:)/(N+1-1+1-N))),'r'); % 2*stdev
  end
'done'
pause
  
%==================================================
%
%       State Space Plots
%
%==================================================
  for nqpts =1:N+1  
      subplot(N+1,1,nqpts);  myplotState(Y_intrp(:,:,nqpts),1,2,3,1);
  end
  figure;  myplotState(y_mean',1,2,3,1);
  figure;  myplotState(y_mean',4,5,6,1);
  figure;  myplotState(Y_ctrl,  1,2,3,1);
  figure;  myplotState(Y_ctrl,  4,5,6,1);

end
