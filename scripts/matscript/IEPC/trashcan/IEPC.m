  % ODE45: explicit Runge-Kutta (4,5) formula, Dormand-PRince pair  
  run ~/matscript/startup
   
  global PR % define global variable for LeoODE to use
  
%======================================== 
%       Specify the Perturbed PR here 
%======================================== 
  alpha   = 0.8;  
  PR      = alpha;
  PRrange = 0.01; % range of the PR
  ndPR    = 10e4;  % # of Perturbations for the given PR
%======================================== 

%========================================
%       Setup here  
%======================================== 
  N              = 5;   % N-th order Polynomial  
  Ny             = 6;   % # of State Variables
  y0(1:Ny)       = 0;   % IC for all State Var
%  y0J(Ny+1:Ny+36)= 1;   % IC for all 36 Jacobian element
  dt             = 10;  % time interval size to sample 
  ndt            = 20*dt; % # of interpolation, for fixed time interval 
%  nt             = dt; % total integrate time, no sampling 
  nt             = 200; % integrate time, with sampling
  numBins        = 10e2; % numBins+1 = numbers of bin to determine which solution is most likely for different parameters
  %y(1)=xa,y(2)=ya,y(3)=za,y(4)=xo,y(5)=yo,y(6)=zo
%========================================


  [qpt,w,P]= lglnodes(N);
  
  [p1 p2 PRend PRpts]=myPRPert(PR,PRrange,ndPR);

  qpt_shift = p1+ p2*qpt;
  
  Leg = leg_poly(N);  
  
  % Legendre for Projection onto Param Space
  Legendre_qdtr(1:N+1,1)=1; % Leg(qdtr,k) 0th order Legendre_qdtr for all six state var
  Legendre_qdtr(1:N+1,2:N+1) = subs(Leg(2:N+1),qpt_shift); % basis function at the shifted quadrature pts for all six State var
  % Legendre for non-quadrature pts that one wants to know
  Legendre_soln(1:ndPR,1) = 1; % Leg(pts,k)
  Legendre_soln(1:ndPR,2:N+1) = subs(Leg(2:N+1),PRpts); 

  % projection matrix
  for k = 1:N+1
      Proj(:,k) = (Legendre_qdtr(:,k) .* w) / int(Leg(k)^2,-1,1);  %Proj(qpts,k)
  end
  t=0; % time zero

  % ODE integration function setup  
  Tol=1e-4; %Tolerence Level
  options = odeset('RelTol',1e-4,'AbsTol',[Tol Tol Tol Tol Tol Tol]);



%========================================
%
%       CONTROL CASE (Perturbed only at the initial time to get ensemble runs)
%     
%    1. Solve the system of equations for all state variable by plugging in few Initial Perturbations and run it for deterministic solutions for comparison 
%========================================
  [T,Y_ctrl] = ode45(@LeoODE,[t nt],y0,options); %time interval [0 1], initial state yi, error: Y size changes bc of numerical precision
  figure
  for i=1:Ny  
      plot(0:nt/(length(T)-1):nt,Y_ctrl(:,i));hold on; 
  end


%========================================
%
%       INTERACTIVE ENSEMBLE CASE (Perturbed only at the initial time to get ensemble runs, and average at given time intervals for the atmosphere)
%     
%    1. Solve the system of equations for all state variable by plugging in few Initial Perturbations and run it for deterministic solutions for comparison 
%========================================
  

%========================================
%
%       PERTURBED CASE (Perturbed at each given time interval)
%   
%    1. Project initial state variables onto Parameter Space by using orthonormal basis getting the initial coefficients 
%    2. Get the Initial State variables expanded in the Parameter Space by summing coeffs and basis functions
%    3. Solve the system of equations for all state variable by plugging in the Initial Perturbed State variables as initial conditions
%      (Perturbations are drawn equally once since it's uniform distribution)
%========================================  
  ii=1;
  while t< nt
        T_intrp=[t:dt/(ndt-1):t+dt]'; % Since T is always different under Runge-Kutta Scheme, so fix the time interval to T_intrp and intrp the Soln onto these time 
        for nqpts=1:N+1 % # of param perturb
            PR = qpt_shift(nqpts);
            [T,Y] = ode45(@LeoODE,[t t+dt],y0,options); % Y(Var,T) time interval to pert [t t+dt], initial state y, error: Y size changes bc of numerical precision
            Y_intrp(:,nqpts,ii:ii+length(T_intrp)-1) = interp1(T,Y,T_intrp)'; % Y(Var,qpts,T) Intrp onto T_intrp time intervals, Y(T,Var), Y_intrp(T,Var,PR)
        end
        yhat = Y_intrp(:,:,end) * Proj % yhat(Var,k)), Y(Var,qpts,T), Proj(qpts,k) the qpts will give exact integration (quadrature),
                                                 % the k-th coeff needs to multiply the kth basis function to get the solution in a polynomial form;
        y = yhat * Legendre_soln'; % IC for perturbed solution, y(Var,npt), yhat(Var,k), Leg(npt,k)
        for nv = 1:6
            botEdge = min(y(nv,:)); 
            topEdge = max(y(nv,:)); % sets the bin edges according to max and min values of Y_intrp solutions
            bins = linspace(botEdge, topEdge, numBins+1);
            [numValInEachBin whichBin] = histc(y(nv,:),bins);
            y0(nv) = mean(bins(numValInEachBin==max(numValInEachBin)));
        end
        ii= ii+length(T_intrp);
        t=t+dt
  end
% Recover the Jacobian matrix and find the eigenvalues of the J^T*J matrix at each time step to find the maximum growth rate
%  for t=1:15
%      J(i,j) = Y_intrp(t,Ny+ii,1);
%%  end
%svds(1)
  figure
  for t = 1:nt 
      yhat0(:,t) = Y_intrp(:,:,t) * Proj(:,1); % yhat(Var,k)), Y(Var,qpts,T), Proj(qpts,k) the qpts will give exact integration (quadrature),
  end
  for i=1:Ny % Var
      plot(0:nt/(size(yhat0,2)-1):nt,yhat0(i,:));hold on; 
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
  figure;  myplotState(YavgPert,1,2,3,1);
  figure;  myplotState(YavgPert,4,5,6,1);
  figure;  myplotState(Y_ctrl,  1,2,3,1);
  figure;  myplotState(Y_ctrl,  4,5,6,1);


%J0 =   
