addpath('/Users/ginochen/matscript/spectralmeth/Collocation')
%addpath('/Users/ginochen/matscript/spectralmeth/Spectral')
%addpath('/Users/ginochen/matscript/spectralmeth/Weakform')
% each ith row is an eqn eval at a ith quad pt, with N+1 quad pts, gives N-1 eqn + 2 bdry eqn  
% each nth col is an nth order basis func, with N+1 basis funcs
%
% phi_n(z_i) = -h''_n(z_i) + kap * del_ni
% [row col] = [ith qpts, nth polyn]
% [  1             0             ...   0             ...  0         
%    phi_0(z_1)    phi_1(z_1)    ...   phi_n(z_1)    ...  phi_N(z_1)
%    phi_0(z_2)    phi_1(z_2)    ...   phi_n(z_2)    ...  phi_N(z_2)
%    ...           ...           ...   ...           ...  ...       
%    phi_0(z_i)    phi_1(z_i)    ...   phi_n(z_i)    ...  phi_N(z_i)
%    ...           ...           ...   ...           ...  ...       
%    phi_0(z_N-1)  phi_1(z_N-1)  ...   phi_n(z_N-1)  ...  phi_N(z_N-1)
%    0             0             ...   0             ...  1             ]
%
%  for BC in question
%  * [ u0 u1 u2 u3 u4 ... uN-1 uN ]' = [0 0 0 0 0 0 0 0 0 ... 0 1]'
%  for other BC
%  * [ u0 u1 u2 u3 u4 ... uN-1 uN ]' = [0 f(z_1) f(z_2) ... f(z_i) ... f(z_N-1) q]'
clear
format long
syms x
i = 1; 
j = 2;
Nmax = 16;
u = [1, 3];
alpha = [1, 5*10^-2, 1D-2];
nx = 100;
kap = u(i)/alpha(j);
xe = -1:2/(nx-1):1;
C =  (exp(kap.*(xe+1))-1)./(exp(2*kap)-1); 
                       % this is our analytic test func
ul=0;
ur=1;
for N = 3:Nmax         % N is exactly N+1, just for notation convenience
                       % C = (exp(kap.*(xGL+1))-1)./(exp(2*kap)-1); 
  N
  xGL = GLpoints2(N);
  xGL = sort(xGL);
  Ca = (exp(kap.*(xGL+1))-1)./(exp(2*kap)-1); 
  h = card(xGL);
  eqns= - diff(h,2) + kap*diff(h,1); 
	               % LHS-adv-diffus-eqn with nth interp func 
  Ain = zeros([N N]);
  Ri  = zeros([N 1]);
  %for i = 2:N-1          % ith collocation pts
  %  Ain(i,1:N)   = subs(eqns,xGL(i)); % 1-nth order card func
  %end
  for n = 1:N   % subs by column will be faster than rows
     Ain(:,n)= subs(eqns(n),xGL); % 1-nth order card func
  end
  Ri(:)      = Ri - Ain(:,1)*ul;
  Ri(:)      = Ri - Ain(:,N)*ur;
  Ain(1,:)   = 0;        % zero b/c BC u = 0 
  Ain(:,1)   = 0;        % zero b/c BC u = 0 
  Ain(1,1)   = 1;        % the 1 rep del_0=1 at u_0, 
  Ain(N,:)   = 0;        % b/c BC u = 1
  Ain(:,N)   = 0;        % b/c BC u = 1
  Ain(N,N)   = 1;        % u_N = 1
  Ri(1)      = ul;             % BC on the left 
  Ri(N)      = ur;             % BC on the right
  Ci         = Ain\Ri;           % Ci = the SOLution on quad pts
  Ceqn       = h*Ci;           % Ceqn = the discrete equation of SOLution 
  Ceqn_cont  = subs(Ceqn,x,xe); % Ceqn for plotting
  hold off;
  plot(xe,C,'k',xe,Ceqn_cont,'r',xGL,Ci,'r*');
%  scatter(xn,Ceqn_disc,'*b'); % random points to calc rmse
  se = sum((C-Ceqn_cont).^2)        %subs(ICfunc,x,xn(i)))^2 + se;
  rmse(N-2) = sqrt(se/(N-2))
pause
end
hold off;
%%!!!!!!!!!!find the rms error on other points other than the qpts
figure
loglog(3:Nmax,rmse,'k')
figure
semilogy(3:Nmax,rmse,'k')
clear

