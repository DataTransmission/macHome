function unplot= sem(xplot,lambda,N,basisflag)
%
% Returns the errors in the solution of
% u_xx - \lambda u = 0
% u(-1) = 0
% u_x(1)= 1
%

%N = 2;            % The truncation order
%basisflag=4;
Np= N+1;          % number of degrees of freedom

switch basisflag
% case 1
%%%%%%%%%%%%Basis functions are Chebyshev polynomials%%%%%%%%%%%
%  b = cheb_poly(N+2);          % the basis function
%  xGL =-cos(pi*[0:1/N:1])';  % Gauss-Lobatto collocation points
  case 2
%%%%%%%%%%%%Basis functions are Legendre polynomials%%%%%%%%%%%
   bL = leg_poly(N);          % the basis function
   b(1:N) = bL(1:N)+bL(2:Np); % basis recombination to enforce Dirichlet BC at xi=-1
%  xGL = GLpoints(Np);  % Gauss-Lobatto collocation points
%  [d in] = sort(xGL);
%  xGL = xGL(in(1:Np)); % sort them in ascending order
% case 3
%  xGL =-cos(pi*[0:1/N:1])';  % Gauss-Lobatto collocation points
%  b =card(xGL);               % Gauss-Lobatto collocation points
  case 4
   xGL = GLpoints(Np);  % Gauss-Lobatto collocation points
   [d in] = sort(xGL);
   xGL = xGL(in(1:Np)); % sort them in ascending order
   b =card(xGL);  % Gauss-Lobatto collocation points
end
db = diff(b);     % 1st derivative of basis function


syms x real;

% Build the equations for the PDE
% Apply the operator to the basis functions
A = int(db'*db+lambda*b'*b,x,-1,1);

% Take care of the right hand side
r = subs(b,x,1)'* 1;

% Apply essential BC on right
if (basisflag>2)
  A(1,2:Np) = 0;  % eliminate one equation
  A(2:Np,1) = 0;  % eliminate one unknown
end

uh=A\r;                % Solve the system

un=b*uh;               % numerical solution

unplot=subs(un,x,xplot);
