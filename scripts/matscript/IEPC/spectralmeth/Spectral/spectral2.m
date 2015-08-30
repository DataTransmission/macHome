function unplot= spectral(xplot,lambda,N,basisflag)
%
% Returns the errors in the solution of
% u_xx - \lambda u = 0
% u(-1) = 0
% u_x(1)= 1
%

%N = 2;            % The truncation order
%basisflag=4;
Np= N+1;          % number of degrees of freedom

K = spectd2(N,basisflag);     % matrix of second derivative
switch basisflag
  case 1
%%%%%%%%%%%%Basis functions are Chebyshev polynomials%%%%%%%%%%%
   b = cheb_poly(N);          % the basis function
   bnorm = pi*ones(size(b));
   bnorm(1)=2*bnorm(1);
  case 2
%%%%%%%%%%%%Basis functions are Legendre polynomials%%%%%%%%%%%
   b  = leg_poly(N);  % the basis function
   bnorm = 2./(1:2:2*N+1);
  case 3
   xGL =-cos(pi*[0:1/N:1])';  % Gauss-Lobatto collocation points
   b =card(xGL);               % Gauss-Lobatto collocation points
  case 4
   xGL = GLpoints(Np);  % Gauss-Lobatto collocation points
   [d in] = sort(eval(xGL));
   xGL = xGL(in(1:Np)); % sort them in ascending order
   b =card(xGL);  % Gauss-Lobatto collocation points
end
db = diff(b);     % 1st derivative of basis function
d2b = diff(db);   % 2nd derivative of basis function



syms x real;

% Build the equations for the PDE
% Apply the operator to the basis functions
pde = d2b-lambda*b;

% Build the system of equations
A = zeros([Np Np]);
A(1,:) = subs(b,x,-1);    % left BC
if (basisflag<3)
%        Modal basis
 for i = 2:Np-1
   A(i,:) = bnorm(i-1)*K(i-1,:); % Galerkin enforcement
   A(i,i-1) = A(i,i-1)-lambda*bnorm(i-1); % Galerkin enforcement
 end
else
%        Nodal basis: skip Dirichlet BC function
 for i = 2:Np-1
   A(i,:) = int( b(i)*pde,x,-1,1); % Galerkin enforcement
 end
end
A(Np,:) = subs(db,x,1);    % right BC

% Take care of the right hand side
r = zeros([Np 1]);
r(Np) = 1;

uh=A\r;                % Solve the system

un=b*uh;               % numerical solution

unplot=subs(un,x,xplot);
