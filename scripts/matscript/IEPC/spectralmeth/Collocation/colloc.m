function em = colloc(N,basisflag)
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
  case 1
%%%%%%%%%%%%Basis functions are Chebyshev polynomials%%%%%%%%%%%
   b = cheb_poly(N);  % the basis function
   xGL =-cos(pi*[0:1/N:1])';  % Gauss-Lobatto collocation points
  case 2
%%%%%%%%%%%%Basis functions are Legendre polynomials%%%%%%%%%%%
   b  = leg_poly(N);  % the basis function
   xGL = GLpoints(Np);  % Gauss-Lobatto collocation points
   [d in] = sort(eval(xGL));
   xGL = xGL(in(1:Np)); % sort them in ascending order
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
lambda=4;
pde = d2b-lambda*b;

% Build the system of equations
A = zeros([Np Np]);
A(1,:) = subs(b,x,xGL(1));    % left BC
for i = 2:Np-1
  A(i,:) = subs(pde,x,xGL(i));  % enforce PDF at point xGL(i)
end
A(Np,:) = subs(db,x,xGL(Np));    % right BC

% Take care of the right hand side
r = zeros([Np 1]);
r(Np) = 1;

uh=A\r;                % Solve the system

un=b*uh;               % numerical solution

slambda=sqrt(lambda);
ue=sinh(slambda*(x+1))/(slambda*cosh(2*slambda));

Nplot=200;
xplot=-1:2/Nplot:1;
unplot=subs(un,x,xplot);
ueplot=subs(ue,x,xplot);
%switch (basisflag)
%  case 1
%    plot(xplot,unplot,'r',xplot,ueplot,'k')
%  case 2
%    plot(xplot,unplot,'b',xplot,ueplot,'k')
%  case 3
%    plot(xplot,unplot,'m',xplot,ueplot,'k')
%  case 4
%    plot(xplot,unplot,'c',xplot,ueplot,'k')
%end
%plot(xplot,unplot,'r',xplot,ueplot,'k')
%if (N<4) 
%  labeltag=strcat('$',num2str(N),'$');
%  text(xplot(3*Nplot/4), unplot(3*Nplot/4),labeltag,'Interpret','Latex','FontSize',16,'HorizontalAlignment','right','Color','red')
%end
err=unplot-ueplot;
 switch (basisflag)
   case 1
     plot(xplot,err,'r');
   case 2
     plot(xplot,err,'b');
   case 3
     plot(xplot(1:4:Nplot),err(1:4:Nplot),'mx','MarkerSize',12);
   case 4
     plot(xplot(1:4:Nplot),err(1:4:Nplot),'cx','MarkerSize',12);
 end
emax = max(err);
erms = sqrt(sum(err.^2))/Nplot;
em = [N erms emax];
