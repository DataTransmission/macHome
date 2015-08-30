function xGL = GLpoints(m)

%
% Usage is xGL = GLpoints(m)
% Find the m Gauss-Lobatto roots of the Legendre polynomial
%

syms x;
Ls = leg_poly(m-1);
xGL=solve((1-x^2)*diff(Ls(m)));
[d in] = sort(eval(xGL));
xGL = xGL(in(m:-1:1));
