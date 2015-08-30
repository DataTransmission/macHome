function xGL = GLpoints(m)

%
% Usage is xGL = GLpoints(m)
% Find the m Gauss-Lobatto roots of the Legendre polynomial
%

syms x;
Ls = leg_poly(m-1);
xGL=roots(sym2poly(collect((1-x^2)*diff(Ls(m)),x)));
[d in] = sort(xGL);
xGL = xGL(in(m:-1:1));
