function [Poly xGL w] = myPoly(polyType,N)

switch polyType
case 'Legendre'
     Poly = leg_poly(N);
     [xGL,w,P]= lglnodes(N);
case 'Chebyshev'
     Poly = cheb_poly(N)
     xGL  =-cos(pi*[0:1/N:1])'; 
     w    = 1./sqrt(1-xGL.^2);
case 'Chebyshev2'
case 'Gegenbauer'
end
