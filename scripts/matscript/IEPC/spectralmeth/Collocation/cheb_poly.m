function Chebp = cheb_poly(N)
%
% return the Chebyshev polynomial of degree 0:N
% usage is Chebp = cheb_poly(N)
% Compute the Chebyshev polynomials using the recurrence relationship.
%
    
    syms x p_1 p_2 p_3 Chebp
    p_3 = 1;		% Legendre polynomial of degree 0
    Chebp(1) = p_3;
    if (N~=0)
      p_2 = p_3;
      p_3 = x;		% Legendre polynomial of degree 1
      Chebp(2) = p_3;
      for k = 2:N 
        p_1 = p_2;
        p_2 = p_3;
        p_3 = ( 2*x*p_2 - p_1 ); % degree k
        Chebp(k+1) = p_3;
      end
    end
