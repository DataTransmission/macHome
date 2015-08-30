function Legp = leg_poly(N)
% return the Legendre polynomial of degree 0:N
% usage is Legp = Legendre_Polynomials(N)
% Compute the Legendre polynomials using the recurrence relationship.
%
    
    syms x p_1 p_2 p_3 Legp
    p_3 = 1;		% Legendre polynomial of degree 0
    Legp(1) = p_3;
    if (N~=0)
      p_2 = p_3;
      p_3 = x;		% Legendre polynomial of degree 1
      Legp(2) = p_3;
      for k = 2:N 
        p_1 = p_2;
        p_2 = p_3;
        p_3 = ( (2*k-1)*x*p_2 - (k-1)*p_1 ) / k; % degree k
        Legp(k+1) = p_3;
      end
    end
