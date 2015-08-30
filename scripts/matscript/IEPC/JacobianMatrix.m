function dJ = JacobianODE(t,J)
%  Steps to find Jacobian Matrix time evolution
%    1. Differentiate Equation of Motion by all state variable x_j and get Stability Matrix A_ij
%    2. Numerically integrate d/dt(J_ij) = A_ij * J_ij  , set J^0_ij = 1
global PR
%==============================================
%           Put the Perturbed PR here
alpha=PR;
%
%==============================================
sigma=10; b=8/3; r=28; k1=-10; k2=-10; tau=0.1;
a1= alpha; a2=alpha/tau;
A = [
     -sigma      sigma   0        -a1                0           0
     r*(1-X(3))  -1      -X(1)    0                  a1          0
     X(2)        X(1)    -b       0                  0           a1
     -tau*a2     0       0        -tau*sigma         tau*sigma   0
     0           tau*a2  0        tau*(r-tau*X(6))   -tau*tau    -tau*tau*X(4)
     0           0       -tau*a2  tau*X(5)           tau*X(4)    -tau*tau*b      ];
A = reshape(A',36,1)
dJ = zeros(36,1);
J = zeros(36,1);
ii=1;
for i=1:6
    for j=1:6
        dJ(ii) = A(ii)*J(6*(j-1)+i);
        ii=ii+1;
    end
end
