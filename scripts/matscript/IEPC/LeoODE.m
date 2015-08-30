function OUT = LeoODE(t,X)
global PR Ny
%==============================================
%           Put the Perturbed PR here
alpha=PR;
%
%==============================================
sigma=10; b=8/3; r=28; k1=-10; k2=-10; tau=0.1;
a1= alpha; a2=alpha/tau; 
dX = zeros(Ny,1);
%dX = zeros(6+36,1);    % a column vector
dX(1) = sigma*(X(2)-X(1)) - a1*(X(4)+k1);
dX(2) = r*X(1) - X(2) - X(1)*X(3) + a1*(X(5)+k1);
dX(3) = X(1)*X(2) - b*X(3) + a1*X(6);
dX(4) = tau*(sigma*(X(5)-X(4)) - a2*(X(1)+k2));
dX(5) = tau*(r*X(4) - tau*X(5) - tau*X(4)*X(6) + a2*(X(2)+k2));
dX(6) = tau*(X(4)*X(5) - tau*b*X(6) - a2*X(3));


%=============================================
%
%      6*6 Jacobian terms X(7:7+35) ODE
%
%=============================================
%     X(1)        X(2)    X(3)     X(4)               X(5)        X(6)
A = [ 
     -sigma      sigma   0        -a1                0           0
     r*(1-X(3))  -1      -X(1)    0                  a1          0
     X(2)        X(1)    -b       0                  0           a1
     -tau*a2     0       0        -tau*sigma         tau*sigma   0
     0           tau*a2  0        tau*(r-tau*X(6))   -tau*tau    -tau*tau*X(4)
     0           0       -tau*a2  tau*X(5)           tau*X(4)    -tau*tau*b      ];     

J = reshape(X(Ny+1:Ny+Ny^2),Ny,Ny);
dJ = A*J;      
OUT = [dX(:);dJ(:)];
%dX(Ny+1:Ny+Ny^2) = A*B;
% this is how Jacobian Matrix Should look like

%         X(1)   X(2)       X(6)

%  J = [ dX(7)  dX(13) ...
%        dX(8)  dX(14) ...
%        dX(9)  dX(15) ...
%        dX(10) dX(16) ...
%        dX(11) dX(17) ...
%        dX(12) dX(18) ... dX(36)]

