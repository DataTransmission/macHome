alpha=PR;
X(1:6)=0;
%
%==============================================
sigma=10; b=8/3; r=28; k1=-10; k2=-10; tau=0.1;
a1= alpha; a2=alpha/tau;
dX = zeros(6,1);
%dX = zeros(6+36,1);    % a colu
A = [ 
     -sigma      sigma   0        -a1                0           0
     r*(1-X(3))  -1      -X(1)    0                  a1          0
     X(2)        X(1)    -b       0                  0           a1
     -tau*a2     0       0        -tau*sigma         tau*sigma   0
     0           tau*a2  0        tau*(r-tau*X(6))   -tau*tau    -tau*tau*X(4)
     0           0       -tau*a2  tau*X(5)           tau*X(4)    -tau*tau*b      ];
