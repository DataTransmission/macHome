%projected basis v= (1,x,x^2,x^3,x^4,x^5)
format long

func=sin(x); % Analytic function to approximate
syms x
v= [1,x,x^2,x^3,x^4,x^5] % v original basis
e(1)= v(1)/sqrt(int(v(1)*v(1),-pi,pi)) % first assume a new orthornormal bases by assuming normal bases of v(1) = e(1)
x1=-pi;x2=pi;
%=================Gram-Schmidt Procedure to find a orthonormal basis with respect to v===================
for j= 2:length(v)
  inp=0;jj=1;
  while(jj<j)
    inp = int(v(j)*e(jj),x1,x2)*e(jj) + inp; %inner products add up
    jj=jj+1;
  end
  e(j)= v(j)- inp; % orthognal bases without normalization
  e(j)= e(j)/sqrt(int(e(j)*e(j),x1,x2)); % orthonormalized bases e(j) to form a orthonormal basis e= (e1,e2...ek)
end

%================Porject vector (sinx) onto bases to get linear coeffs and form the vector's polynomial form=======
polyn=0;
for j=1:length(v)
  coeff(j)=int(func*e(j),x1,x2); % project func onto each bases e(j) to get the coeff
  polyn=coeff(j)*e(j) + polyn; % linear combination of bases e(j) to form polyn
end
x=[-pi:0.01:pi];
plot(x,eval(polyn))
