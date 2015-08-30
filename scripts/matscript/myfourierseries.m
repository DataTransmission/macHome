syms x;
f=sin(x)+x^3; %test function
nmin=-35
nmax=35
n(1:abs(nmin))=nmin:-1;
n(abs(nmin)+1:2*nmax)=1:nmax;
xmin=-pi;
xmax=pi;
%basis=sin(n.*pi*x*2/(xmax-xmin)); % odd fourier series
%basis=cos(n.*x/2) % even fourier series
basis=exp(i.*n.*pi*x*2/(xmax-xmin)); % odd + even fourier series
%basis=cos(n.*x/2)-i*sin(n.*x/2);
%xmin=0;
a0=1/(xmax-xmin)*int(f,x,xmin,xmax);
a=1/(xmax-xmin)*int(f.*basis,x,xmin,xmax);

syms x1;
%ibasis=sin(n.*pi*x1*2/(xmax-xmin));
%ibasis=cos(n.*x1/2);
ibasis=exp(-i.*n.*pi*x1*2/(xmax-xmin));
%ibasis=cos(n.*x1/2) + i*sin(n.*x1/2);
u=a0+sum(a.*ibasis);
n1=10;
xrange=xmin:(xmax-xmin)/n1:xmax;
u = subs(u,x1,xrange);

plot(xrange,real(u)); %the imaginery part is all 0, so no need to plot
hold on;
plot(xrange,subs(f,x,xrange),'r');
hold off;
