clear
format long;
k = [8,16,24,32,48,64,80,100,200,400,800,10000];
%===========Analytically derived spectral coeff=======
ks =1:1:k(length(k));
syms x n;
%u=sin(3*x);
u = (x^2-pi^2)^2;
uc = int(u*exp(i*n*x),x,0,2*pi)/(2*pi); % uc is the coefficient at different n
%uc = -24*(-1).^ks./(ks.^4); % derived by hand, should be same as above
auc=abs(subs(uc,n,ks)); % subs n with ks and get uc from n=8 to n=1000
hold off;
loglog(ks,auc,'r') % uc should become smaller with bigger n, which represents the error getting smaller
%plot(log(ks),auc,'r')
hold on;



%===========FFT (spectral coeff)======================
% fft will generate symmetric array, 
% the 1st is a_0, the 1st half is c_n by +n, the 2nd half is conjugate(c_n) by -n 
for i=1:length(k)
  x = 0:2*pi/k(i):2*pi-2*pi/k(i);
%  u=sin(3*x);
 u=(x.^2-pi^2).^2;
  fftu=fft(u,k(i))/k(i);  %/k(i) to normalize
  kd = 1:1:k(i)/2; %starting from 1:8 and 1:16 and 1:32 ... to 1:1000
%  ifftu = ifft(fftu1,2*k(i)+1)
  loglog(kd,abs(fftu(2:k(i)/2+1)))
%plot(log(kd),fftu)
pause
end
%plot(ifftu)
%==========Plot=========================

%subplot(1,2,2)
%plot(rmsu,'g');
%ylabel('uc')
%xlabel('N')
%set(gca,'xticklabel',{'8' '16' '24' '32' '48' '64'});
%legend('continuous case spectral coeff')
%
%title('The rms of Spectral Coefficients')
%
%figure;
%plot(logfftu,'b');hold on
%plot(loguc,'g');
%ylabel('uc')
%xlabel('N')
%set(gca,'xticklabel',{'8' '16' '24' '32' '48' '64'});
%itle('The log10(rms of Spectral Coefficients)')
