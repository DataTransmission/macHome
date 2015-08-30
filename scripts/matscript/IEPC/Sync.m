% psi:  analytic soln
% s:    arbitrary signal
% shat: Hilbert Transform of s
% A:    amplitude
% phi:  phase

% psi(t)  = s(t) + j*shat(t) = A(t)*exp(j*phi(t))
% shat(t) = pi^-1 P.V. int(s(tau)/(t-tau),tau,-inf,inf)
% P.V. means integral taken in the sense of the Cauchy principal value

% shat is the convolution of s and 1/(pi*t)
% Fourier transform Shat of shat = product of fourier transform of s and 1/(pi*t)
