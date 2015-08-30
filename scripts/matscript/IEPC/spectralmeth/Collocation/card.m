function lagint = card(G)

[M n] = size(G);

syms x lagint;
for n = 1:M
% ratj is the interpolation polynomial n.
 lagint(n) = 1; %lagrange interp poly
 for m = 1:n-1
  lagint(n) = lagint(n) * (x-G(m))/(G(n)-G(m));
 end
 for m = n+1:M
  lagint(n) = lagint(n) * (x-G(m))/(G(n)-G(m));
 end
end
