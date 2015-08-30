
a = [ 15 150 300 ];
b = [ 1  2   3 ];
c = [ 2 3 3 ];

for i =1:1000
  a(i) - b(i)
end


r = regress(a,b);
