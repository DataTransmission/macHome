function map = cmap(casenumber,rows)
%BLUE2RED  colormap showing a gradient from blue to white to red.

if nargin == 1 % no rows input
   rows = 157;  %this should remain odd # for the [1 1 1] to be at center
end
irow = linspace(1,rows,5); % the number of boundary rows
nrow = numel(irow(1):irow(2)); % the number of rows between two boundaries
%levels = linspace(0,100,rows)'; %level for Johna

if nargin ==0
   casenumber = 1;
end
% 0.49 1    0.83  (am) aquamarine 
% 0.3  0.8  0     (av) avocado 
% 0    0.8  0.8   (aq) aqua/turquoise 
% 0    0    1     (b)  blue 
% 0.25 0.25 0.9   (bc) blue cobalt
% 0    0    0.8   (bd) blue dark 
% 0    0    0     (bk) black 
% 0.5  0.5  1     (bl) blue light
% 0.8  0.5  0     (br) brown 
% 0    1    1     (c)  cyan 
% 1    0.62 0.4   (co) copper 
% 1    1    0.6   (cr) cream 
% 0.5  0    0     (rd) red dark
% 1    0.5  0     (yd) yellow dark 
% 0    1    0     (g)  green 
% 0.9  0.75 0     (gl) gold 
% 0.5  0.5  0.5   (gy) grey 
% 1    0    1     (m)  magenta 
% 1    0.6  0     (or) orange 
% 1    0    0     (r)  red 
% 1    0.5  0.5   (p)  peach
% 1    1    0     (y)  yellow
% 1    1    1     (w)  white 
%
switch casenumber
  case 1 % b c w y r
    rvec = [0 0 1 1 1];
    gvec = [0 1 1 1 0];
    bvec = [1 1 1 0 0];
  case 2 % w c g y r
    rvec = [1 0 0 1 1];
    gvec = [1 1 1 1 0];
    bvec = [1 1 0 0 0];
  case 3 % b m w y r
    rvec = [0 1 1 1 1];
    gvec = [0 0 1 1 0];
    bvec = [1 1 1 0 0]; 
  case 4 % m b w y r
    rvec = [1 0 1 1 1];
    gvec = [0 0 1 1 0];
    bvec = [1 1 1 0 0];
  case 5 % w y dy dr ddr 
    rvec = [1 1 1.0 0.5 0.2];
    gvec = [1 1 0.5 0.0 0.0];
    bvec = [1 0 0.0 0.0 0.0];
  case 6 % w k
    rvec = [1 0.75 0.5 0.25 0];
    gvec = [1 0.75 0.5 0.25 0];
    bvec = [1 0.75 0.5 0.25 0];
  case 7 % c w y r dr
    rvec = [0 1 1 1 0.5];
    gvec = [1 1 1 0 0.0];
    bvec = [1 1 0 0 0.0];
  case 8 % c y r r dr
    rvec = [0 1 1 1 0.5];
    gvec = [1 1 0 0 0.0];
    bvec = [1 0 0 0 0.0];
  case 9 % lb y or r rd
    rvec = [0 1 1.0 1 0.5];
    gvec = [1 1 0.6 0 0.0];
    bvec = [1 0 0.0 0 0.0];
end

    for i = 1:numel(irow)-1
       r(irow(i):irow(i+1),1)  = linspace(rvec(i),rvec(i+1),nrow);
       g(irow(i):irow(i+1),1)  = linspace(gvec(i),gvec(i+1),nrow);
       b(irow(i):irow(i+1),1)  = linspace(bvec(i),bvec(i+1),nrow);
    end

map = [r g b];
%map_level = [levels 1e2*r 1e2*g 1e2*b] % output the corresponding row index in the first column for Johna

