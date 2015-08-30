%The original coast starts from -180 to 180. This file splits the longitude
%at given longitude(splitlon) and paste whatever on left to the right.
function mycoastshift(nlong)
load coast
splitlon = nlong % split longitude at 0
index=find(long<splitlon)
longg=long;
latgg=lat;
longg(index)=NaN;
latgg(index)=NaN;
clear index;
index=find(long>splitlon)
lonll=long;
latll=lat;
lonll(index)=NaN;
latll(index)=NaN;
clear index:
lonll = lonll + 360; %paste left half to right, right half to left
long=cat(1,longg,lonll);
lat=cat(1,latgg,latll);
plot(long,lat)

