%The original coast starts from -180 to 180. This file splits the longitude
%at given longitude(nlong) and paste whatever on left to the right.
%function mycoastshift(nlong)
% split longitude at nlong degrees
%load coast
%iless=find(long<nlong); %find the longitude indices less than nlong
% preserve longitude indices greater than nlong, and remove the rest of the indices
%longg=long;
%latgg=lat;

%long(iless)=long(iless) + 360;

%longg(iless)=[];
%latgg(iless)=[];

%igreater=find(long>nlong);
%lonll=long;
%latll=lat;
%lonll(igreater)=[];
%latll(igreater)=[];
%lonll = lonll + 360; % shift the degrees by 360 
%long=cat(1,longg,lonll);
%lat=cat(1,latgg,latll);
%plot(long,lat)

