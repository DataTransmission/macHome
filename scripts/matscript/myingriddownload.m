ntyear = 12


sstyearlymean=loaddap('http://iridl.ldeo.columbia.edu/expert/SOURCES/.CAC/.sst/T/12/boxAverage/dods')
sst3monthmean=loaddap('http://iridl.ldeo.columbia.edu/expert/SOURCES/.CAC/.sst/T/3/boxAverage/dods')

for 1:ntyear
sstyearlymean( - sst3monthmean =
%Time
%grid: /T (months since 1960-01-01) ordered (Jan-Mar 1970) to (Jan-Mar 2003) by 3. N= 133 pts :grid
%Longitude
%grid: /X (degree_east) ordered (124E) to (70W) by 2. N= 84 pts :grid
%Latitude
%grid: /Y (degree_north) ordered (29S) to (29N) by 2. N= 30 pts :grid
%mtoolbox
