function fout=readgrd(filename)

% 'r':open file for read
fid = fopen(filename,'r') ;
% scan to get the floating numbers from the txt file
ftmp=textscan(fid,'%f');
% tansfer the structure/cell data to matrix form
fout=cell2mat(ftmp);
% close the file
fclose(fid);
