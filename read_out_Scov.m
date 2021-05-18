function data = read_out_Scov(filename)
%
% Reads binary file containing outputs from my_welch_Scov
% (as written from write_out_Scov)
%
% INPUT filename : name of file
%
% OUTPUT data : structured object with the data
%           contains: S,SV,pSV,ny,dt,y,ffq
%

% open file
fid = fopen(filename,'r','ieee-le.l64');

% read parameters and frequency axis
data.ny = fread(fid,1,'int');
data.dt = fread(fid,1,'float64');
data.y =  fread(fid,data.ny,'float64');
nf = fread(fid,1,'int');
data.ffq = fread(fid,nf,'float64');

% initialize containers
data.S = zeros(data.ny*3,data.ny*3,nf);
data.SV = data.S;
data.pSV = data.S;

% read S
dum = fread(fid,(data.ny*3)^2*nf,'float64');
data.S = reshape(dum,[data.ny*3 data.ny*3 nf]);
dum = fread(fid,(data.ny*3)^2*nf,'float64');
data.S = data.S + 1i*reshape(dum,[data.ny*3 data.ny*3 nf]);

% read SV
dum = fread(fid,(data.ny*3)^2*nf,'float64');
data.SV = reshape(dum,[data.ny*3 data.ny*3 nf]);

% read pSV
dum = fread(fid,(data.ny*3)^2*nf,'float64');
data.pSV = reshape(dum,[data.ny*3 data.ny*3 nf]);
dum = fread(fid,(data.ny*3)^2*nf,'float64');
data.pSV = data.pSV + 1i*reshape(dum,[data.ny*3 data.ny*3 nf]);

% close file
fclose(fid);
