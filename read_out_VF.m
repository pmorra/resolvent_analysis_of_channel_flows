function data = read_out_VF(filename)
%
% Reads binary file containing outputs from my_welch_Scov
% (as written from write_out_Scov)
%
% INPUT filename : name of file
%
% OUTPUT data : structured object with the data
%           contains: V,F,ny,dt,y,nt
%

% open file
fid = fopen(filename,'r','ieee-le.l64');

% read parameters and frequency axis
data.ny = fread(fid,1,'int');
data.dt = fread(fid,1,'float64');
data.y =  fread(fid,data.ny,'float64');
data.nt = fread(fid,1,'int');

% initialize containers
data.V = zeros(data.nt,data.ny*3);
data.F = data.V;

% read V
dum = fread(fid,(data.ny*3)*data.nt,'float64');
data.V = reshape(dum,[data.nt data.ny*3]);
dum = fread(fid,(data.ny*3)*data.nt,'float64');
data.V = data.V + 1i*reshape(dum,[data.nt data.ny*3]);

% read F
dum = fread(fid,(data.ny*3)*data.nt,'float64');
data.F = reshape(dum,[data.nt data.ny*3]);
dum = fread(fid,(data.ny*3)*data.nt,'float64');
data.F = data.F + 1i*reshape(dum,[data.nt data.ny*3]);


% close file
fclose(fid);
