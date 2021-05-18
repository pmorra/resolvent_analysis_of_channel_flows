function write_out_VF(filename,V,F,ny,dt,y,nt)
%
% Write binary file containing the outputs from my_welch_Scov
%
% INPUT filename: name of file
%       V   : V matrix (nt,nx)
%       F   : F matrix (nt,nx)
%       ny  : number of points along y
%       dt  : sampling delta t
%       y   : y axis 
%

% open file
fid = fopen(filename,'w','ieee-le.l64');

% write parameters and frequency axis
fwrite(fid,ny,'int');
fwrite(fid,dt,'float64');
fwrite(fid,y,'float64');
fwrite(fid,nt,'int');

% write V
dum = reshape(real(V),[length(V(:)) 1]);
fwrite(fid,dum,'float64');
dum = reshape(imag(V),[length(V(:)) 1]);
fwrite(fid,dum,'float64');

% write F
dum = reshape(real(F),[length(F(:)) 1]);
fwrite(fid,dum,'float64');
dum = reshape(imag(F),[length(F(:)) 1]);
fwrite(fid,dum,'float64');

% close file
fclose(fid);

end
