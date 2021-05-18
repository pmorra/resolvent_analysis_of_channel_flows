function write_out_Scov(filename,S,SV,pSV,ny,dt,y,ffq)
%
% Write binary file containing the outputs from my_welch_Scov
%
% INPUT filename: name of file
%       S   : S matrix (output from my_welch_Scov)
%       SV  : SV matrix (output from my_welch_Scov)
%       pSV : pSV matrix (output from my_welch_Scov)
%       ny  : number of points along y
%       dt  : sampling delta t
%       y   : y axis 
%       ffq : frequency axis
%

% open file
fid = fopen(filename,'w','ieee-le.l64');

% write parameters and frequency axis
fwrite(fid,ny,'int');
fwrite(fid,dt,'float64');
fwrite(fid,y,'float64');
fwrite(fid,length(ffq),'int');
fwrite(fid,ffq,'float64');

% write S
dum = reshape(real(S),[length(S(:)) 1]);
fwrite(fid,dum,'float64');
dum = reshape(imag(S),[length(S(:)) 1]);
fwrite(fid,dum,'float64');

% write SV
dum = reshape(SV,[length(SV(:)) 1]);
fwrite(fid,dum,'float64');

% write pSV
dum = reshape(real(pSV),[length(pSV(:)) 1]);
fwrite(fid,dum,'float64');
dum = reshape(imag(pSV),[length(pSV(:)) 1]);
fwrite(fid,dum,'float64');

% close file
fclose(fid);

end
