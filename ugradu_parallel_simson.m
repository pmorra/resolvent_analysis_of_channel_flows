function [fx,fy,fz,varargout] = ugradu_parallel_simson(field_dir,times,varargin)

rem_ref = false;
do_DFT = false;
get_freq = false; my_alf = []; my_bet = [];

if nargin > 3 && nargin < 6
  printf('\n Wrong number of inputs. STOP.\');
  return
elseif nargin >= 6
  if ~isempty(varargin{2}) && ~isempty(varargin{3}) && ~isempty(varargin{4})
    rem_ref = true;
    u_ref = varargin{2}; v_ref = varargin{3}; w_ref = varargin{4};
    NNy = size(u_ref,2);
    physref = cat(3,permute(u_ref(:,NNy:-1:1,:),[1 3 2]),...
                    permute(v_ref(:,NNy:-1:1,:),[1 3 2]),...
                    permute(w_ref(:,NNy:-1:1,:),[1 3 2]));
    fouref = phys2fou(physref);
  else
    rem_ref = false; fouref = 0;
  end
end

if nargin >= 7 && ~varargin{5} || nargin < 7
  do_DFT = false;
elseif nargin >= 7 && varargin{5}
  do_DFT = true;
  get_freq = false;
  if nargin >= 9
    get_freq = true;
    my_alf = varargin{6};
    my_bet = varargin{7};
  end
end

nt = length(times); my_num = nargin;

% Initialize output
fx = cell(nt,1); fy = fx; fz = fx;
if nargout > 3;
  ux = cell(nt,1); uy = ux; uz = ux;
  n_out = nargout;
else
  ux = []; uy = []; uz = [];
  n_out = 3;
end

if nargin >= 3 && ~isempty(varargin{1});
  nproc = varargin{1}; parpool(nproc);
else
  parpool;
end

parfor it = 1:nt
%for it = 1:1
  if mod(floor(it/nt*100),20) < 1e-2; disp([num2str(it/nt*100),'%']); end
  % filename
  filename = [field_dir,'/field.',num2str(times(it)),'.u'];
  
  % Read DNS data
  [vel0,xF,yF,zF,Lx,Ly,Lz,maxt,Re,~,dstar] = readdns(filename,false);
  
  % Retrieve NNx NNy NNz
  [~,NNx,NNy,NNz]=fou2phys(vel0,0,0);
  
  % Spectral indeces along streamwise direction
  kxvec = linspace(0,2*pi/Lx*(NNx-1),NNx);
 
  % Spectral indeces along spanwise direction
  kzvec = linspace(0,2*pi/Lz*((NNz)/2-1),(NNz)/2);
  kzvec = [kzvec -fliplr(kzvec(2:end))];
  
  % Remove reference value (if given)
  if rem_ref;  vel = vel0 - fouref;  else;  vel = vel0;  end
  
  % Compute the gradients
  [dxfou,dyfou,dzfou] = gradfield(vel,yF,kxvec,kzvec);
  
% 3/2 rule for physical space
  Nx_extra = 0;%0.5*NNx;
  Nz_extra = 0;%0.5*NNz;
    
  % Move to physical space
  vel   = fou2phys(  vel,Nx_extra,Nz_extra);
  dxvel = fou2phys(dxfou,Nx_extra,Nz_extra);
  dyvel = fou2phys(dyfou,Nx_extra,Nz_extra);
  dzvel = fou2phys(dzfou,Nx_extra,Nz_extra);
  
  % Decouple the velocity matrix into the three u,v,w components
  [vel,x,y,z] = mat2vel(vel,xF,yF,zF); 
  dxvel = mat2vel(dxvel,xF,yF,zF);     
  dyvel = mat2vel(dyvel,xF,yF,zF);     
  dzvel = mat2vel(dzvel,xF,yF,zF);     
  
  % Compute the gradient
  if ~do_DFT
    fx{it} = vel.u.*dxvel.u + vel.v.*dyvel.u + vel.w.*dzvel.u; 
    fy{it} = vel.u.*dxvel.v + vel.v.*dyvel.v + vel.w.*dzvel.v;
    fz{it} = vel.u.*dxvel.w + vel.v.*dyvel.w + vel.w.*dzvel.w;
    if n_out > 3
      ux{it} = vel.u;
      uy{it} = vel.v;
      uz{it} = vel.w;
    end
  elseif do_DFT
    % x
    dum = vel.u.*dxvel.u + vel.v.*dyvel.u + vel.w.*dzvel.u;
    if get_freq
      dum_f = fft( fft( dum,[],1 ),[],3);
      fx{it} = dum_f(my_alf,:,my_bet);
    else
      fx{it} = fft( fft( dum,[],1 ),[],3);
    end
    % y
    dum = vel.u.*dxvel.v + vel.v.*dyvel.v + vel.w.*dzvel.v; 
    if get_freq
      dum_f = fft( fft( dum,[],1 ),[],3);
      fy{it} = dum_f(my_alf,:,my_bet);
    else
      fy{it} = fft( fft( dum,[],1 ),[],3);
    end
    % z
    dum = vel.u.*dxvel.w + vel.v.*dyvel.w + vel.w.*dzvel.w;
    if get_freq
      dum_f = fft( fft( dum,[],1 ),[],3);
      fz{it} = dum_f(my_alf,:,my_bet);
    else
      fz{it} = fft( fft( dum,[],1 ),[],3);
    end
    % Output fields as well if asked
    if n_out > 3
      % x
      dum = vel.u;
      if get_freq
        %dum_f = zeros(size(dum,1),size(dum,2),size(dum,3));
        dum_f = fft( fft( dum,[],1 ),[],3);
        ux{it} = dum_f(my_alf,:,my_bet);
      else
        ux{it} = fft( fft( dum,[],1 ),[],3);
      end
      % y
      dum = vel.v; 
      if get_freq
        dum_f = fft( fft( dum,[],1 ),[],3);
        uy{it} = dum_f(my_alf,:,my_bet);
      else
        uy{it} = fft( fft( dum,[],1 ),[],3);
      end
      % z
      dum = vel.w;
      if get_freq
        dum_f = fft( fft( dum,[],1 ),[],3);
        uz{it} = dum_f(my_alf,:,my_bet);
      else
        uz{it} = fft( fft( dum,[],1 ),[],3);
      end
    end
  end
end
varargout{1} = ux; varargout{2} = uy; varargout{3} = uz;
delete(gcp('nocreate'));
