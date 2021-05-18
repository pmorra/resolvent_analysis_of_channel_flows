function H = resolvent_4var(U,omegaref,alf,bet,Re,ny)
%
% Computes the resolvent operator H based on the given parameters
%
% INPUT U : reference profile (used to build linearized system)
%       omegaref : reference omega
%       alf      : reference alpha
%       bet 	 : reference beta
%       Re  	 : Reynolds number
%       ny  	 : number of points along y-direction
%
% OUTPUT H : resolvent operator
%

% operator for derivative in y
[~,D] = chebdif(ny,1); %D = -D;

% operator for double derivative in y
D2= D*D;

% dU/dy
dU = D*U;

% useful building blocks
II = eye(size(D));
ZZ = zeros(size(D));

% equivalent frequency k
k = sqrt(alf^2+bet^2);

% phase speed c
c = omegaref/alf;

% direct operator L = (-i*om*I -A) 
L = [1i*alf*diag(U-c)-1/Re*(D2-k^2*II)                           diag(dU)                                 ZZ  1i*alf*II;
                                    ZZ  1i*alf*diag(U-c)-1/Re*(D2-k^2*II)                                 ZZ          D;
                                    ZZ                                 ZZ  1i*alf*diag(U-c)-1/Re*(D2-k^2*II)  1i*bet*II;
                             1i*alf*II                                  D                          1i*bet*II         ZZ];

% number of points in y (should be np == ny)
np = length(U);
if np ~= ny;
  error('length(U) not equal to ny!')
  return
end


% u = 0 (BC)
L(1,:)  = 0;  L(1,1)   = 1; 
L(np,:) = 0;  L(np,np) = 1;

% w = 0 (BC)
L(2*np+1,:) = 0;  L(2*np+1,2*np+1) = 1; 
L(3*np,:)   = 0;  L(3*np,3*np)     = 1;

% v = 0 (BC)
L(np+1,:) = 0; L(np+1,np+1) = 1; 
L(2*np,:) = 0; L(2*np,2*np) = 1;

% compute resolvent
R = inv(L);
B = [ II ZZ ZZ;
      ZZ II ZZ;
      ZZ ZZ II;
      ZZ ZZ ZZ];
H = R*B;





