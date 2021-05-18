function [H] = resolvent_generalized_Cess(u,alpha,beta,L,omega,Re)

zi = sqrt(-1);

ny = length(u);
N = ny-2;

[A,B,C,D,Q,y]=oss(N,alpha,beta,Re,L,u);

% Resolvent
dum = -C*((zi*omega*eye(2*N)+A)\B);
H = dum; 
% H = zeros(3*ny,3*ny);
% 
% id1 = 1:ny;
% id2 = 1:N;
% 
% % first line-block
% H(id1(2:end-1),     id1(2:end-1)) = dum(id2,     id2);
% H(id1(2:end-1),  ny+id1(2:end-1)) = dum(id2,  N+id2);
% H(id1(2:end-1),2*ny+id1(2:end-1)) = dum(id2,2*N+id2);
% 
% % second line-block
% H(ny+id1(2:end-1),     id1(2:end-1)) = dum(N+id2,    id2);
% H(ny+id1(2:end-1),  ny+id1(2:end-1)) = dum(N+id2,  N+id2);
% H(ny+id1(2:end-1),2*ny+id1(2:end-1)) = dum(N+id2,2*N+id2);
% 
% % third line-block
% H(2*ny+id1(2:end-1),     id1(2:end-1)) = dum(2*N+id2,    id2);
% H(2*ny+id1(2:end-1),  ny+id1(2:end-1)) = dum(2*N+id2,  N+id2);
% H(2*ny+id1(2:end-1),2*ny+id1(2:end-1)) = dum(2*N+id2,2*N+id2);
% 
% iy = 1:ny; iy = iy(end:-1:1);
% ii = [iy,ny+iy,2*ny+iy];
% 
% H = H(:,ii);
% H = H(ii,:);
end


function [A,B,C,D,Q,y]=oss(N,alpha,beta,Re,L,u)
% build the Orr-Sommerfeld/Squire system using Chebyshev collocation
%        inputs:
% N: number of inner points
% alpha, beta: streamwise and spanwise wavenumbers
% Re: Reynolds number, based on deltastar
% L,u: boxheight and base flow profile
%        outputs:
% A,B,C,D: state space operators
% Q, y: inner product matrix and collocation points

%%% parameters
zi = sqrt(-1);

%%%% differentiation matrices and base flow
scale=2/L;
[y,DM] = chebdif(N+2,2); %Chebyshev collocation
DM(:,:,1)=DM(:,:,1)*scale;    
DM(:,:,2)=DM(:,:,2)*scale^2;    
y=(y(2:end-1)+1)/scale; 
up=DM(:,:,1)*u; % differentiate base flow
upp=DM(:,:,2)*u;

%%%% implement homogeneous boundary conditions 
sel=2:N+1;
D1=DM(sel,sel,1);
D2=DM(sel,sel,2);
[y,D4]=cheb4c(N+2); % fourth order differentiation matrix
D4=D4*(2/L)^4;

%%%% laplacian
I=eye(N);Z=zeros(N,N);
k2=alpha^2+beta^2;
delta=(D2-k2*I); % Laplacian
delta2=(D4-2*k2*D2+k2*k2*I); % square of the Laplacian

%%%% Cess function
f = f_cess(Re,u);
fp = DM(:,:,1)*f;
fpp = DM(:,:,2)*f;

%%%% build oss matrix
LOS=-zi*alpha*diag(u(sel))*delta+zi*alpha*diag(upp(sel))+delta2/Re;     
LOS=LOS+1/Re*diag(f(sel))*delta2;        % adding to LOS
LOS=LOS+2/Re*diag(fp(sel))*delta*D1;     % the contribution
LOS=LOS+1/Re*diag(fpp(sel))*(D2+k2*I);   % from Cess
LSQ=-zi*alpha*diag(u(sel))+delta/Re;
LSQ=LSQ+1/Re*diag(f(sel))*delta;         % adding to LSQ
LSQ=LSQ+1/Re*diag(fp(sel))*D1;           % the contribution from Cess
LC= -zi*beta*diag(up(sel));
A=[delta\LOS,Z ; LC,LSQ]; % the dynamic matrix for OSS

%%%% input,  output, and inner product  operators
B=[delta\(-zi*alpha*D1),delta\(-k2*I),delta\(-zi*beta*D1) ; ...
   zi*beta*I,Z,-zi*alpha*I];
C=[zi*alpha*D1/k2,-zi*beta*I/k2 ; ...
   I,zeros(N,N) ; ...
   zi*beta*D1/k2,zi*alpha*I/k2];
D=[Z,I,Z; zi*beta*I,Z,-zi*alpha*I]; 
Q=enermat(N,L,DM,k2);
end

function [x, D4] = cheb4c(N)

%  The function [x, D4] =  cheb4c(N) computes the fourth 
%  derivative matrix on Chebyshev interior points, incorporating 
%  the clamped boundary conditions u(1)=u'(1)=u(-1)=u'(-1)=0.
%
%  Input:
%  N:     N-2 = Order of differentiation matrix.  
%               (The interpolant has degree N+1.)
%
%  Output:
%  x:      Interior Chebyshev points (vector of length N-2)
%  D4:     Fourth derivative matrix  (size (N-2)x(N-2))
%
%  The code implements two strategies for enhanced 
%  accuracy suggested by W. Don and S. Solomonoff in 
%  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
%  The two strategies are (a) the use of trigonometric 
%  identities to avoid the computation of differences 
%  x(k)-x(j) and (b) the use of the "flipping trick"
%  which is necessary since sin t can be computed to high
%  relative precision when t is small whereas sin (pi-t) cannot.
   
%  J.A.C. Weideman, S.C. Reddy 1998.

    I = eye(N-2);                   % Identity matrix.
    L = logical(I);                 % Logical identity matrix.

   n1 = floor(N/2-1);               % n1, n2 are indices used 
   n2 = ceil(N/2-1);                % for the flipping trick.

    k = [1:N-2]';                   % Compute theta vector.
   th = k*pi/(N-1);                 

    x = sin(pi*[N-3:-2:3-N]'/(2*(N-1))); % Compute interior Chebyshev points.

    s = [sin(th(1:n1)); flipud(sin(th(1:n2)))];   % s = sin(theta)
                               
alpha = s.^4;                       % Compute weight function
beta1 = -4*s.^2.*x./alpha;          % and its derivatives.
beta2 =  4*(3*x.^2-1)./alpha;   
beta3 = 24*x./alpha;
beta4 = 24./alpha;
    B = [beta1'; beta2'; beta3'; beta4'];

    T = repmat(th/2,1,N-2);                
   DX = 2*sin(T'+T).*sin(T'-T);     % Trigonometric identity 
   DX = [DX(1:n1,:); -flipud(fliplr(DX(1:n2,:)))];   % Flipping trick. 
DX(L) = ones(N-2,1);                % Put 1's on the main diagonal of DX.

   ss = s.^2.*(-1).^k;              % Compute the matrix with entries
    S = ss(:,ones(1,N-2));          % c(k)/c(j)
    C = S./S';                      

    Z = 1./DX;                      % Z contains entries 1/(x(k)-x(j)).
 Z(L) = zeros(size(x));             % with zeros on the diagonal.

    X = Z';                         % X is same as Z', but with 
 X(L) = [];                         % diagonal entries removed.
    X = reshape(X,N-3,N-2);

    Y = ones(N-3,N-2);              % Initialize Y and D vectors.
    D = eye(N-2);                   % Y contains matrix of cumulative sums,
                                    % D scaled differentiation matrices.
for ell = 1:4
          Y = cumsum([B(ell,:); ell*Y(1:N-3,:).*X]); % Recursion for diagonals
          D = ell*Z.*(C.*repmat(diag(D),1,N-2)-D);   % Off-diagonal
       D(L) = Y(N-2,:);                              % Correct the diagonal
DM(:,:,ell) = D;                                     % Store D in DM
end

   D4 = DM(:,:,4);                  % Extract fourth derivative matrix
end

function Q=enermat(N,L,DM,k2)
% build the energy measure matrix for OSS
% 	inputs: 
% N: number of inner points
% L: domain height
% DM: differentiation matrices
% k2=alpha^+beta^2 
% 	output:
% Q: energy measure matrix for collocation points

%%%% integration weights
n=0:1:N+1;
j=0:1:N+1;
b=ones(1,N+2);
b([1 N+2])=0.5;
c=2*b;
b=b/(N+1);
S=cos(n(3:N+2)'*j*(pi/(N+1)));
IWT=L/2*diag(b.*[(2+(c(3:N+2).*((1+(-1).^n(3:N+2))./(1-n(3:N+2).^2)))*S)]);

%%%% build energy measure matrix
QvT=0.125*(DM(:,:,1)'*IWT*DM(:,:,1)/k2+IWT); % for v
QetaT=0.125* IWT/k2; % for eta
Q=[QvT(2:N+1,2:N+1),zeros(N,N);zeros(N,N),QetaT(2:N+1,2:N+1)];
end

function out = f_cess(Re,u)

k = 0.426;%0.4;
A = 25.4;%16;

[y,D] = chebdif(length(u),1);
dUdy = D(:,:,1)*u;
Ret = sqrt(Re* (abs(dUdy(1)) + abs(dUdy(end)))/2 );
yp = Ret*(abs(y)-1);

dum = (1-y.^2).^2 .* (1+2*y.^2).^2 .* (1-exp(yp/A)).^2;
out = 0.5*sqrt( 1+k.^2.*Ret.^2/9.*dum )-0.5;

end

