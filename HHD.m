function [divfree,curlfree,gradp,residual] = HHD(f,alpha,beta,ny,sc)
% differential by chebyshev collacation weights
[~,DM] = chebdif(ny,2);
D1 = sc*DM(:,:,1);
D2 = DM(:,:,2);

% build operators: divergence, gradient, laplacian
div = [1i*alpha*eye(ny,ny), D1, 1i*beta*eye(ny,ny)];
k2 = alpha^2 +beta^2;
grad = [1i*alpha*eye(ny,ny); D1; 1i*beta*eye(ny,ny)];
lap = D2-k2*eye(ny,ny);

% A*p = b : system for Poisson equation
b = div*f;
A = lap; 

% boundary conditions
A(1,:)   = D1(1,:);
A(end,:) = D1(end,:); 
b(1)   =  f(ny+1,1);
b(end) = -f(2*ny,1);

% solve: p is in the Helmotz decomposition w = u + grad(p)
p = A\b;

% extract the divergence free part of the Helmotz decomposition
fdf = f - grad*p;

divfree = fdf;
curlfree = p;
gradp = grad*p;
residual = div*fdf;
end