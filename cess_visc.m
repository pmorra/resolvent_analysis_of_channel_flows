function out = cess_visc(Re,u)

k = 0.426;
A = 25.4;

[y,D] = chebdif(length(u),1);
dUdy = D(:,:,1)*u;
Ret = sqrt(Re* (abs(dUdy(1)) + abs(dUdy(end)))/2 );
yp = Ret*(abs(y)-1);

dum = (1-y.^2).^2 .* (1+2*y.^2).^2 .* (1-exp(yp/A)).^2;
out = 0.5*sqrt( 1+k.^2.*Ret.^2/9.*dum )-0.5;

end