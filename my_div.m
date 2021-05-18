function div = my_div(alpha,beta,N);

% Compute the divergence

zi = sqrt(-1);

[~,D] = chebdif(N,1); %D = -D;

div = [zi*alpha*eye(N,N), D, zi*beta*eye(N,N)];

end