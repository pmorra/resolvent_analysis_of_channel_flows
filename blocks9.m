function IP = blocks9(ny)
% indeces to extract pieces of Pff
ii = 1:ny; U = ones(ny,ny); Z = zeros(3*ny,3*ny);
IP(:,:,1) = Z; IP(      ii,      ii, 1 ) = U;                               
IP(:,:,2) = Z; IP(      ii,   ny+ii, 2 ) = U; IP(   ny+ii,    ii, 2 ) = U;  
IP(:,:,3) = Z; IP(      ii, 2*ny+ii, 3 ) = U; IP( 2*ny+ii,    ii, 3 ) = U; 
IP(:,:,4) = Z; IP(   ny+ii,   ny+ii, 4 ) = U;                               
IP(:,:,5) = Z; IP(   ny+ii, 2*ny+ii, 5 ) = U; IP( 2*ny+ii, ny+ii, 5 ) = U;  
IP(:,:,6) = Z; IP( 2*ny+ii, 2*ny+ii, 6 ) = U; 
end