function A = vander_poly(van_degree,v,w)
%VANDER Vandermonde matrix.
%   A = VANDER(V) returns the Vandermonde matrix whose columns
%   are powers of the vector V, that is A(i,j) = v(i)^(n-j).

%   Copyright 1984-2001 The MathWorks, Inc. 
%   $Revision: 5.8 $  $Date: 2001/04/15 12:02:43 $

n = length(v);
v = v(:);
w=w(:);

A=w.*ones(n,1);
A = repmat(A,1,van_degree+1);

for j = van_degree:-1:1
    A(:,j) = v.*A(:,j+1);
end

A=(fliplr(A))';
