function C_1 = dC_split(C,t0)
% Given a bernstein polynomial, splits it into two polynomials at t0;
% Input parameters
% C: This is a bernstein polynomial's control points (aka coefficients)
% t0: This is the point at which we wish to split the polynomial
% Output parameters
% C_1: Coefficients of resulting polynomials (rows 1:d are first poly, rows
% d+1:2*d are second poly)
[d,n] = size(C);
C_1 = zeros(2*d,n);
% C_2 = zeros(d,n);
for j=1:d
    cdiv_mat=zeros(n);
    cdiv_mat(1,:)=C(j,:);
    for i=2:n
        cdiv_mat(i,:)=cdiv_mat(i-1,:)*(1-t0) + circshift(cdiv_mat(i-1,:),-1)*t0;
    end
    cdiv_mat=triu(fliplr(cdiv_mat));

    C_1(j,:) = cdiv_mat(:,end);
    C_1(j+d,:) = flipud(diag(cdiv_mat));
end
