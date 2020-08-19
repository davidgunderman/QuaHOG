function val = dC_eval(C,t0)
% Given a bernstein polynomial, evaluates it i at t0;
% Input parameters
% C: This is a bernstein polynomial's control points (aka coefficients)
% t0: This is the point at which we wish to split the polynomial
% Output parameters
% val: row vector of values of the bernstein polynomial defined by C at t0
[d,n] = size(C);
val=zeros(length(t0),d);
% if length(t0)>1
%     1
% end
for j=1:d
%     tempdivmat=repmat(C(j,:),1,length(t0));
    cdiv_mat=C(j,:);
    for i=2:n
        cdiv_mat=cdiv_mat.*(1-t0') + circshift(cdiv_mat,-1,2).*t0';
    end
%     cdiv_mat=triu(fliplr(cdiv_mat));

    val(:,j) = cdiv_mat(:,1);
end
