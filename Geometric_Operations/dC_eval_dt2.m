function val = dC_eval_dt2(C,t0)
% Given a bernstein polynomial, evaluates its second derivative at t0;
% Input parameters
% C: This is a bernstein polynomial's control points (aka coefficients)
% t0: This is the point at which we wish to split the polynomial
% Output parameters
% val: row vector of values of the bernstein polynomial defined by C at t0
t0=reshape(t0,1,1,length(t0));
[d,n] = size(C);
val=zeros(length(t0),d);
for j=1:d
    cdiv_mat=zeros(n,n,length(t0));
    cdiv_mat(1,:,:)=repmat(C(j,:),1,1,length(t0));
    for i=2:n
        cdiv_mat(i,:,:)=cdiv_mat(i-1,:,:).*(1-t0) + circshift(cdiv_mat(i-1,:,:),-1,2).*t0;
    end
%     cdiv_mat=triu(fliplr(cdiv_mat));
    if n>2
%     val(:,j) = (n-1)*(cdiv_mat(end-1,2,:)-cdiv_mat(end-1,1,:));
        val(:,j) = (n-1)*(n-2)*(cdiv_mat(end-2,3,:)-2*cdiv_mat(end-2,2,:)+cdiv_mat(end-2,1,:));
    else
        val=0;
    end
end
