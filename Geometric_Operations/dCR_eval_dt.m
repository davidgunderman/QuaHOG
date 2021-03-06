function val = dCR_eval_dt(C,t0)
% Given a bernstein polynomial, evaluates it i at t0;
% Input parameters
% C: This is a bernstein polynomial's control points (aka coefficients)
% t0: This is the point at which we wish to split the polynomial
% Output parameters
% val: row vector of values of the bernstein polynomial defined by C at t0
t0=reshape(t0,1,1,length(t0));
[d,n] = size(C);
d=d-1;
val=zeros(length(t0),d);
cwv_mat=zeros(n,n,length(t0));
cwv_mat(1,:,:)=repmat(C(d+1,:),1,1,length(t0));
wv=dC_eval(C(d+1,:),t0);
wdv=dC_eval_dt(C(d+1,:),t0);
nv=dC_eval(C(1:3,:),t0);
ndv=dC_eval_dt(C(1:3,:),t0);
val=(wv.*ndv-wdv.*nv)./(wv.^2);
% for i=2:n
%     cwv_mat(i,:,:)=cwv_mat(i-1,:,:).*(1-t0) + circshift(cwv_mat(i-1,:,:),-1,2).*t0;
% end
% for j=1:d
%     cdiv_mat=zeros(n,n,length(t0));
%     cdiv_mat(1,:,:)=repmat(C(j,:),1,1,length(t0));
%     for i=2:n
%         cdiv_mat(i,:,:)=cdiv_mat(i-1,:,:).*(1-t0) + circshift(cdiv_mat(i-1,:,:),-1,2).*t0;
%     end
% %     cdiv_mat=triu(fliplr(cdiv_mat));
% %     ptmat=cdiv_mat./cwv_mat;
% 
%     val(:,j) = (n-1)*(cdiv_mat(end-1,2,:)./cwv_mat(end-1,2,:)-cdiv_mat(end-1,1,:)./cwv_mat(end-1,1,:));
% end
% wt= (cwv_mat(end-1,2,:).*cwv_mat(end-1,1,:))./(cwv_mat(end,1,:).^2);
% val=(val).*reshape(wt,length(t0),1);
end