function J = TdCjacobian(ctrlpoints,ubar,d)
    pt=TdCReval_dt(ctrlpoints,ubar,d);
    J=reshape(pt(1,1,:).*pt(2,2,:)-pt(1,2,:).*pt(2,1,:),size(ubar,1),1);
end