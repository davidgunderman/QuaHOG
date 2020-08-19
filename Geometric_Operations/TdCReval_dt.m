function pt = TdCReval_dt(ctrlpoints,ubar,d)
    ubar(:,3)=1-ubar(:,1)-ubar(:,2);
    pt=zeros(2,d,size(ubar,1));
    partialsHOM=TdCeval_dt(ctrlpoints,ubar,d+1);
    evalsHOM=TdCeval(ctrlpoints,ubar,d+1);
    for jj=1:2
        for ii=1:size(ubar,1)
            pt(jj,:,ii)=(evalsHOM(ii,d+1).*partialsHOM(jj,1:d,ii)-evalsHOM(ii,1:d).*partialsHOM(jj,d+1,ii))./(evalsHOM(ii,d+1).^2);
        end
    end
end