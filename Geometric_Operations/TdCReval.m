function pt = TdCReval(ctrlpoints,ubar,d)
    tt=TdCeval(ctrlpoints,ubar,d+1);
    pt=tt(:,1:(end-1))./tt(:,end);
end

            