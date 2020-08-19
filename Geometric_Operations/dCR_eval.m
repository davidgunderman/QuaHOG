function vals= dCR_eval(CPs,t)
    tt=dC_eval(CPs,t);
    if size(CPs,1)==3
    vals=tt(:,1:2)./tt(:,3);
    elseif size(CPs,1)==4
    vals=tt(:,1:3)./tt(:,4);
    end
%     if any(t==1)
%         vals(t==1,:)=CPs(1:2,end-1)';
%     end
end