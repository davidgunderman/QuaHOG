function outcurves= correctOrder(incurves)
len=size(incurves,1)/3;
for i=1:(len-1)
    pt=incurves((3*(i-1)+1):(3*i),end);
    dummy=[nan , nan, nan]';
    dummyend=[nan , nan, nan]';
    ctr=i;
    while (~all(abs(pt-dummy)<1e-14) && ~all(abs(pt-dummyend)<1e-14)) && ctr<(len)
        ctr=ctr+1;
        dummy=incurves((3*(ctr-1)+1):(3*ctr),1);
        dummyend=incurves((3*(ctr-1)+1):(3*ctr),end);
    end
    if all(abs(pt-dummy)<1e-14)
        swppt=incurves((3*(ctr-1)+1):(3*ctr),:);
    elseif all(abs(pt-dummyend)<1e-14)
        swppt=fliplr(incurves((3*(ctr-1)+1):(3*ctr),:));
    else
        error('curves are not closed');
    end
        incurves((3*(ctr-1)+1):(3*ctr),:)=incurves((3*(i)+1):(3*(i+1)),:);
        incurves((3*(i)+1):(3*(i+1)),:)=swppt;
end
outcurves=incurves;
end
        