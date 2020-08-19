function Inside = isContainedPt(pt, Element2,rat,tol)
% Returns true if pt is inside Element2
    d=2;
    if rat
        Element2W=Element2((d+1):(d+1):end,:);
        Element2C=Element2;
        Element2C((d+1):(d+1):end,:)=[];
        for i=1:size(Element2W)
            Element2C(((i-1)*d+1):(i*d),:)=Element2C(((i-1)*d+1):(i*d),:)./Element2W(i,:);
        end
    end
    numc2=size(Element2,1)/(d+rat);
    n2=size(Element2,2);
    if rat
        firstcps=reshape(Element2C(:,1),2,numc2)';
        outside=[pt(1),min(firstcps(:,2)-rand(1))];
        LineGuess=[linspace(pt(1),outside(1),n2); linspace(pt(2),outside(2),n2); ones(1,n2);];
    else
        firstcps=reshape(Element2(:,1),2,numc2)';
        outside=[pt(1),min(firstcps(:,2)-rand(1))];
        LineGuess=[linspace(pt(1),outside(1),n2); linspace(pt(2),outside(2),n2)];
    end
    counter=0;
for j=1:numc2
    temps=zeros(0,0);
    if rat
        if min(Element2C((2*(j-1)+1),:))<pt(1) && max(Element2C((2*(j-1)+1),:))>pt(1)
           [~,~,temps,~]=Ratbern_inter_recurse_v2(LineGuess,Element2((3*(j-1)+1):(3*j),:),[0 1], [0 1],tol);
        end
    else
        if min(Element2((((d+1)*(j-1))+1),:))<pt(1) && max(Element2((((d+1)*(j-1))+1),:))>pt(1)
           [~,~,temps,~]=bern_inter_recurse_v2(LineGuess,Element2((d*(j-1)+1):(d*j),:),[0 1], [0 1], tol);
        end
    end
    counter=counter+size(temps,1);
end
Inside=logical(mod(counter,2));
end