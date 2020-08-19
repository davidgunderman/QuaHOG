function plotNURBScurve(d,p,Xi,P,w)
    if (size(Xi,2)-p-1~=length(w))
        error("Length of weight vector must be p+1 less than length of knot vector")
    else
        numScatter=10000;
        vals=zeros(numScatter,d);
        for i=1:numScatter
            for di=1:d
                totval=0;
                for j=1:(length(Xi)-p-1)
                    [~,Bvals]=ComputeSplineBasis((i-1)/numScatter,p,(Xi/max(Xi)));
                    vals(i,di)=vals(i,di)+w(j)*P(di,j)*Bvals(j);
                    totval=totval+w(j)*Bvals(j);
                end
                vals(i,di)=vals(i,di)/totval;
            end
        end
        if d==2
            plot(vals(:,1),vals(:,2),'r');
        elseif d==3
            plot3(vals(:,1),vals(:,2),vals(:,3),'r');
        else
            error("Dimension must be 2 or 3")
        end
    end
end