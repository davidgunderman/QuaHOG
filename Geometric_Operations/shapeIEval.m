function eval = shapeIEval(pt,integrand,shape,rat,tol)
    global evalCounter;
    eval = zeros(size(pt,1),1);
    for i=1:size(pt,1)
        if isContainedPt(pt(i,:),shape,rat,tol)
            eval(i)=integrand(pt(i,1),pt(i,2));
        else
            eval(i)=0;
        end
    end
    evalCounter=evalCounter+1;
    scatter(pt(:,1),pt(:,2),'k.')
    hold on
end