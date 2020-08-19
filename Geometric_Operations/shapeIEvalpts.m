function pts = shapeIEvalpts(pt,shape,rat,tol)
    eval = zeros(size(pt,1),1);
    for i=1:size(pt,1)
        eval(i)= isContainedPt(pt(i,1:2),shape,rat,tol);
    end
    pts=pt(logical(eval),:);
end