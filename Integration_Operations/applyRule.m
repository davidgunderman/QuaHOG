function I = applyRule(xq,yq,wq,funct)
    I=wq'*funct(xq,yq);
end