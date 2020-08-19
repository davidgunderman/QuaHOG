function avgWeight= avgWeight(curves)
    l=0; w=0;
    for i=1:length(curves)
        weights= curves{i}(3:3:end,:);
        l=l+length(weights(:));
        w=w+sum(weights(:));
    end
    avgWeight=w/l;
end