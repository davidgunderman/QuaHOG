function KnotNum= which_knot(xi,Xi)
    l = 0; %Storage for correct knot
    m = [1,length(Xi)]; %Number of knots left to search
    if xi==1.0
        KnotNum = m(2);
    else
        while m(1) ~= m(2)
            if Xi(ceil((m(1)+m(2))/2))>xi
              m = [m(1),floor((m(1)+m(2))/2)];
            else
              m = [ceil((m(1)+m(2))/2),m(2)];
            end
        end
        KnotNum = m(2);
    end
end