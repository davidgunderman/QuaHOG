function I = gauss1D(funct,bound1,bound2,pts)
    I=0;
    [x, w] = GaussLegendreQuad1D(pts,bound1,bound2);
    for i=1:size(x,1)
        I=I+ w(i)*funct(x(i));
    end
    I=I;
end