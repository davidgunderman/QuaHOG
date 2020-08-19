function quadpts = gauss2Dbdspts(bds,pts)
%gauss2Dbds Returns tensor product quadature pts and wts over the rectangle
%defined by bds with gaussian quadrature of order pts
%   Inputs:
%   funct       integrand function given as @(x,y)...
%   bds         [xmin xmax, ymin, ymax] boundaries of rectangular domain
%   pts         Order of Gaussian quadrature to be used (same in both dirs)
%
%   Outputs:
%   I           value of integral
    [x, wx] = GaussLegendreQuad1D(pts,bds(1),bds(2));
    [y, wy] = GaussLegendreQuad1D(pts,bds(3),bds(4));
    [s,t]=meshgrid(x,y);
    ss=s(:); tt=t(:);
    w2=wx'*wy;
    quadpts= [ss,tt,w2(:)];
end