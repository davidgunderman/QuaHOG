function [xquad,yquad,wquad] = SPECTRAL_quads(Cellintel,pts,kg)
% Calculates quadrature points for the region defined by Cellintel using
% the "SPECTRAL" method from the paper "spectral mesh-free quadrature for
% planar regions bounded by rational parametric curves."
% Input:
% Cellintel      Oriented boundaries of the region
% pts            Number of intermediate quadrature points to use
% kg             Number of antiderivative quadrature points to use
%
% Output:
% [xquad, yquad] Quadrature points locations
% wquad          Quadrature weights
   xquad=zeros(0,1); yquad=zeros(0,1); wquad=zeros(0,1);
   for i=1:length(Cellintel)
       cp=Cellintel{i}(1:3:end,:);
       for j=1:3:(size(Cellintel{i},1))
            [xg, wg ]=GaussLegendreQuad1D(pts,0,1);
            x=[xg';wg];
            [xc,yc,yp]= RatGreensFunctionQuad(x(1,:),Cellintel{i}(j:(j+2),:));
            clear xq; clear wq; clear yq;
            for ii=1:length(xc)
                    [xxq,wwq]=GaussLegendreQuad1D(kg,min(min(cp(:))),xc(ii));
                    xq(ii,:)=xxq;
                    wq(ii,:)=wwq.*yp(ii).*x(2,ii);
            end
            yq=repmat(yc,1,kg);
            xquad=[xquad; xq(:)];
            yquad=[yquad; yq(:)];
            wquad=[wquad; wq(:)];
       end
   end
end

function [x,y,yp] = RatGreensFunctionQuad(t,Side)
% Finds x(t),y(t), and y'(t) on a Bezier Curve defined in Side
% Input:
% t         Parameter values at which evaluation should occur
% Side      Rational Bezier Curve
%
% Output:
% [x,y,yp]  x,y, and dy/ds values for "Side" at t-values
    for i=1:length(t)
        xy=dCR_eval(Side,t(i));
        x(i,:)=xy(:,1)';y(i,:)=xy(:,2)';
        yp(i,:) = dCR_eval_dt(Side(2:3,:),t(i));
    end
end


