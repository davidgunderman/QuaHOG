function I = sgPolygonIntegrate(Cellintel,funct,pts,pp,kg)
   % Funct is the x antiderivative of some function defined in the interior
   % of Cellintel. Cellintel contains the boundaries of the region of
   % interest, ordered counterclockwise and expressed as Bezier Curves
   % pts is the number of total quadrature points in each direction, pp is
   % the number of total poles (kg*m+2*m) and kg is the gauss quad pts
   % green's theorem antiderivative
   I=0;
   I2=0;
   
   for i=1:length(Cellintel)
       cp=Cellintel{i}(1:3:end,:);
        fanti=@(a,b) gauss1D(@(x)funct(x,b),min(min(cp(:))),a,kg);
        mfanti=@(a,b) arrayfun(fanti,a,b);
       for j=1:3:(size(Cellintel{i},1))
            newfunct= @(t) RatGreensFunction(t,Cellintel{i}(j:(j+2),:),mfanti);
            rts=roots(fliplr(BernsteinToMonomial(Cellintel{i}(j+2,:))));
            sgl=repmat(rts,ceil(pp/(size(Cellintel{1},2)-1)),1);
            sgl((length(sgl)+1):(pts),:)=[Inf];
            [fout,sglout ] = transf( @(x)newfunct(x) , sgl , [0,1] );
%             I=I+rfejer(sglout,fout);
            if isempty(rts)
                I=I+gauss1D(newfunct,0,1,pp+pts);
            else
                x = rfejer(sglout);
                if any(x(2,:)<0)
                    I=I+gauss1D(newfunct,0,1,pp+pts);
                elseif all(abs(Cellintel{i}(j+2,:)-1)<1e-15)
                    I=I+gauss1D(newfunct,0,1,pp+pts);
                    else
                I = I + fout(x(1,:))*x(2,:)';
                end
            end
       end
   end
end

function I = RatGreensFunction(t,Side,funct)
    xy=dCR_eval(Side,t);
    x=xy(:,1);y=xy(:,2);
    yp = dCR_eval_dt(Side(2:3,:),t);
    I=(funct(x,y).*yp)';
end

% a=5; % horizontal radius
% b=10; % vertical radius
% x0=Centroid(1); % x0,y0 ellipse centre coordinates
% y0=Centroid(2);
% t=-pi:0.01:pi;
% x=s(4)*cos(t);
% y=s(1)*sin(t);
% xy=v*[x;y] +[x0;y0];
% hh=plot(xy(1,:),xy(2,:),'k','Linewidth',2)