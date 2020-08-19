function I = sgPolygonIntegrate(Cellintel,funct,pts,pp)
   % Funct is the x antiderivative of some function defined in the interior
   % of Cellintel. Cellintel contains the boundaries of the region of
   % interest, ordered counterclockwise and expressed as Bezier Curves
   I=0;
   I2=0;
   for i=1:length(Cellintel)
       for j=1:3:(size(Cellintel{i},1))
            newfunct= @(t) RatGreensFunction(t,Cellintel{i}(j:(j+2),:),funct);
            sgl=repmat(roots(fliplr(BernsteinToMonomial(Cellintel{i}(j+2,:)))),ceil(pp/2),1);
            sgl((length(sgl)+1):(pts),:)=[Inf];
            [fout,sglout ] = transf( @(x)newfunct(x) , sgl , [0,1] );
            x = rfejer(sglout);
            I = I + fout(x(1,:))*x(2,:)';
%             I=I+gauss1D(newfunc,0,1,pts);
       end
   end
end

function I = RatGreensFunction(t,Side,funct)
    xy=dCR_eval(Side,t);
    x=xy(:,1);y=xy(:,2);
    yp = dCR_eval_dt(Side(2:3,:),t);
    I=(funct(x,y).*yp)';
end