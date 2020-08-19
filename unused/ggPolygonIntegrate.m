function I = ggPolygonIntegrate(Cellintel,funct,pts,kg)
   % Funct is the x antiderivative of some function defined in the interior
   % of Cellintel. Cellintel contains the boundaries of the region of
   % interest, ordered counterclockwise and expressed as Bezier Curves
   
   I=0;
   I2=0;
   for i=1:length(Cellintel)
       cp=Cellintel{i}(1:3:end,:);
        fanti=@(a,b) gauss1D(@(x)funct(x,b),min(min(cp(:))),a,kg);
        mfanti=@(a,b) arrayfun(fanti,a,b);
       for j=1:3:(size(Cellintel{i},1))
            newfunct= @(t) RatGreensFunction(t,Cellintel{i}(j:(j+2),:), mfanti);
            I=I+gauss1D(newfunct,0,1,pts);
       end
   end
end

function I = RatGreensFunction(t,Side,funct)
    xy=dCR_eval(Side,t);
    x=xy(:,1);y=xy(:,2);
    yp = dCR_eval_dt(Side(2:3,:),t);
    I=(funct(x,y).*yp)';
end