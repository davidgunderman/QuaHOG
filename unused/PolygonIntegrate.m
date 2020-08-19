function I = PolygonIntegrate(Cellintel,funct,pts)
   % Funct is the x antiderivative of some function defined in the interior
   % of Cellintel. Cellintel contains the boundaries of the region of
   % interest, ordered counterclockwise and expressed as Bezier Curves
   I=0;
   for i=1:length(Cellintel)
       for j=1:2:(size(Cellintel{i},1))
          newfunct= @(t) GreensFunction(t,Cellintel{i}(j:(j+1),:),funct);
          Inow=gauss1D(newfunct,0,1,pts);
          I=I+Inow;
       end
   end
end

function I = GreensFunction(t,Side,funct)
    xy=dC_eval(Side,t);
    x=xy(:,1);y=xy(:,2);
    yp = dC_eval_dt(Side(2,:),t);
    I=(funct(x,y).*yp)';
%     [x y yp I]
end