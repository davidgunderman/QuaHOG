clear all
funct = @(x,y) funct2(x,y);
% load treble_clef_polynomial;
load guitar_polynomial
shape = guitar;
for j=1:length(guitar)
    SO=shape{j};
    SO(3:3:end,2:3)= 1+.3*(rand(length(SO)/3,2)-1);
    shape{j}=SO;
    plot_rat_bern_poly(SO,2,.001,"b.")
end
truev = ggPolygonIntegrate(shape,funct,50);
% cp=guitar{1}(1:3:end);
% truev=0;
% errg=zeros(0,length(shape));
% errs=zeros(0,length(shape));
% for ii=1:length(shape)
    for i=1:15
        errg(i)=ggPolygonIntegrate(shape, funct,i)-truev;
       errs(i)=sgPolygonIntegrate(shape,funct,i,5)-truev;
    end
% end
errs(errs==0)=1e-17;
errg(errg==0)=1e-17;
semilogy(1:15,abs(errg),'b.')
hold on
semilogy(1:15,abs(errs),'k.')
legend({"gauss","rational"})

function I = funct2(x,y)
% Function should be the moment of interest
%     if ~isempty(x)
%     scatter(x(:),y(:),'k.')
%     hold on
%     end
    I= ones(length(x),1);
%     I= x;
end