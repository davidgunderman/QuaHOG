% This file's purpose is to compare quadtree, gauss-green, and szego-green...
% quadrature methods using a small suite of test cases defined by two
% domains and various integrands:
% Domains:
shape=0;
testIntegrands=1;
% 1. A circular region whose boundary is defined by four rational curves
lroffset=[.3957106819596820 -1.5728593603867];
Circletemp=load("PolyCircle.mat"); Circle= Circletemp.PolyCircle;
Circle(1:2:end,:)=Circle(1:2:end,:)+lroffset(1);
Circle(2:2:end,:)=Circle(2:2:end,:)+lroffset(2);
% 2. A region defined by the intersection of two circular regions.
offset=.1364728595817;
C1=Circle; C1(2:2:end,:)=C1(2:2:end,:)+offset;
C2=Circle; C2(2:2:end,:)=C2(2:2:end,:)-offset;
InterCircles=boolElsv2(C1,C2,true); InterCircles=InterCircles{1};
% plot_rat_bern_poly(InterCircles{1},2,.001,'k')
% Integrands:
% 1. Monomials up to 6th degree: (1, x, y, ...)
monfuncts={@(x,y) ones(size(x));
           @(x,y) x;
           @(x,y) y;
           @(x,y) x.^2;
           @(x,y) x.*y;
           @(x,y) y.^2;
           @(x,y) x.^3;
           @(x,y) x.^2.*y;
           @(x,y) x.*y.^2;
           @(x,y) y.^3;
           @(x,y) x.^4;
           @(x,y) x.^3.*y;
           @(x,y) x.^2.*y.^2;
           @(x,y) x.*y.^3;
           @(x,y) y.^4;
           @(x,y) x.^5;
           @(x,y) x.^4.*y;
           @(x,y) x.^3.*y.^2;
           @(x,y) x.^2.*y.^3;
           @(x,y) x.*y.^4;
           @(x,y) y.^5;
           @(x,y) x.^6;
           @(x,y) x.^5.*y;
           @(x,y) x.^4.*y.^2;
           @(x,y) x.^3.*y.^3;
           @(x,y) x.^2.*y.^4;
           @(x,y) x.*y.^5;
           @(x,y) y.^6;}
% 2. Three polynomials of degree 2 (bilinear), 4 (biquadratic), and 6
% (bicubic)
polyfuncts={@(x,y) (2*x.^2 +x.*y - y +2);
            @(x,y) (2*x.^2.*y.^2 +.3*x.^2.*y - y.^4 + 3*x +2);
            @(x,y) (x.^5 - 5*y.^3.*x.^3 + .2*x.^2 + 2*y.*x.^2 +3);}
% 3. A rational function of degree 4 and an exponential function.
otherfuncts={@(x,y) (y.^3 - (x.^3.*y.^2) - (x.*y) -3)./((x.^2).*(y.^2) +10);
             @(x,y) 10*(exp( - x.^2 ) + 2*y);
             @(x,y) sqrt((x+10).^2+(x+10).*(y+10).^2 +x)}
addpath("../Rational_Quadrature/Matlab/Src",...
"../Rational_Quadrature/Matlab/Tests",...
"../Rational_Quadrature/Matlab/ThirdPartySupportingCode")
d=2;
gaussOrders=[2:15]
% Test cases: There are 3 functions that we will consider, two of which have
% known antiderivatives. These functions were taken from
if shape==0
    shapeObject=Circle;
elseif shape==1
    shapeObject=InterCircles;
end

if testIntegrands==0
    integrands=monfuncts;
elseif testIntegrands==1
    integrands=polyfuncts;
elseif testIntegrands==2
    integrands=otherfuncts;
end
% figure(7)
% plot_bern_poly([C1;C2],2,.001,{},{},false)
% hold on
% plot_bern_poly(shapeObject,2,.001,{},{},false)
figure(8)
% plot_bern_poly([C1;C2],2,.001,{},{},false)
% hold on
plot_bern_poly(shapeObject,2,.001,{},{},false)

RationalOn=0;
numIntegrands=length(integrands);
elemSize=size(shapeObject,1)/3;

int2evals=zeros(length(numIntegrands),1);
int2errs=zeros(length(numIntegrands),1);
ggevals=int2evals;
sgevals=int2evals;
global evalCounter;
global scatterEvals;
scatterEvals=0;
for i=1:numIntegrands
    figure(1)
    field = @(x,y) field2(x,y,integrands{i});
%     if shape==0
%         truev = integral2(field, -1+lroffset(1),1+lroffset(1), @(x)-sqrt(1-(x-lroffset(1)).^2)+lroffset(2), @(x)sqrt(1-(x-lroffset(1)).^2)+lroffset(2),'AbsTol',1e-17,'RelTol',1e-18);
%     elseif shape==1
%         truev = integral2(field, -sqrt(1-offset^2)+lroffset(1), sqrt(1-offset^2)+lroffset(1), @(x)-sqrt(1-(x-lroffset(1)).^2)+offset+lroffset(2),@(x)sqrt(1-(x-lroffset(1)).^2)-offset+lroffset(2),'AbsTol',1e-17,'RelTol',1e-18);
%     end
% %     fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,15);
%     mfanti=@(a,b) arrayfun(fanti,a,b);
%     SO{1}=shapeObject; truev=RatPolygonIntegrate(SO,mfanti,24,16);
%     for jj=1:14
%         evalCounter=0;
%         if jj==7
%             scatterEvals=1;
%         else
%             scatterEvals=0;
%         end
%         if shape==0
%             int2errs(jj) = truev-integral2(field, -1+lroffset(1),1+lroffset(1), @(x)-sqrt(1-(x-lroffset(1)).^2)+lroffset(2), @(x)sqrt(1-(x-lroffset(1)).^2)+lroffset(2),'RelTol',10^(-jj));
%             int2evals(jj)=evalCounter;
%         elseif shape==1
%             int2errs(jj) = truev-integral2(field, -sqrt(1-offset^2)+lroffset(1), sqrt(1-offset^2)+lroffset(1), @(x)-sqrt(1-(x-lroffset(1)).^2)+offset+lroffset(2),@(x)sqrt(1-(x-lroffset(1)).^2)-offset+lroffset(2),'RelTol',10^(-jj+3));
%             int2evals(jj)=evalCounter;
%         end
%     end
%     int2errs(int2errs==0)=1e-17;
%     figure(1)
%     semilogy(int2evals,abs(int2errs),'k.','MarkerSize',36)
    hold on
%     plot_rat_bern_poly(shapeObject,2,.001,'k');
    % Intersect each element, store moment of each material
%     fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,7);
%     mfanti=@(a,b) arrayfun(fanti,a,b);
    SO{1}=shapeObject; 
%     truev=PolygonIntegrate(SO,mfanti,7);
%     IntersectionI1=zeros(nElemMesh1,1);

%     s=plot_bern_poly(Intersection{1},2,.001,{},{},true);
    %         plot_rat_bern_poly(Intersection{1},2,.1,'r');
    RationalOn=1;
    ggevals=zeros(length(gaussOrders),1);
    ggerrs=zeros(length(gaussOrders),1);
    sgevals=zeros(length(gaussOrders),1);
    sgerrs=zeros(length(gaussOrders),1);
    for jj=2:length(gaussOrders)
        j=gaussOrders(jj);
        if testIntegrands==0
            kk=2*invTri(i)+6;
            kg=max(ceil(invTri(i)+1),2);
        elseif testIntegrands==1
            kk=4*i+6;
            kg=max(i,2);
        else
            kg=4;
            kk=j;
        end
        evalCounter=0;
        fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,kg);
    %     intxw= @(a,b) gaussXW(@(x)field(x,b),0,a,j);
        mfanti=@(a,b) arrayfun(fanti,a,b);
    %     mintxw= @(a,b) arrayfun(intxw,a,b);
        if jj==2
            scatterEvals=1;
        else
            scatterEvals=0;
        end
        ggerrs(jj)=PolygonIntegrate(SO,mfanti,j);
        ggevals(jj)=evalCounter;
%         evalCounter=0;
%         sgerrs(jj)=sgPolygonIntegrate(SO,mfanti,j-1,kk)-truev;
%         sgevals(jj)=evalCounter;
    end
    figure(1)
    figure;
    semilogy(ggevals,abs(ggerrs-ggerrs(end)),'b.','MarkerSize',36);
end
    
%     semilogy(sgevals,abs(sgerrs),'g.','MarkerSize',36)
%     xlim([0,max(max(ggevals)*3,3*min(int2evals))])
% end
% Error = (sum(IntersectionI)-integral2(field,-2*.26180283,2*.26180283,-2*.26180283,2*.26180283))./integral2(field,-2*.26180283,2*.26180283,-2*.26180283,2*.26180283,'AbsTol',0);
% Error=abs(sum(IntersectionI1(:,gaussOrders),1)-truev);
% figure
% nplot=1000;
% epts=1.5;
% [x,y]= meshgrid([-epts:(epts/nplot):epts]+3,[(-epts:(epts/nplot):epts)']+3); x=x(:); y=y(:); xp=x; yp=y;
% xp(x.^2+y.^2>1)=nan; yp(x.^2+y.^2>1)=nan;
% surf(reshape(xp,2*nplot+1,2*nplot+1),reshape(yp,2*nplot+1,2*nplot+1),zeros(2*nplot+1,2*nplot+1),field(reshape(xp,2*nplot+1,2*nplot+1),reshape(yp,2*nplot+1,2*nplot+1)),'edgecolor','none');
% view([0 90])
% hold on
% for i=1:nElemMesh1
% %     Mesh1{i}=Mesh1{i};
%     plot_bern_poly(Mesh1{i},2,.001,{},{'k'},false)
% end
% plot_bern_poly(shapeObject,2,.001,{},{'k'},false)
% for i=1:nElemMesh1
% %     Mesh1{i}=Mesh1{i};
%     plot_rat_bern_poly(Mesh2{i},2,.001,'k')
% end
% plot_rat_bern_poly(shapeObject,2,.001,'b')
printError= floor(log(Error)/log(10));
title({sprintf('Background function: $5y^3 + x^2 + 2y +3$, Error $\\approx 10^{%d}$',printError),sprintf('Quadrature points per side of intersection: $%d^2$',gaussOrders(end-1)-1)},'interpreter','latex','FontSize',16)
% title({sprintf('Background function: $\\frac{y^3 - x^3 y^2 - xy -3}{x^2y^2 + 100}$, Error $\\approx 10^{%d}$',printError),sprintf('Quadrature points per side of intersection: $%d^2$',gaussOrders(end-1)-1)},'interpreter','latex','FontSize',16)
axis off
colorbar
title({sprintf('Background function: $\\frac{y^3 - x^3 y^2 - xy -3}{x^2y^2 + 100}$')},'interpreter','latex','FontSize',16)

function xw = gaussXW(bound1,bound2,pts)
    % 15 point gauss quadrature weights and nodes on interval [-1,1]
    gaussQuad=load("gaussQuad");
    w=gaussQuad.wv{pts-1};
    x=gaussQuad.abc{pts-1};
    scale=bound2-bound1;
    w=w*scale/2;
    x=(scale/2)*(x+1)+bound1;
    xw=[x w]
end

function kk = invTri(i)
    kk=floor(real((sqrt(1/4+2*(i-1))-1/2)));
end