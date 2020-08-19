% This file's purpose is to compare quadtree, gauss-green, and szego-green...
% quadrature methods using a small suite of test cases defined by two
% domains and various integrands:
% Domains:
clear all;
shape=0;
testIntegrands=0;
% 1. A circular region whose boundary is defined by four rational curves
lroffset=[.3957106819596820 -1.5728593603867];
Circletemp=load("Circle.mat"); Circle= Circletemp.Circle1;
Circle(1:3:end,:)=Circle(1:3:end,:)+lroffset(1).*Circle(3:3:end,:);
Circle(2:3:end,:)=Circle(2:3:end,:)+lroffset(2).*Circle(3:3:end,:);
% 2. A region defined by the intersection of two circular regions.
offset=.1364728595817;
C1=Circle; C1(2:3:end,:)=C1(2:3:end,:)+offset.*C1(3:3:end,:);
C2=Circle; C2(2:3:end,:)=C2(2:3:end,:)-offset.*C1(3:3:end,:);
C3=Circle; C3(1:3:end,:)=C3(1:3:end,:)+2*offset.*C3(3:3:end,:);
InterCirclestemp=RatboolEls(C1,C2,true); InterCircles=InterCirclestemp{1};
InterCirclesShift=InterCircles; InterCirclesShift(2:3:end,:)=InterCirclesShift(2:3:end,:)+(offset).*InterCirclesShift(3:3:end,:);
WeirdShapetemp=RatboolEls(InterCirclesShift,C3,false); WeirdShape=WeirdShapetemp{1};
WeirdShape(1:3:end,:)=WeirdShape(1:3:end,:)-.5*WeirdShape(3:3:end,:);
WeirdShape(2:3:end,:)=WeirdShape(2:3:end,:)+1.5*WeirdShape(3:3:end,:);
Circles{1}=Circle;

weirdObjectt=load("weirdObject"); weirdObject= weirdObjectt.fig1;
for kk=2:6
        Circles{kk}=zeros(12,kk+2);
    for iii=1:3:12
        transformtemp=Circles{kk-1}(iii:(iii+2),:);
%         transformtemp(1:2,:)=transformtemp(1:2,:).*transformtemp(3,:);
        transformtemp2=zeros(3,kk+2);
        transformtemp2(:,1)=transformtemp(:,1);
        transformtemp2(:,kk+2)=transformtemp(:,kk+1);
        for kkk=2:(kk+1)
            transformtemp2(:,kkk)=(kkk-1)/(kk+1)*transformtemp(:,kkk-1)+(1-((kkk-1)/(kk+1)))*transformtemp(:,kkk);
        end
%         transformtemp2(1:2,:)=transformtemp2(1:2,:)./transformtemp2(3,:);
        Circles{kk}(iii:(iii+2),:)=transformtemp2;
    end
end

% plot_rat_bern_poly(InterCircles,2,.001,'k')
% plot_rat_bern_poly(WeirdShape,2,.001,'k')
% Integrands:
% 1. Monomials up to 6th degree: (1, x, y, ...)
counter=1;
for i=0:5
    for j=0:i
        a=i; b=j;
        monfuncts(counter) = {@(x,y) x.^a.*y.^(a-b)};
        counter=counter+1;
    end
end
% 2. Three polynomials of degree 2 (bilinear), 4 (biquadratic), and 6
% (bicubic)
polyfuncts={@(x,y) (2*x.^2 +x.*y - y +2);
            @(x,y) (2*x.^2.*y.^2 +.3*x.^2.*y - y.^4 + 3*x +2);
            @(x,y) (x.^5 - 5*y.^3.*x.^3 + .2*x.^2 + 2*y.*x.^2 +3);};
% 3. A rational function of degree 4 and an exponential function.
otherfuncts={@(x,y) (y.^3 - (x.^3.*y.^2) - (x.*y) -3)./((x.^2).*(y.^2) +10);
             @(x,y) 10*(exp( - x.^2 ) + 2*y);
             @(x,y) sqrt((x+10).^2+(x+10).*(y+10).^2 +x)};
addpath("../Rational_Quadrature/Matlab/Src",...
"../Rational_Quadrature/Matlab/Tests",...
"../Rational_Quadrature/Matlab/ThirdPartySupportingCode")
d=2;
gaussOrders=[3:3];
% Test cases: There are 3 functions that we will consider, two of which have
% known antiderivatives. These functions were taken from


if testIntegrands==0
    integrands=monfuncts;
elseif testIntegrands==1
    integrands=polyfuncts;
elseif testIntegrands==2
    integrands=otherfuncts;
end
monfunctmat=zeros(length(integrands),6);
for jjj=1:1
if shape==0
    shapeObject=Circles(1);
elseif shape==1
    shapeObject=InterCircles;
elseif shape==2
    shapeObject=WeirdShape;
elseif shape==3
%     shapeObject=weirdObject;
end
% figure(7)
% plot_rat_bern_poly([C1;C2],2,.001,"k.")
% hold on
% plot_rat_bern_poly(shapeObject,2,.001,"b.")
% figure(8)
% plot_rat_bern_poly([C1;C2],2,.001,"k.")
% hold on
% plot_rat_bern_poly(shapeObject,2,.001,"b.")

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

for i=1:length(integrands)
%     figure(1)
    field = @(x,y) field2(x,y,integrands{i});
    if shape==0
        truev = integral2(field, -1+lroffset(1),1+lroffset(1), @(x)-sqrt(1-(x-lroffset(1)).^2)+lroffset(2), @(x)sqrt(1-(x-lroffset(1)).^2)+lroffset(2),'AbsTol',1e-17,'RelTol',1e-18);
    elseif shape==1
        truev = integral2(field, -sqrt(1-offset^2)+lroffset(1), sqrt(1-offset^2)+lroffset(1), @(x)-sqrt(1-(x-lroffset(1)).^2)+offset+lroffset(2),@(x)sqrt(1-(x-lroffset(1)).^2)-offset+lroffset(2),'AbsTol',1e-17,'RelTol',1e-18);
    elseif shape==2
        fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,15);
        mfanti=@(a,b) arrayfun(fanti,a,b);
        SO{1}=shapeObject; truev=ggPolygonIntegrate(SO,mfanti,15);
    elseif shape==2
    fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,15);
    mfanti=@(a,b) arrayfun(fanti,a,b);
    SO{1}=shapeObject; truev=ggPolygonIntegrate(SO,mfanti,15);
    elseif shape==3
    fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,15);
    mfanti=@(a,b) arrayfun(fanti,a,b);
    SO{1}=shapeObject; truev=ggPolygonIntegrate(SO,mfanti,15);
    end
    for jj=1:14
        evalCounter=0;
%         if jj==2
%             scatterEvals=1;
%         else
%             scatterEvals=0;
%         end
        if shape==0
            int2errs(jj) = truev-integral2(field, -1+lroffset(1),1+lroffset(1), @(x)-sqrt(1-(x-lroffset(1)).^2)+lroffset(2), @(x)sqrt(1-(x-lroffset(1)).^2)+lroffset(2),'RelTol',10^(-jj));
            int2evals(jj)=evalCounter;
        elseif shape==1
            int2errs(jj) = truev-integral2(field, -sqrt(1-offset^2)+lroffset(1), sqrt(1-offset^2)+lroffset(1), @(x)-sqrt(1-(x-lroffset(1)).^2)+offset+lroffset(2),@(x)sqrt(1-(x-lroffset(1)).^2)-offset+lroffset(2),'RelTol',10^(-jj+3));
            int2evals(jj)=evalCounter;
        end
    end
    int2errs(int2errs==0)=1e-17;
%     figure(1)
%     semilogy(int2evals,abs(int2errs),'k.','MarkerSize',36)
%     hold on
%     plot_rat_bern_poly(shapeObject,2,.001,'k');
    % Intersect each element, store moment of each material
%     fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,7);
%     mfanti=@(a,b) arrayfun(fanti,a,b);
    SO{1}=shapeObject; 

%     s=plot_bern_poly(Intersection{1},2,.001,{},{},true);
    %         plot_rat_bern_poly(Intersection{1},2,.1,'r');
    RationalOn=1;
    ggevals=zeros(length(gaussOrders),1);
    ggerrs=zeros(length(gaussOrders),1);
    sgevals=zeros(length(gaussOrders),1);
    sgerrs=zeros(length(gaussOrders),1);
    for jj=1:length(gaussOrders)
        j=gaussOrders(jj);
        if testIntegrands==0
            kk=(jjj+1)*invTri(i)+3*(jjj+1);
            kg=ceil(invTri(i)+1);
        elseif testIntegrands==1
            kk=4*i+6;
            kg=i;
        else
            kg=4;
            kk=j;
        end
        evalCounter=0;
        fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,kg);
    %     intxw= @(a,b) gaussXW(@(x)field(x,b),0,a,j);
        mfanti=@(a,b) arrayfun(fanti,a,b);
    %     mintxw= @(a,b) arrayfun(intxw,a,b);
        if jj==3
            scatterEvals=1;
        else
            scatterEvals=0;
        end
        evalCounter=0;
        ggerrs(jj)=ggPolygonIntegrate(SO,mfanti,j)-truev;
        ggevals(jj)=evalCounter;
        evalCounter=0;
        sgerrs(jj)=sgPolygonIntegrate(SO,mfanti,j-1,kk)-truev;
        sgevals(jj)=evalCounter;
    end
    monfunctmat(i,jjj)=[ggerrs(find(abs(ggevals)==sgevals(1),1))/sgerrs(1)]
%     figure(i)
%     semilogy(ggevals,abs(ggerrs),'b.','MarkerSize',36)
%     hold on
%     semilogy(sgevals,abs(sgerrs),'g.','MarkerSize',36)
%     xlim([0,max(max(ggevals)*3,min(int2evals))])
end

end
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
% printError= floor(log(Error)/log(10));
% title({sprintf('Background function: $5y^3 + x^2 + 2y +3$, Error $\\approx 10^{%d}$',printError),sprintf('Quadrature points per side of intersection: $%d^2$',gaussOrders(end-1)-1)},'interpreter','latex','FontSize',16)
% % title({sprintf('Background function: $\\frac{y^3 - x^3 y^2 - xy -3}{x^2y^2 + 100}$, Error $\\approx 10^{%d}$',printError),sprintf('Quadrature points per side of intersection: $%d^2$',gaussOrders(end-1)-1)},'interpreter','latex','FontSize',16)
% axis off
% colorbar
% title({sprintf('Background function: $\\frac{y^3 - x^3 y^2 - xy -3}{x^2y^2 + 100}$')},'interpreter','latex','FontSize',16)
figure;
if isempty(find(monfunctmat(:,1)==0))
    lasti=size(monfunctmat,1)
else
    lasti=find(monfunctmat(:,1)==0);
end
for jjj=1:1
    for i=1:lasti
        ertemp=monfunctmat(find(invTri(1:length(monfuncts))==i-1),1);
        er(i)=mean(ertemp(ertemp~=0 & ertemp~=Inf & ~isnan(ertemp)));
    end
plot(0:(lasti-1),log(abs(er))/log(10));
end

% function xw = gaussXW(bound1,bound2,pts)
%     % 15 point gauss quadrature weights and nodes on interval [-1,1]
%     gaussQuad=load("gaussQuad");
%     w=gaussQuad.wv{pts-1};
%     x=gaussQuad.abc{pts-1};
%     scale=bound2-bound1;
%     w=w*scale/2;
%     x=(scale/2)*(x+1)+bound1;
%     xw=[x w]
% end
