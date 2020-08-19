% In this file, we show the concept of solution transfer by integrating a
% function over some simple curved polygons with known solutions:
clear all;
close all;
addpath("../Rational_Quadrature/Matlab/Src",...
"../Rational_Quadrature/Matlab/Tests",...
"../Rational_Quadrature/Matlab/ThirdPartySupportingCode")
d=2;
gaussOrders=[3:8]
% Test cases: There are 3 functions that we will consider, two of which have
% known antiderivatives. These functions were taken from
% P.E. Farrell and J.R. Maddison. Conservative interpolation between 
% volume meshes by local Galerkin projection. Computer Methods in Applied 
% Mechanics and Engineering, 200(1-4):89–100, Jan 2011. 
RationalOn=0;

f = @(x,y) 5*y.^3 + x.^2 + 2*y +3; Afx= @(x,y) 5*y.^3.*x + ( 1/3 )*x.^3 + 2.*y.*x + 3*x; truev=10.21017612416683; 10.210176124869363
% f2= @(x,y) 2*x.^2 +y; Af2x = @(x,y) x.^2 + x.*y; truev=1.5707963267948966;
% f3 = @(x,y) 5*y.^3.*x.^3 + .2*x.^2 + 2*y.*x.^2 +3; Afx= @(x,y) 5*y.^3.*x + ( 1/3 )*x.^3 + 2.*y.*x + 3*x; truev=9.5818575934488694;

% g = @(x,y) 10*(exp( - x.^2 ) + 2*y); Agy = @(x,y) y.* exp( - x.^2 ) + y.^2;truev=25.17848513103563;
% h = @(x,y) sin( x ) + cos( y ); Agx = @(x,y) -cos( x ) + x.*cos( y ); truev=2.764919374768382;
% h2= @(x,y) (y.^3 - (x.^3.*y.^2) - (x.*y) -3)./((x.^2).*(y.^2) +100); truev=-0.09420855381254023;
% truev=26.85095435968714;
% Numintpoints=0;
field = @(x,y) field2(x,y,f);
% We will test the method by shaping an object into a 20 element Butterfly
% mesh defined by cubic boundaries, and a 16 element jittered cartesian 
% mesh defined by cubics:

meonet=load("Butterfly.mat"); Mesh1=meonet.refined;
nElemMesh1=length(Mesh1);

% Randomly rotate butterfly mesh
RotateMat= @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];
for i=1:nElemMesh1
    for j=1:2:size(Mesh1{i},1)
        Mesh1{i}(j:(j+1),:)=RotateMat(.175829069218)*Mesh1{i}(j:(j+1),:);
    end
end

RotateMat= @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];
for i=1:nElemMesh1
    for j=1:(size(Mesh1{i},1)/2)
        Mesh2{i}((3*j-2):(3*j-1),:)=1.15816729386*RotateMat(.175829069218)*Mesh1{i}((2*j-1):(2*j),:);
        Mesh2{i}(3*j,:)=ones(length(Mesh1{i}(2*j,:)),1);
    end
end
% We will find correct integrals for three different boundaries: an
% approximated circle, a cut hyperbola, and the letters H and M

PolyCircletemp=load("PolyCircle.mat"); PolyCircle=PolyCircletemp.PolyCircle;
RatCircletemp=load("Circle.mat"); RatCircle= RatCircletemp.Circle1;
Squaretemp=load("Square.mat"); Square=Squaretemp.Square;
% ShapeObject=RatCircle;
ShapeObject=PolyCircle*.89;
% ShapeObject(2:3:end,:)=
% ShapeObject(1:2,2)= -ShapeObject(1:2,2);
% ShapeObject(7:8,2)= 0;
% ShapeObject=.1672386*Square+3;
% ShapeObject(:,4)=0;

if RationalOn
    nEdShape=size(ShapeObject,1)/(d+1);
else
    nEdShape=size(ShapeObject,1)/d;
end

for i=1:nElemMesh1
    plot_bern_poly(Mesh1{i},2,.001,{},{},false)
    hold on
end
plot_bern_poly(ShapeObject,2,.001,{},{},false)

% for i=1:nElemMesh1
% %     Mesh1{i}=Mesh1{i};
%     plot_rat_bern_poly(Mesh2{i},2,.001,'k')
% end
% plot_rat_bern_poly(ShapeObject,2,.001,'k');
% Intersect each element, store moment of each material
fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,10);
mfanti=@(a,b) arrayfun(fanti,a,b);
SO{1}=ShapeObject; truev=PolygonIntegrate(SO,mfanti,10);
IntersectionI1=zeros(nElemMesh1,1);
for i=1:nElemMesh1
    clear Intersection
    if RationalOn
        Intersection=RatboolEls(ShapeObject,Mesh2{i},true);
    else
        Intersection=boolElsv2(ShapeObject,Mesh1{i},true);
    end
    if ~isempty(Intersection{1})
        s=plot_bern_poly(Intersection{1},2,.001,{},{},true);
%         plot_rat_bern_poly(Intersection{1},2,.1,'r');
        for j=gaussOrders
            fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,j);
            mfanti=@(a,b) arrayfun(fanti,a,b);
            if RationalOn
                IntersectionI1(i,j)=RatPolygonIntegrate(Intersection,mfanti,j);
            else
                IntersectionI1(i,j)=PolygonIntegrate(Intersection,mfanti,j);
            end
        end
    else
        IntersectionI1(i,:)=0;
    end
end
figure
% semilogy(gaussOrders,abs(sum(IntersectionI1(:,gaussOrders),1)-sum(IntersectionI1(:,gaussOrders(end)))))
semilogy(gaussOrders,abs(sum(IntersectionI1(:,gaussOrders),1)-truev))
% Error = (sum(IntersectionI)-integral2(field,-2*.26180283,2*.26180283,-2*.26180283,2*.26180283))./integral2(field,-2*.26180283,2*.26180283,-2*.26180283,2*.26180283,'AbsTol',0);
Error=abs(sum(IntersectionI1(:,gaussOrders),1)-truev);
figure
nplot=1000;
epts=1.5;
[x,y]= meshgrid([-epts:(epts/nplot):epts]+3,[(-epts:(epts/nplot):epts)']+3); x=x(:); y=y(:); xp=x; yp=y;
% xp(x.^2+y.^2>1)=nan; yp(x.^2+y.^2>1)=nan;
surf(reshape(xp,2*nplot+1,2*nplot+1),reshape(yp,2*nplot+1,2*nplot+1),zeros(2*nplot+1,2*nplot+1),field(reshape(xp,2*nplot+1,2*nplot+1),reshape(yp,2*nplot+1,2*nplot+1)),'edgecolor','none');
view([0 90])
hold on
% for i=1:nElemMesh1
% %     Mesh1{i}=Mesh1{i};
%     plot_bern_poly(Mesh1{i},2,.001,{},{'k'},false)
% end
% plot_bern_poly(ShapeObject,2,.001,{},{'k'},false)
% for i=1:nElemMesh1
% %     Mesh1{i}=Mesh1{i};
%     plot_rat_bern_poly(Mesh2{i},2,.001,'k')
% end
% plot_rat_bern_poly(ShapeObject,2,.001,'b')
printError= floor(log(Error)/log(10));
title({sprintf('Background function: $5y^3 + x^2 + 2y +3$, Error $\\approx 10^{%d}$',printError),sprintf('Quadrature points per side of intersection: $%d^2$',gaussOrders(end-1)-1)},'interpreter','latex','FontSize',16)
% title({sprintf('Background function: $\\frac{y^3 - x^3 y^2 - xy -3}{x^2y^2 + 100}$, Error $\\approx 10^{%d}$',printError),sprintf('Quadrature points per side of intersection: $%d^2$',gaussOrders(end-1)-1)},'interpreter','latex','FontSize',16)
axis off
colorbar
title({sprintf('Background function: $\\frac{y^3 - x^3 y^2 - xy -3}{x^2y^2 + 100}$')},'interpreter','latex','FontSize',16)
