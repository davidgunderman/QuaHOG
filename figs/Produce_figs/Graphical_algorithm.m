close all

Circletemp=load("Circle.mat"); Circle= Circletemp.Circle1;
exampleCurve=Circle(1:3,:);
EC{1}=exampleCurve;
exampleIntegrand=@(x,y) x.^1.*y.^(2);
antiDerivative=@(x,y) .5*x.^2.*y.^2;
CPs=exampleCurve(1:2,:)./exampleCurve(3,:);
xx0=linspace(-1.1,1.1,100);
yy0=linspace(-1.1,1.1,100);
[xp0,yp0]=meshgrid(xx0,yy0);
figure('Position',[10 10 1200 300])
subplot(1,5,1);
xpeval=xp0;
ypeval=yp0;
xpeval(sqrt(xp0.^2 + yp0.^2)>1)=nan;
ypeval(sqrt(xp0.^2 + yp0.^2)>1)=nan;
pcolor(xpeval,ypeval,exampleIntegrand(xpeval,ypeval));
shading interp
hold on
plot_rat_bern_poly(Circle,2,.001,"k.")
xlim([-1.1 1.1]); ylim([-1.1 1.1]);
axis square
axis off
xmin=min(CPs(1,:))-.1; xmax=max(CPs(1,:))+.1;
ymin=min(CPs(2,:))-.1; ymax=max(CPs(2,:))+.1;
xx1=linspace(xmin,xmax,100);
yy1=linspace(ymin,ymax,100);
[xp1,yp1]=meshgrid(xx1,yy1);
xx2=linspace(0,1,5);
yy2=linspace(0,1,5);
[xp2,yp2]=meshgrid(xx2,yy2);
% xp=xp(:); yp=yp(:);
subplot(1,5,2)
pcolor(xp1,yp1,exampleIntegrand(xp1,yp1));
axis square
shading interp
hold on
plot_rat_bern_poly(exampleCurve,2,.001,'k.');
axis off
xlim([xmin/1.05 xmax*1.05]);
ylim([ymin/1.05 ymax*1.05]);
subplot(1,5,3)
% quiver(xp2,yp2,zeros(5,5),antiDerivative(xp2,yp2),'k','LineWidth',2,'MaxHeadSize',14)
quiverc(xp2,yp2,zeros(5,5),antiDerivative(xp2,yp2),1.8)
set(gca,'Color',[1 1 1]);
set(gcf,'Color',[0.9400    0.9400    0.9400]);
hold on
plot_rat_bern_poly(exampleCurve,2,.001,'k.');
axis square
axis off
axis off
xlim([xmin/1.3 xmax*1.3]);
ylim([ymin/1.3 ymax*1.3]);
subplot(1,5,4)
quiverc(xp2,yp2,zeros(5,5),antiDerivative(xp2,yp2),1.8)
set(gca,'Color',[1 1 1]);
set(gcf,'Color',[0.9400    0.9400    0.9400]);
axis square
hold on
plot_rat_bern_poly(exampleCurve,2,.001,'k.');
axis off
xy=[0.9945    0.1050
    0.7071    0.7071
    0.1050    0.9945]
scatter(xy(:,1),xy(:,2),500,'r.')
xlim([xmin/1.3 xmax*1.3]);
ylim([ymin/1.3 ymax*1.3]);
subplot(1,5,5)
newx=zeros(3,3);
for i=1:3
    [xw, w] = GaussLegendreQuad1D(3,0,xy(i,1));
    newx(i,:)=xw;
end
pcolor(xp1,yp1,exampleIntegrand(xp1,yp1));
shading interp
axis square
axis off
hold on
plot_rat_bern_poly(exampleCurve,2,.001,'k.');
ss=scatter(newx(:),repmat(xy(:,2),3,1),500,[.5 0 0],'.')
xlim([xmin/1.05 xmax*1.05]);
ylim([ymin/1.05 ymax*1.05]);