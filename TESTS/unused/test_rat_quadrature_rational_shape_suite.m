% This file's purpose is to compare quadtree, gauss-green, and szego-green...
% quadrature methods using a small suite of test cases defined by two
% domains and various integrands:
% Domains:
close all;
clear all;
shape=3;
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

testFig=load("testFig3"); testFigs= testFig.testFig;

% plot_rat_bern_poly(InterCircles,2,.001,'k')
% plot_rat_bern_poly(WeirdShape,2,.001,'k')
% Integrands:
% 1. Monomials up to 6th degree: (1, x, y, ...)
counter=1;
for i=0:5
    for j=0:i
        a=i-j; b=j;
        monfuncts(counter) = {@(x,y) x.^a.*y.^(b)};
        bla(counter,:)=[a,b];
        counter=counter+1;
    end
end
monfuncts(7)= {@(x,y) (2*x.^2 +x.*y - y +2)};
monfuncts(8)= {@(x,y) sqrt((x+10).^2+(x+10).*(y+10).^2 +x)};
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
gaussOrders=[2:25];
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
monfunctquad=zeros(length(integrands),6);

% if shape==0
    shapeObject{1}=Circle;
%     SO{1}=Circle;
% elseif shape==1
%     shapeObject=InterCircles;
% elseif shape==2
%     shapeObject=WeirdShape;
% elseif shape==3
%     shapeObject=weirdObject;
% end
% figure(7)
% plot_rat_bern_poly([C1;C2],2,.001,"k.")
% hold on
% plot_rat_bern_poly(shapeObject,2,.001,"b.")
% figure(8)
% plot_rat_bern_poly([C1;C2],2,.001,"k.")
% hold on
% plot_rat_bern_poly(shapeObject,2,.001,"b.")

for jjj=1:5
    
    SO=generateTestFigures(jjj);
    if jjj==4
        M=3;
    else
        M=2;
    end
    componentsizes=cellfun(@length,SO);
%     SO{1}=SOS{2};
%     shapeObject=SO;
%     shapeObject{1}=shapeObject{1}+1e-3*rand(size(shapeObject{1}));
%     shapeObject{2}=shapeObject{2}+1e-5*rand(size(shapeObject{2}));
RationalOn=0;
numIntegrands=length(integrands);
elemSize=size(shapeObject,1)/3;
int2evals={};
int2errs={};
int2evals(1:length(integrands))={zeros(length(numIntegrands),1)};
int2errs(1:length(integrands))={zeros(length(numIntegrands),1)};
ggevals=int2evals;
sgevals=int2evals;
global evalCounter;
global scatterEvals;
scatterEvals=0;
for i=8:8%length(integrands)
%     figure(1)
    field = @(x,y) field2(x,y,integrands{i});
    if shape==0
        truev = integral2(field, -1+lroffset(1),1+lroffset(1), @(x)-sqrt(1-(x-lroffset(1)).^2)+lroffset(2), @(x)sqrt(1-(x-lroffset(1)).^2)+lroffset(2),'AbsTol',1e-17,'RelTol',1e-18);
    elseif shape==1
        truev = integral2(field, -sqrt(1-offset^2)+lroffset(1), sqrt(1-offset^2)+lroffset(1), @(x)-sqrt(1-(x-lroffset(1)).^2)+offset+lroffset(2),@(x)sqrt(1-(x-lroffset(1)).^2)-offset+lroffset(2),'AbsTol',1e-17,'RelTol',1e-18);
    elseif shape==2
%         fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,15);
%         mfanti=@(a,b) arrayfun(fanti,a,b);
        SO{1}=shapeObject; truev=ggPolygonIntegrate(SO,field,15,15);
    elseif shape==2
    fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,15);
    mfanti=@(a,b) arrayfun(fanti,a,b);
    SO{1}=shapeObject; truev=ggPolygonIntegrate(SO,field,40,40);
    elseif shape==3
    fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,40);
    mfanti=@(a,b) arrayfun(fanti,a,b);
     truev=ggPolygonIntegrate(SO,field,40,40);
    end
    for jj=1:14
        evalCounter=0;
        if jj==100000
            scatterEvals=1;
        else
            scatterEvals=0;
        end
        if shape==0
            int2errs{i}(jj) = truev-integral2(field, -1+lroffset(1),1+lroffset(1), @(x)-sqrt(1-(x-lroffset(1)).^2)+lroffset(2), @(x)sqrt(1-(x-lroffset(1)).^2)+lroffset(2),'RelTol',10^(-jj));
            int2evals{i}(jj)=evalCounter;
        elseif shape==1
            int2errs(jj) = truev-integral2(field, -sqrt(1-offset^2)+lroffset(1), sqrt(1-offset^2)+lroffset(1), @(x)-sqrt(1-(x-lroffset(1)).^2)+offset+lroffset(2),@(x)sqrt(1-(x-lroffset(1)).^2)-offset+lroffset(2),'RelTol',10^(-jj+3));
            int2evals(jj)=evalCounter;
        end
    end
%     int2errs(int2errs==0)=1e-17;
%     figure(1)
%     semilogy(int2evals,abs(int2errs),'k.','MarkerSize',36)
%     hold on
%     plot_rat_bern_poly(shapeObject,2,.001,'k');
    % Intersect each element, store moment of each material
%     fanti=@(a,b) gauss1D(@(x)field(x,b),0,a,7);
%     mfanti=@(a,b) arrayfun(fanti,a,b);
%     SO{1}=shapeObject; 

%     s=plot_bern_poly(Intersection{1},2,.001,{},{},true);
    %         plot_rat_bern_poly(Intersection{1},2,.1,'r');
    RationalOn=1;
    ggevals=zeros(length(gaussOrders),1);
    ggerrs=zeros(length(gaussOrders),1);
    sgevals=zeros(length(gaussOrders),1);
    sgerrs=zeros(length(gaussOrders),1);
    for jj=8:length(gaussOrders)
        j=gaussOrders(jj);
        if testIntegrands==0
            kk=M*ceil(invTri(i-1))+M*(3);
            kg=ceil(invTri(i-1)+1);
            if i==7
                kk=M*ceil(invTri(i-2))+M*(3);
                kg=ceil(invTri(i-2)+1);
            elseif i==8
                kk=M*2+M*3;
                kg=j;
            end
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
        evalCounter=0;
        ggerrs(jj)=ggPolygonIntegrate(SO,field,jj,kg)-truev;
        ggevals(jj)=evalCounter;
        evalCounter=0;
%         if jj==9
%             scatterEvals=1;
%         else
%             scatterEvals=0;
%         end
        sgerrs(jj)=sgPolygonIntegrate(SO,field,jj,kk,kg)-truev;
        scatterEvals=0;
        sgevals(jj)=evalCounter;
        evalCounter=0;
    end
%     sgerrs;
%     ggerrs;
%     mpr=find(abs(sgerrs)<(3*min(abs(sgerrs(abs(sgerrs)>0)))),1);
%     mpg=find(abs(ggerrs)<(3*min(abs(sgerrs(abs(sgerrs)>0)))),1);
%     monfunctmat(i,jjj)=[ggerrs(find(abs(ggevals)>sgevals(mpr),1))/abs(sgerrs(mpr))];
%     monfunctquad(i,jjj)=sgevals(mpr)/ggevals(mpg);
%     figure(i)
    ggev{i,jjj}=ggevals;
    gger{i,jjj}=ggerrs;
%     semilogy(ggevals,abs(ggerrs),'bx')
%     hold on
    Quadpoints=(invTri(i-1)*kg*+3*invTri(i-1))*M*(sum(componentsizes)/3);
%     Quadpoints=(kk+1)*kg*(sum(componentsizes)/3);
    sgindex=find(sgevals>Quadpoints,1);
%     semilogy(sgevals(sgindex),abs(sgerrs(sgindex)),'ro')
    if i==7
        sger
    end
    if i<8
        sgev{i,jjj}=sgevals(1);
        if ~isempty(sgindex)
            sger{i,jjj}=sgerrs(sgindex);
        else
            sger{i,jjj}=sgerrs(1);
        end
    else
        sgev{i,jjj}=sgevals;%(sgindex);
        sger{i,jjj}=sgerrs;%(sgindex);
    end
%     xlim([0,max(max(ggevals),min(int2evals))])
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
% figure;
% if isempty(find(monfunctmat(:,1)==0))
%     lasti=size(monfunctmat,1)
% else
%     lasti=find(monfunctmat(:,1)==0);
% end
% for jjj=1:1
%     for i=1:invTri(lasti)
%         ertemp=monfunctmat(find(invTri(0:length(monfuncts))==i-1),1);
%         er(i)=mean(ertemp(ertemp~=0 & ertemp~=Inf & ~isnan(ertemp)));
%     end
% plot(0:(invTri(lasti)-1),log(abs(er))/log(10));
% end
% for i=1:invTri(lasti)
%     er(i)=mean(monfunctquad(find(invTri(0:length(monfuncts))==i-1),1));
% end
% plot(0:(invTri(lasti)-1),er)
figure('Position',[10 10 500 700])
for jjj=1:5
    ggert=gger(:,jjj);
    ggevt=ggev(:,jjj);
    sgert=sger(:,jjj);
    sgevt=sgev(:,jjj);
%     for i=1:3
%         subplot(5,5,(jjj-1)*5 +i+1)
%         if jjj==1
%             if i==1
%                 title({"$0$-th degree","monomial"},'Interpreter','Latex','FontName','times')
%             elseif i==2
%                 title({"$1$-st degree","monomials"},'Interpreter','Latex','FontName','times')
%             elseif i==3
%                 title({"$2$-nd degree","monomials"},'Interpreter','Latex','FontName','times')
%             end
%         end    
%         set(gca,'YTickLabel',[]);
%         yyaxis('right')
%         set(gca,'YColor','black','FontName','times');
%         gMer=mean(abs(cat(3,ggert{invTri(0:5)==i-1})),3);
%         gMev=mean(abs(cat(3,ggevt{invTri(0:5)==i-1})),3);
%         sMer=mean(abs(cat(3,sgert{invTri(0:5)==i-1})),3);
%         sMev=mean(abs(cat(3,sgevt{invTri(0:5)==i-1})),3);
%         sMer(sMer==0)=1e-15;
%         semilogy(gMev,gMer,'bx','MarkerSize',4)
%         hold on
%         semilogy(sMev,sMer,'r.','MarkerSize',15)
% %         if i==1 && j==3
% %              ylabel("Absolute Error",'Interpreter','Latex','FontName','times','FontSize',16)
% %         end
%         if i==2 && jjj==5
%             xlabel({"Number of quadrature points"},'Interpreter','Latex','FontName','times','FontSize',16)
%         end
%         axis square
%         ylim([1e-17 1])
%         yticks([1e-15 1e-10 1e-5 1])
%         if i<4
%             set(gca,'YTickLabel',[]);
%         end
%         xtop=max(gMev);
%         xticks([0 floor(xtop/2) xtop]);
%     end
%     subplot(5,5,(jjj-1)*5 +5)
%     
%     if jjj==1
%         title("$p_1(x,y)$",'Interpreter','Latex','FontName','times')
%     end
%     set(gca,'YTickLabel',[]);
%     yyaxis('right')
%     set(gca,'YColor','black','FontName','times');
%     gMer=abs(ggert{7});
%     gMev=ggevt{7};
%     sMer=abs(sgert{7});
%     sMer(sMer==0)=1e-15;
%     sMev=sgevt{7};
%     semilogy(gMev,gMer,'bx','MarkerSize',4)
%     hold on
%     semilogy(sMev,sMer,'r.','MarkerSize',15)
%     axis square
%         ylim([1e-17 1])
%         yticks([1e-15 1e-10 1e-5 1])
% %     set(gca,'YTickLabel',[]);
%     xtop=max(gMev);
%     xticks([0 floor(xtop/2) xtop]);
%     yyaxis('right')
%     set(gca,'YColor','black','FontName','times');
%     if jjj==3
%         ylabel("Absolute Error",'Interpreter','Latex','FontName','times')
%     end
    subplot(5,2,(jjj-1)*2 +2)
    if jjj==1
        title("$f_2(x,y)$",'Interpreter','Latex','FontName','times')
    end
    
    gMer=abs(ggert{8});
    gMev=ggevt{8};
    sMer=abs(sgert{8});
    sMev=sgevt{8};
    sMer(sMer==1e-17)=1e-15;
    semilogy(gMev,gMer,'bx','MarkerSize',4)
    hold on
    semilogy(sMev,sMer,'r.','MarkerSize',15)
    axis square
    yyaxis('right')
    set(gca,'yscale','log')
    ylim([1e-17 1])
    yticks([1e-15 1e-10 1e-5 1])
%     set(gca,'YTickLabel',[]);
    set(gca,'YColor','black','FontName','times');
    if jjj==3
        ylabel("Absolute Error",'Interpreter','Latex','FontName','times','FontSize',16)
    end
    if jjj==5
        xlabel({"Number of ","Integration Points"},'Interpreter','Latex','FontName','times','FontSize',12)
    end
    yyaxis('left')
    set(gca,'YTickLabel',[]);
    xtop=max(gMev);
    xlim([0 2*xtop/5])
    xticks([0 floor(xtop/5) floor(2*xtop/5)]);
end
for jjj=1:5
    subplot(5,2,(jjj-1)*2 +1)
    ss=generateTestFigures(jjj);
    for ii=1:length(ss)
        plot_rat_bern_poly(ss{ii},2,.01,{})
    end
    ll=get(gca,'Children')
    for i=1:length(ll)
        set(ll(i),'SizeData',5)
    end
    if jjj==3 || jjj==2
        axis equal
    else
        axis ij
        axis square
    end
    axis off
    if jjj==1
        title("Region",'Interpreter','Latex')
    end
end
% subplot(5,4,12)
% ylabel("Absolute Error",'FontSize',16,'Interpreter','Latex','FontName','times')
% subplot(5,4,19)
% xlabel("Number of Quadrature Points",'FontSize',16,'Interpreter','Latex','FontName','times')
% figure(7)
% p=3

figure('Position',[10 10 500 700])
for jjj=1:1
    ggert=gger(:,8);
    ggevt=ggev(:,8);
    sgert=sger(:,8);
    sgevt=sgev(:,8);
    for i=1:6
        subplot(3,2,i)
        if jjj==1
            if i==1
                title("$0$th Moment",'Interpreter','Latex','FontName','times')
            elseif i==2
                title("$1$st Moments",'Interpreter','Latex','FontName','times')
            elseif i==3
                title("$2$nd Moments",'Interpreter','Latex','FontName','times')
            elseif i==4
                title("$3$st Moments",'Interpreter','Latex','FontName','times')
            elseif i==5
                title("$4$th Moments",'Interpreter','Latex','FontName','times')
            elseif i==6
                title("$5$th Moments",'Interpreter','Latex','FontName','times')
            end
        end    
        set(gca,'YTickLabel',[]);
        yyaxis('right')
        set(gca,'YColor','black','FontName','times');
        gMer=mean(abs(cat(3,ggert{invTri(0:21)==i-1})),3);
        gMev=mean(abs(cat(3,ggevt{invTri(0:21)==i-1})),3);
        sMer=mean(abs(cat(3,sgert{invTri(0:21)==i-1})),3);
        sMev=mean(abs(cat(3,sgevt{invTri(0:21)==i-1})),3);
        iMer=mean(abs(cat(3,int2errs{invTri(0:21)==i-1})),3);
        iMev=mean(abs(cat(3,int2evals{invTri(0:21)==i-1})),3);
        sMer(sMer==0)=2e-16;
        semilogy(gMev,gMer,'bx','MarkerSize',15)
        hold on
        semilogy(sMev,sMer,'r.','MarkerSize',35)
        semilogy(iMev,iMer,'ko','MarkerSize',15)
%         if i==1 && j==3
%             ylabel("Absolute Error",'Interpreter','Latex','FontName','times')
%         end
%         if i==2 && j==5
%             xlabel({"Number of quadrature points"},'Interpreter','Latex','FontName','times')
%         end
        axis square
        ylim([1e-16 1])
        yticks([1e-15 1e-10 1e-5 1])
        if i==1 | i==3 | i==5
            set(gca,'YTickLabel',[]);
        end
        xtop=floor(1.1*max(min(iMev),min(sMev)));
        xlim([0 xtop]);
        xticks([0 floor(xtop/2) xtop]);
    end
end
subplot(3,2,4)
ylabel("Absolute Error",'FontSize',16,'Interpreter','Latex','FontName','times')
subplot(3,2,5)
xlabel("         Number of Quadrature Points",'FontSize',16,'Interpreter','Latex','FontName','times')
figure(7)
p=3


for ii=1:2
    for i=(1:size(shapeObject{ii},1)/3)
%         plotNURBScurve(2,p,[zeros(1, p+1) ones(1, p+1)],shapeObject{ii}(((i-1)*3+1):(3*i-1),:),shapeObject{ii}((3*i),:))
        plot_rat_bern_poly(shapeObject{ii}(((i-1)*3+1):(3*i),:),2,.001,{})
    end
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

