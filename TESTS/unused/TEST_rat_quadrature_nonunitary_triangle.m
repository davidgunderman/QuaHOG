close
clearvars


               
for jj=1:5
    weightrat=2^(1-jj);
    shapeObject=[0 0 0;
           0 .5 1;
           1 1 1;
           0 weightrat/2 weightrat;
           1 weightrat/2 0;
           1 weightrat/2 weightrat;
           weightrat weightrat/2 0;
           0 0 0;
           weightrat weightrat/2 1;];
%        plot_rat_bern_poly(shapeObject,2,.01,'k')
    ctrlpoints(1,:,:)=[0    0  0
                       weightrat/2 weightrat/2 0
                       weightrat 0  0];
    ctrlpoints(2,:,:)=[0 .5   1
                       0  weightrat/2  0
                       0  0   0];
    ctrlpoints(3,:,:)=[1    1     1
                       weightrat/2   weightrat/2    0
                       weightrat 0   0];

    testfunct= @(x,y) x;
    field = @(x,y) field2(x,y,testfunct);
    fieldpt = @(pt) field2(pt(:,1),pt(:,2),testfunct);
    Tfield = @(x,y) fieldpt(TdCReval(ctrlpoints,[x,y,1-x-y],2)).*TdCjacobian(ctrlpoints,[x,y,1-x-y],2);
    
    color=.3+.7*rand(3,1);
    color(floor(3*rand(1)+1))=0;
    global evalCounter;
    for i=1:60
        evalCounter=0;
        SO{1}=shapeObject; gv(i)=ggPolygonIntegrate(SO,field,i,i);
        ge(i)=evalCounter; evalCounter=0;
        tv(i)=gaussT(Tfield,i);
        te(i)=evalCounter; 
%         evalCounter=0;
%         sv(i)=sgPolygonIntegrate(SO,field,0,i,i);
%         se(i)=evalCounter;
    end
    evalCounter=0;
    sv(jj)=sgPolygonIntegrate(SO,field,0,8,1);
    se(jj)=evalCounter;
%     ger=abs(gv+tv(60));
    ter=abs(tv-tv(60));
    ser=abs(sv+tv(60));
%     ger(ger<=1e-16)=1e-16;
    ter(ter<=1e-16)=1e-16;
    ser(ser<=1e-16)=1e-16;
    loglog(te,ter,'color',color');
    hold on
%     evalCounter=0;
    ser(jj)=abs(sv(jj)+tv(60));
%     loglog(se,ser,'--','color',color');
%     scatter(se,ser,100,'r','.');
    scatter(se(jj),ser(jj),100,color','.','HandleVisibility','off');
end
legend({"Weight=$1$","Weight=$.5^1$","Weight=$.5^2$","Weight=$.5^3$","Weight=$.5^4$"},'Interpreter','latex')
