clear vars; close all;
% Tests shapeIEval on a circle. Plots the function f(x,y) inside a circle 0
% outside

global evalCounter;
f = @(x,y) ones(size(x));
% Circletemp=load("PolyCircle.mat"); Circle= Circletemp.PolyCircle;
% plot_bern_poly(Circle,2,.01,{},{},false);
Circletemp=load("Circle.mat"); Circle= Circletemp.Circle1;

Ifun = @(x,y) shapeIEval([x y],f,Circle,1);
arrayIfun = @(x,y) arrayfun(Ifun,x,y);
for i=1:2
    figure;
    plot_rat_bern_poly(Circle,2,.01,[0 0 0]);
    hold on
    evalCounter=0;
    tic
    int2val(i)=integral2(arrayIfun,-1.0110105618234,1.017296712386571,-1.016727386,1.065823757183,'RelTol',2^(-i));
    numEvals(i)=evalCounter;
    timings(i)=toc;
end

field = @(x,y) field2(x,y,f);
jj=1; kk=14; kg=2; SO{1}=Circle;
sgval=sgPolygonIntegrate(SO,field,jj,kk,kg);