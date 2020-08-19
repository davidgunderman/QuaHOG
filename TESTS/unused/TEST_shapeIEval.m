clear vars; close all;
% Tests shapeIEval on a circle. Plots the function f(x,y) inside a circle 0
% outside

f = @(x,y) 1-x.^2;
% Circletemp=load("PolyCircle.mat"); Circle= Circletemp.PolyCircle;
% plot_bern_poly(Circle,2,.01,{},{},false);
Circletemp=load("Circle.mat"); Circle= Circletemp.Circle1;
% plot_rat_bern_poly(Circle,2,.01,[0 0 0]);
hold on
[PTSx,PTSy]=meshgrid(-1:.02:1,-1:.02:1);
PTSx=PTSx(:); PTSy=PTSy(:);
% scatter(PTSx,PTSy,40,'k.');
evals(1:length(PTSx))=0;
for i=1:length(PTSx)
    evals(i)=shapeIEval([PTSx(i), PTSy(i)],f,Circle,1);
end
nn=sqrt(length(PTSx));
surf(reshape(PTSx,nn,nn),reshape(PTSy,nn,nn),reshape(evals,nn,nn))
    