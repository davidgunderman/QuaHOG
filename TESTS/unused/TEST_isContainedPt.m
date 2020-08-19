clear vars; close all;

% Circletemp=load("PolyCircle.mat"); Circle= Circletemp.PolyCircle;
% plot_bern_poly(Circle,2,.01,{},{},false);
Circletemp=load("Circle.mat"); Circle= Circletemp.Circle1;
plot_rat_bern_poly(Circle,2,.01,[0 0 0]);
hold on
[PTSx,PTSy]=meshgrid(-1:.11:1,-1:.11:1);
PTSx=PTSx(:); PTSy=PTSy(:);
scatter(PTSx,PTSy,40,'k.');
inn(1:length(PTSx))=false;
for i=1:length(PTSx)
    inn(i)=isContainedPt([PTSx(i), PTSy(i)],Circle,1);
end
scatter(PTSx(inn),PTSy(inn),40,'r.')
    