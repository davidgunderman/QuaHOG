function hh = visBezElTrim(Element,curves)
NEls=length(Element);
oldC=Element;
for j=1:NEls
curves{j}=correctOrder(curves{j});
st2=zeros(0,0);
for i=1:3:size(curves{j},1)
    [st]=dCR_eval(curves{j}((i:(i+2)),:),0:.05:1);
    [st]=dCR_eval(curves{j}((i:(i+2)),:),0:.05:1);
    st2=[st2;st(1:(end-1),:)];
    xyzcur=dCR_eval_3d(Element{j},st);
    plot3(xyzcur(:,1),xyzcur(:,2),xyzcur(:,3),'k');
    hold on
end
    options.output=false;
    options.mlim=.02;
    options.maxit=20;
    options.dhmax=.3;
    [p,t]=mesh2d(st2,[],[],options);
    xyzt=dCR_eval_3d(Element{1},p);
    trimesh(t,xyzt(:,1),xyzt(:,2),xyzt(:,3),'FaceColor','b','LineStyle','none','FaceAlpha',.1,'Edgecolor','none')

end
% clear xyzcur;
% for j=1:NEls
%     pts=zeros(0,0);
%     xycur=zeros(0,0);
% for i=1:(size(curves{j},1)/3)
%     xyzcur((4*(i-1)+1):(4*i),:)=[dCR_eval_3d(Element{j},curves{j}((3*(i-1)+1):(3*i-1),:)')'; curves{j}(3*i,:)];
%     pts=[pts;dCR_eval(xyzcur((4*(i-1)+1):(4*i),:),0:.05:1)];
% end
%     fill3(pts(:,1),pts(:,2),pts(:,3),'b');
%     hold on
% end
end