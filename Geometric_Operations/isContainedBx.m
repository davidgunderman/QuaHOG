function InsideBx = isContainedBx(Bx, Element2,rat,tol,pts)
% Returns 0 if rectangle Bx is inside Element2
% Returns 1 if rectangle Bx intersects Element2
% Returns 2 if rectange Bx is outside Element2
% Inside=zeros(1,4);
% bdpts = [Bx(1) Bx(3);
%         Bx(1) Bx(4);
%         Bx(2) Bx(3);
%         Bx(2) Bx(4);];
 sd{1}=[Bx(1) Bx(1);
        Bx(3) Bx(4);];
    sd{2}=[Bx(2) Bx(2);
        Bx(3) Bx(4);];
    sd{3}=[Bx(1) Bx(2);
        Bx(3) Bx(3);];
    sd{4}=[Bx(1) Bx(2);
        Bx(4) Bx(4);];
% [x, ~] = GaussLegendreQuad1D(pts,Bx(1),Bx(2));
% [y, ~] = GaussLegendreQuad1D(pts,Bx(3),Bx(4));
for i=1:4
    p=size(Element2,2)-1;
    csd=[linspace(sd{i}(1,1),sd{i}(1,2),p+1);
        linspace(sd{i}(2,1),sd{i}(2,2),p+1);
        ones(1,p+1);];
    for k=1:3:size(Element2,1) 
        [~,~,temps,tempt]=Ratbern_inter_recurse_v2(csd,Element2(k:(k+2),:),[0 1], [0 1], 10^(-7));
        if ~isempty(temps)
            break
        end
    end
    if ~isempty(temps)
        break
    end
end
% [s,t]=meshgrid(x,y);
% ss=s(:); tt=t(:);
% ckpts=[ss,tt];
% ckpts=[ckpts;bdpts];
if rat
    ptDom=dCR_eval(Element2(1:3,:),.5);
else
    ptDom=dC_eval(Element2(1:2,:),.5);
end
ptDom=[Bx(1) Bx(3)];
% for pti =1:size(ckpts,1)
%     Inside(pti) = isContainedPt(ckpts(pti,:),Element2,rat,tol);
% end
Inside=isContainedPt(ptDom,Element2,rat,tol);
DomInside=(ptDom(1)>Bx(1) && ptDom(1)<Bx(2) &&...
            ptDom(2)>Bx(3) && ptDom(2)<Bx(4)); % point on Element2 is inside box
if isempty(temps)&& Inside
    InsideBx =0; %Bx inside Element2
elseif isempty(temps) && ~Inside
    InsideBx = 2; %Bx outside Element2
else
    InsideBx = 1; %Bx intersects Element2
end
% scatter(ckpts(:,1),ckpts(:,2),'k.');
% hold on
end