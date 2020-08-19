function [oR1,oR2,otbds,osbds] = Ratbern_inter_recurse_v2(R1,R2,itbds,isbds,tol)
% Given two bernstein polynomials, finds their intersection points.
% Input parameters
% C1: This is a d by n matrix consisting of control points of all possible
% intersection polynomials coming from the first original polynomial
% C1: This is a d by n matrix consisting of control points of all possible
% intersection polynomials coming from the first original polynomial
% Output parameters
% C1: A d by n matrix consisting of control points for the bernstein
% polynomial of degree n of dimension d
% C2: A d by n matrix consisting of control points for th
% p: A list of numerically calculated intersection points (empty if none)
d=size(R1,1)-1;
W1=R1(end,:); W2=R2(end,:);
C1=R1(1:d,:)./W1; C2=R2(1:d,:)./W2;
% Number of polynomials to check
% w=0;
% h=0;
% figure('Units','normalized','Position', [.5+(w) .6+h .0675 .125])
% if w<(.5-.0675)
%     w=w+.0675;
% else
%     w=.5;
%     h=h-.125;
% end
% hold off   
% plot_bern_poly(C1(1,:),C1(2,:),.01);
% plot_bern_poly(C2(1,:),C2(2,:),.01);
% drawnow
% d1=point_to_line(padarray(C1,1,0,'post')',padarray(C1(:,1),1,0,'post')',padarray(C1(:,end),1,0,'post')')
d1=sum(point_to_line(C1',C1(:,1)',C1(:,end)'));
d2=sum(point_to_line(C2',C2(:,1)',C2(:,end)'));
p1_scale=max(abs(C1(:,end)-C1(:,1)));
p2_scale=max(abs(C2(:,end)-C2(:,1)));
otbds=zeros(0,2);
osbds=zeros(0,2);
oR1=zeros(0,size(C1,2));
oR2=zeros(0,size(C1,2));
if (d1>tol)||(d2>tol)
    r1s=dC_split(R1,.5); r2s=dC_split(R2,.5);
    w1s=r1s([3,6],:); w2s=r2s([3,6],:);
    c1s=[r1s(1:2,:)./w1s(1,:);r1s(4:5,:)./w1s(2,:)];
    c2s=[r2s(1:2,:)./w2s(1,:);r2s(4:5,:)./w2s(2,:)];
    intersect_bool_mat=[did_intersect_box(c1s(1:d,:),c2s(1:d,:)) did_intersect_box(c1s((d+1):(2*d),:),c2s(1:d,:));
    did_intersect_box(c1s(1:d,:),c2s((d+1):(2*d),:)) did_intersect_box(c1s((d+1):(2*d),:),c2s((d+1):(2*d),:));];
    if intersect_bool_mat(1,1)
        [temp1,temp2,tempt,temps]=Ratbern_inter_recurse_v2(r1s(1:(d+1),:),r2s(1:(d+1),:), ...
            [itbds(1),mean(itbds)],[isbds(1),mean(isbds)],tol);
        oR1=[oR1;temp1]; oR2=[oR2;temp2];
        otbds=[otbds;tempt]; osbds=[osbds;temps];
    end
    if intersect_bool_mat(1,2)
        [temp1,temp2,tempt,temps]=Ratbern_inter_recurse_v2(r1s((d+2):(2*(d+1)),:),r2s(1:(d+1),:), ...
            [mean(itbds),itbds(2)],[isbds(1),mean(isbds)],tol);
        oR1=[oR1;temp1]; oR2=[oR2;temp2];
        otbds=[otbds;tempt]; osbds=[osbds;temps];
    end
    if intersect_bool_mat(2,1)
        [temp1,temp2,tempt,temps]=Ratbern_inter_recurse_v2(r1s(1:(d+1),:),r2s((d+2):(2*(d+1)),:), ...
            [itbds(1),mean(itbds)],[mean(isbds),isbds(2)],tol);
        oR1=[oR1;temp1]; oR2=[oR2;temp2];
        otbds=[otbds;tempt]; osbds=[osbds;temps];
    end
    if intersect_bool_mat(2,2)
        [temp1,temp2,tempt,temps]=Ratbern_inter_recurse_v2(r1s((d+2):(2*(d+1)),:),r2s((d+2):(2*(d+1)),:), ...
            [mean(itbds),itbds(2)],[mean(isbds),isbds(2)],tol);
        oR1=[oR1;temp1]; oR2=[oR2;temp2];
        otbds=[otbds;tempt]; osbds=[osbds;temps];
    end
else 
    p=newton_intersect(C1(:,1),C1(:,end),C2(:,end),C2(:,1));
    if all(p<(1) & p>=0) && all(~isnan(p))
        oR1=C1;
        oR2=C2;
        otbds=[1 itbds(1)+(itbds(2)-itbds(1))*p(1)];
        osbds=[1 isbds(1)+(isbds(2)-isbds(1))*p(2)];
    end
end
