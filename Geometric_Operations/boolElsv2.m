function Cellintel = boolElsv2(Element1,Element2,InorUn)

%% Important constants defining size of Elements
    d=2;
    numc1=size(Element1,1)/d;
    numc2=size(Element2,1)/d;
    n=size(Element1,2);
    n2=size(Element2,2);
    
    
    intlist=zeros(0,5); % To store all of the intersections
    numinters=0;
    for i=1:numc1
        for j=1:numc2
            if did_intersect_box(Element1(((d*(i-1))+1):((d*i)),:),Element2(((d*(j-1))+1):((d*j)),:))
                [~,~,temps,tempt]=bern_inter_recurse_v2(Element1(((d*(i-1))+1):((d*i)),:),Element2(((d*(j-1))+1):((d*j)),:),[0 1], [0 1], 10^(-22));
%                 temps = unique(temps)
                if ~isempty(temps)
%                     scatter(dC_eval(Element1((d*(i-1)+1),:),temps(:,2)),dC_eval(Element1((d*i),:),temps(:,2)),100,"r.");
                    intlist((numinters+1):(numinters+length(temps(:,2))),1:4)=[temps(:,2), tempt(:,2), repmat(i,length(temps(:,2)),1), repmat(j,length(temps(:,2)),1)];
                    numinters=length(temps(:,2))+numinters;
                end
            end
        end
    end
    size(intlist);
    inclusion=0;
%     if isempty(intlist)

%     end
    constintlist=intlist; % This is a list of all of the intersections
    addedintel=0; % This is the number of segments on the current region of intersections
    
    % These variables define which intersection to start at
    startint=1;
    if numinters>0
        startelem=orient(...
        Element1(((d*(intlist(startint,3)-1))+1):((d*intlist(startint,3))),:),...
        Element2(((d*(intlist(startint,4)-1))+1):((d*intlist(startint,4))),:),...
        intlist(startint,1),...
        intlist(startint,2)...
        );
        if InorUn
            startelem=~startelem;
        end
        startc=intlist(startint,startelem+3);
        starts=intlist(startint,startelem+1); 
        %     if startelem Element1(
        % These variables change during the walking process
        element = startelem; % This is the current polygon
        k=startc; % This is the number of the current curve
        temps=starts; % This is the time on the current curve
        startnuminters=numinters;% This is the total number of intersections between the two polygons
        intlist(:,5)=0; % None of the intersections have been used yet
        ll=1; % This is the number of regions that are part of the intersection
        addingcurves=true; % This is a boolean to conrol whether curves are 
                         % being added or we are travelling between regions

    
    %% This is the directional walking method
    % Keep going until we have gone through all of the intersections and returned to an intersection
    while numinters>0
        % Find the intersections that haven't been used yet along the 
        % current curve in order of increasing parameter value
        [~,idxsort1]=sort(intlist(:,1+element));
        intlist=intlist(idxsort1,:);
        finder=find((~intlist(:,5) & intlist(:,3+element)==k & intlist(:,1+element)>temps),1);

        
            % If there are intersections, use them to add the next curve
        if ~isempty(finder) 
            if ~element
                Splitter1=Element1(((d*(k-1))+1):((d*k)),:);
            else
                Splitter1=Element2(((d*(k-1))+1):((d*k)),:);
            end
            splitpoints=[temps,intlist(finder,1+element)];
            splitpoints(2)=((splitpoints(2)-splitpoints(1))./(1-splitpoints(1)));
            for i=1:2
                Splitter1((d*(i-1)+1):(d*(i+1)),:)=dC_split(Splitter1((d*(i-1)+1):(d*(i)),:),splitpoints(i));
            end        
            % Check to see if we've made it back to the starting intersection of
            % this region. If so, then store the region, stay on the same curve, 
            % and get to the next intersection point.
            if intlist(finder,3+~element)==startc && ~element==startelem && numinters<startnuminters && intlist(finder,1+~element)==starts
                intlist(finder,5)=1; % Mark the current intersection as having been used
                temps=intlist(finder,1+element);
                numinters=numinters-1;
                if ~element
                    Intel((2*addedintel+1):(2*addedintel+(d)),:)=[Splitter1((d+1):(2*d),:),repmat(Splitter1((d+1):(2*d),end),1,max([n n2])-n)]; % Intersection hit
                else
                    Intel((2*addedintel+1):(2*addedintel+(d)),:)=[Splitter1((d+1):(2*d),:),repmat(Splitter1((d+1):(2*d),end),1,max([n n2])-n2)]; % Intersection hit
                end
                Cellintel{ll}=Intel; % Start a new region and
                ll=ll+1;             % add one to the region counter
                Intel=zeros(0,0);    % Create an object to store the current region
                addedintel=0;        % Reset the counter for the number of curves on this region
                startelem=element;   % Reset the start element to be this one
                startnuminters=numinters; % Reset the remaining intersections to be the current amount
    %             intlist(idxsort1,:)=intlist; % unsort the intersection list.
                addingcurves =false; % don't add curves until we get to the next intersection
            elseif addingcurves
                intlist(finder,5)=1; % Mark the current intersection as having been used
                element=~element;  %Switch elements
                k=intlist(finder,3+element);  % Set the curve# to where the inters is on the other quad
                temps=intlist(finder,1+element);
                numinters=numinters-1;
                if element
                    Intel((2*addedintel+1):(2*addedintel+(d)),:)=[Splitter1((d+1):(2*d),:),repmat(Splitter1((d+1):(2*d),end),1,max([n n2])-n)]; % Intersection hit
                else
                    Intel((2*addedintel+1):(2*addedintel+(d)),:)=[Splitter1((d+1):(2*d),:),repmat(Splitter1((d+1):(2*d),end),1,max([n n2])-n2)]; % Intersection hit
                end
%                 Intel((2*addedintel+1):(2*addedintel+(d)),:)=Splitter1((d+1):(2*d),:); % Intersection hit
                addedintel=addedintel+1; %We've added a curve to the intersection polygon
            else
                temps=intlist(finder,1+element);
                starts=intlist(finder,1+element);
                addingcurves=true;
                startc=intlist(finder,3+element);
                startelem=element;
            end
            intlist(idxsort1,:)=intlist;
            
        % If there are no intersections, add to the end of the current curve
        elseif element 
            if addingcurves
                Splitter1=dC_split(Element2(((d*(k-1))+1):((d*k)),:),temps);
                Intel((2*addedintel+1):(2*addedintel+(d)),:)=[Splitter1((d+1):(2*d),:),repmat(Splitter1((d+1):(2*d),end),1,max([n n2])-n2)]; % Intersection hit
%                 Intel((2*addedintel+1):(2*addedintel+(d)),:)=Splitter1((d+1):(2*d),:);
                addedintel=addedintel+1; %We've added a curve
            end
            k=mod(k,numc2)+1;
            temps=0; %We're at the start of the next curve on element 2
            intlist(idxsort1,:)=intlist;
        else
            if addingcurves
                Splitter1=dC_split(Element1(((d*(k-1))+1):((d*k)),:),temps);
                Intel((2*addedintel+1):(2*addedintel+(d)),:)=[Splitter1((d+1):(2*d),:),repmat(Splitter1((d+1):(2*d),end),1,max([n n2])-n)]; % Intersection hit
%               Intel((2*addedintel+1):(2*addedintel+(d)),:)=Splitter1((d+1):(2*d),:);
                addedintel=addedintel+1; %We've added a curve
            end
            k=mod(k,numc2)+1;
            temps=0; %We're at the start of the next curve on element 1
            intlist(idxsort1,:)=intlist;
        end
        % Debugging
%         plot_bern_poly(Intel,2,.0001,"");
%         [element,k]
%         Intel
    end
    if ~isempty(Intel)
        Cellintel{ll}=Intel;
    end
    intlist=constintlist;
 
    else
        inclusion=isContained(Element1,Element2);
        if inclusion==1
            Cellintel{1}=Element1;
        elseif inclusion==2
            Cellintel{1}=Element2;
        else
            Cellintel{1}=zeros(0,size(Element1,2));
        end
    end
end

function orientation = orient(Element1,Element2,s,t)
crossp= cross([dC_eval_dt(Element1,s) 0],[dC_eval_dt(Element2,t) 0]);
orientation = crossp(3)<0;
% c1vals=dC_eval(Element1,[s+10^(-13) s-10^(-13)]);
% c2vals=dC_eval(Element2,t+10^(-13));     
% % scatter(dC_eval(Element1(1,:),s),dC_eval(Element1(2,:),s),500,'k.');
% % hold on
% orientation=(-(c1vals(1,1)-c1vals(2,1))*(c2vals(2)-c1vals(2,2)) + (c1vals(1,2)-c1vals(2,2))*(c2vals(1)-c1vals(2,1)))>0;
% % Orientation is found using scaled signed distance from point to line equation:
% % D* = (x2-x1)(y-y1) + (y2-y1)(x-x1) 
end

function [oC1,oC2,otbds,osbds] = bern_inter_recurse_v2(C1,C2,itbds,isbds,tol)
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
[d,n]=size(C1);
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
% p1_scale=1;
% p2_scale=1;
otbds=zeros(0,2);
osbds=zeros(0,2);
oC1=zeros(0,size(C1,2));
oC2=zeros(0,size(C1,2));
if (d1>tol*(p1_scale))||(d2>tol*(p2_scale))
    c1s=dC_split(C1,.5);
    c2s=dC_split(C2,.5);
    intersect_bool_mat=[did_intersect_box(c1s(1:d,:),c2s(1:d,:)) did_intersect_box(c1s((d+1):(2*d),:),c2s(1:d,:));
    did_intersect_box(c1s(1:d,:),c2s((d+1):(2*d),:)) did_intersect_box(c1s((d+1):(2*d),:),c2s((d+1):(2*d),:));];
    if intersect_bool_mat(1,1)
        [temp1,temp2,tempt,temps]=bern_inter_recurse_v2(c1s(1:d,:),c2s(1:d,:), ...
            [itbds(1),mean(itbds)],[isbds(1),mean(isbds)],tol);
        oC1=[oC1;temp1]; oC2=[oC2;temp2];
        otbds=[otbds;tempt]; osbds=[osbds;temps];
    end
    if intersect_bool_mat(1,2)
        [temp1,temp2,tempt,temps]=bern_inter_recurse_v2(c1s((d+1):(2*d),:),c2s(1:d,:), ...
            [mean(itbds),itbds(2)],[isbds(1),mean(isbds)],tol);
        oC1=[oC1;temp1]; oC2=[oC2;temp2];
        otbds=[otbds;tempt]; osbds=[osbds;temps];
    end
    if intersect_bool_mat(2,1)
        [temp1,temp2,tempt,temps]=bern_inter_recurse_v2(c1s(1:d,:),c2s((d+1):(2*d),:), ...
            [itbds(1),mean(itbds)],[mean(isbds),isbds(2)],tol);
        oC1=[oC1;temp1]; oC2=[oC2;temp2];
        otbds=[otbds;tempt]; osbds=[osbds;temps];
    end
    if intersect_bool_mat(2,2)
        [temp1,temp2,tempt,temps]=bern_inter_recurse_v2(c1s((d+1):(2*d),:),c2s((d+1):(2*d),:), ...
            [mean(itbds),itbds(2)],[mean(isbds),isbds(2)],tol);
        oC1=[oC1;temp1]; oC2=[oC2;temp2];
        otbds=[otbds;tempt]; osbds=[osbds;temps];
    end
else 
    p=newton_intersect(C1(:,1),C1(:,end),C2(:,end),C2(:,1));
    if all(p<(1) & p>0) && all(~isnan(p))
        oC1=C1;
        oC2=C2;
        otbds=[1 itbds(1)+(itbds(2)-itbds(1))*p(1)];
        osbds=[1 isbds(1)+(isbds(2)-isbds(1))*p(2)];
    end
end
end

function d = point_to_line(pt, v1, v2)
% pt should be nx3
% v1 and v2 are vertices on the line (each 1x3)
% d is a nx1 vector with the orthogonal distances
v1 = repmat(v1,size(pt,1),1);
v2 = repmat(v2,size(pt,1),1);
a = v1 - v2;
b = pt - v2;
% d = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));
d = sum(cross2d(a,b).^2,2) ./ sum(a.^2,2);
end

function cr = cross2d(a,b)
cr = a(:,1).*b(:,2)-a(:,2).*b(:,1);
end

function t = newton_intersect(a,d,c,b)
% Gives parameter values of intersection of two lines: a+(d-a)t and
% c+(b-c)s for 0<t<1 and 0<s<1.
% Input parameters:
% a,d,b,c: d by 1 vectors representing start/end points of each line
% output parameters:
% t : d by 1 vector representing "time" of interesction for each line
% determ = a(2)*b(1)-a(1)*b(2)-a(2)*c(1)+a(1)*c(2)+b(2)*d(1)-c(2)*d(1)-b(1)*d(2)+c(1)*d(2);
% Test the bad line intersection algorithm currently used in c++
% t=zeros(2,1);
% t(1)=(1.0/determ)*((b(2)-c(2))*(b(1)-a(1))+(c(1)-b(1))*(b(2)-a(2)));
% t(2)=(1.0/determ)*((a(2)-d(2))*(b(1)-a(1))+(d(1)-a(1))*(b(2)-a(2)));
% if determ==0
%     t=1.5*ones(2,1);
% end

[t]=[d-a,b-c]\(b-a);
% if any(t==Inf | t==-Inf | isnan(t))
%     t=[.5,.5];
% end
end

% This functions assumes that the two elements have NO intersections !!!!!
% It returns 0 if the two elements are totally separated
% It returns 1 if element 1 is totally contained within element 2
% It returns 2 if element 2 is totally contained within element 1
% Returns 1 if element1 is totally contained within element 2
% Returns 0 if elements are separate
function whichContained = isContained(Element1, Element2)
    d=2;
    numc1=size(Element1,1)/d;
    numc2=size(Element2,1)/d;
    n=size(Element1,2);
    n2=size(Element2,2);
    LineGuess=[linspace(dC_eval(Element1(1,:),.5) , dC_eval(Element2(1,:),.5),n); linspace(dC_eval(Element1(2,:),.5)  , dC_eval(Element2(2,:),.5),n)];
    Linestart=[1 .5];
    Lineend=[1 .5];
    line1s=0;
    line2s=1;
for j=1:numc1
    [~,~,temps,tempt]=bern_inter_recurse_v2(LineGuess,Element1((d*(j-1)+1):(d*j),:),[0 1], [0 1], 10^(-8));
    for i=1:size(temps,1)
        if temps(i,2)>line1s
            line1s=temps(i,2);
            Linestart=[j tempt(i,2)];
        end
    end
end

for j=1:numc2
    [~,~,temps,tempt]=bern_inter_recurse_v2(LineGuess,Element2((d*(j-1)+1):(d*j),:),[0 1], [0 1], 10^(-8));
    for i=1:size(temps,1)
        if (temps(i,2)<line2s)&& (temps(i,2)>line1s)
            line2s=temps(i,2);
            Lineend=[j tempt(i,2)];
        end
    end
end
Linedt=dC_eval_dt(LineGuess,line2s);
Element1dt=dC_eval_dt(Element1((d*Linestart(1)-1):(d*Linestart(1)),:),Linestart(2));
Element2dt=dC_eval_dt(Element2((d*Lineend(1)-1):(d*Lineend(1)),:),Lineend(2));

E1inE2=cross2d(Element1dt,Linedt)>0;
E2inE1=cross2d(Element2dt,Linedt)<0;
if E1inE2 && ~E2inE1
    whichContained=2;
%     plot_bern_poly(Element1,2,.001,"",{"b.","b.","b.","b."},true)
elseif ~E1inE2 && E2inE1
    whichContained=1;
%     plot_bern_poly(Element2,2,.001,"",{"b.","b.","b.","b."},true)
else
    whichContained=0;
%     plot_bern_poly(Element1,2,.001,"",{"r.","r.","r.","r."},true)
%     plot_bern_poly(Element2,2,.001,"",{"r.","r.","r.","r."},true)
end
end