% Plot a 2-D polygon with d-th order edges using control points in c and
% resolution tstep
function h1 = plot_bern_poly(C1,d,tstep,label,colorarr,filling)
% Input parameters
%   cx      n by (d+1) matrix with (d+1) weights for each of n x-coordinate 
%          bernstein polynomials. The first column are the vertices of the polygon.
%   cxy      n by (d+1) matrix with (d+1) weights for each of n y-coordinate 
%          bernstein polynomials. The first column are the vertices of the polygon.
%   tstep  resolution to use for plotting
cx=C1(1:d:end,:);
cy=C1(2:d:end,:);
cz=C1(3:d:end,:);
n= size(cx,1);
% colorarr=["b.","r.","k.","k."];
% colorarr=["b.","r.","k.","k."];
for i=1:n
    plotx(((i-1)*(1/tstep+1)+1):((i)*(1/tstep+1)))=bernsteinpoly(cx(i,:),0:tstep:1);
    ploty(((i-1)*(1/tstep+1)+1):((i)*(1/tstep+1)))=bernsteinpoly(cy(i,:),0:tstep:1);
    if filling ==false
        if d==2
            if ~isempty(colorarr)
                h1{i}= plot(bernsteinpoly(cx(i,:),0:tstep:1),bernsteinpoly(cy(i,:),0:tstep:1),'Color',colorarr{1});
            else
                h1{i} =plot(bernsteinpoly(cx(i,:),0:tstep:1),bernsteinpoly(cy(i,:),0:tstep:1)); %'k.'
            end
            hold on
       elseif d==1
           scatter(0:tstep:1,bernsteinpoly(cx(i,:),0:tstep:1));
           hold on
       else
            if ~isempty(colorarr)
           plot3(bernsteinpoly(cx(i,:),0:tstep:1),bernsteinpoly(cy(i,:),0:tstep:1),bernsteinpoly(cz(i,:),0:tstep:1),'LineWidth',1,'Color',colorarr{1});
            else
           scatter3(bernsteinpoly(cx(i,:),0:tstep:1),bernsteinpoly(cy(i,:),0:tstep:1),bernsteinpoly(cz(i,:),0:tstep:1),'k.');
            end
           hold on
        end 
    end
end
if ~isempty(label)
    text(plotx(1),ploty(1),label)
end

% text(bernsteinpoly(cx(n,:),.5)-.2*(max(cx(:))-min(cx(:))),bernsteinpoly(cy(n,:),.5)-.1*(max(cy(:))-min(cy(:))),label,'FontSize',12);
if filling
%     for i=1:6
% %         plot([0 cx(i,1)],[0 cy(i,1)],'k','Linewidth',2)
% % %         fill([0:(cx(i,1)/100):cx(i,1),plotx((((i-1)*(1001)+1):((i)*(1001)))),0:(cx(i,end)/100):cx(i,end)],[zeros(1,101),ploty((((i-1)*(1001)+1):((i)*(1001)))),0:(cy(i,end)/100):cy(i,end)],[1 0 0],'FaceAlpha',.25);
% %         fill(plotx((((i-1)*(1001)+1):((i)*(1001)))),ploty((((i-1)*(1001)+1):((i)*(1001)))),[1 0 0],'FaceAlpha',.25);
        h1 = fill(plotx,ploty,[0.2,0.8,0.11],'FaceAlpha',1,'EdgeColor',[0.2,0.8,0.11]);
% %         fill([0 2 2],[0 0 -1.5],[1 0 0],'FaceAlpha',.25);
% %     scatter(cx(:),cy(:),300,'k')
%     end
end
% hold on
% xlim auto
% ylim auto
end
