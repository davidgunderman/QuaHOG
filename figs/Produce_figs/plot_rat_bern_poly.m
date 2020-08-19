function valstot = plot_rat_bern_poly(Ccells,d,tstep,color,orientations)

% if nargin>4
% if ~all(~orientations)
%     for jj=1:length(Ccells)
%         C1=Ccells{jj};
%         n = size(C1,1)/(d+1);
%         valstot=zeros(0,2);
%         for i=1:n
%             vals=dCR_eval(C1((3*(i-1)+1):(3*i),:),0:tstep:1);
%             valstot=[valstot;vals];
%         end
%         if orientations(jj)==1
%             fill(valstot(:,1),valstot(:,2),[0 0 .5],'LineStyle','none','FaceAlpha',.2);
%             hold on
%         else
%             fill(valstot(:,1),valstot(:,2),[1 1 1],'LineStyle','none','FaceAlpha',1);
%             hold on
%         end
%     end
% end
% end
ctr=0;
for jj=1:length(Ccells)
    C1=Ccells{jj};
    n = size(C1,1)/(d+1);
    for i=1:1:n
        vals=dCR_eval(C1(((d+1)*(i-1)+1):((d+1)*i),:),0:tstep:1);
        if isempty(color)
            if d==2
            plot(vals(:,1),vals(:,2),'LineWidth',2);
            elseif d==3
               plot3(vals(:,1),vals(:,2),vals(:,3),'LineWidth',1);
            end
%             text(vals(floor(length(vals)/2),1),vals(floor(length(vals)/2),2),...
%                 sprintf("$c_{%d}$",ctr+1),'Interpreter','Latex','FontName',...
%                 'times','FontSize',12);
%             ctr=ctr+1;
        else
            if d==2
            plot(vals(:,1),vals(:,2),'LineWidth',2,'Color',color);
            elseif d==3
           plot3(vals(:,1),vals(:,2),vals(:,3),'LineWidth',1,'Color',color);
            end
        end
        hold on
    end
end
end