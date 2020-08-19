function [xq,yq,wq]=quadtreeIntegratePts(CPW_my_format,intlvls,lvls,tol,orientations)
xq=zeros(0,1); yq=zeros(0,1); wq=zeros(0,1);
    for ii=1:length(CPW_my_format)
    p=2;
    bdbx=[min(CPW_my_format{ii}(1:3:end,:)./CPW_my_format{ii}(3:3:end,:),[],'all')-1e-3*rand(1),...
          max(CPW_my_format{ii}(1:3:end,:)./CPW_my_format{ii}(3:3:end,:),[],'all')+1e-3*rand(1),...
          min(CPW_my_format{ii}(2:3:end,:)./CPW_my_format{ii}(3:3:end,:),[],'all')-1e-3*rand(1),...
          max(CPW_my_format{ii}(2:3:end,:)./CPW_my_format{ii}(3:3:end,:),[],'all')+1e-3*rand(1)
          ];
    pts=recurse_quad_pts(CPW_my_format{ii},bdbx,p,intlvls,lvls,tol);
    xq=[xq; pts(:,1)]; yq=[yq; pts(:,2)]; wq=[wq; orientations(ii)*pts(:,3)];
    end
end