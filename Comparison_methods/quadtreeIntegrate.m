function [xq,yq,wq]=quadtreeIntegratePts(CPW_my_format,0,lvls,1e-8)
 p=size(CPW_my_format,2)-1;
            bdbx=[min(CPW_my_format{ii}(1:3:end,:)./CPW_my_format{ii}(3:3:end,:),[],'all')-1e-3*rand(1),...
                  max(CPW_my_format{ii}(1:3:end,:)./CPW_my_format{ii}(3:3:end,:),[],'all')+1e-3*rand(1),...
                  min(CPW_my_format{ii}(2:3:end,:)./CPW_my_format{ii}(3:3:end,:),[],'all')-1e-3*rand(1),...
                  max(CPW_my_format{ii}(2:3:end,:)./CPW_my_format{ii}(3:3:end,:),[],'all')+1e-3*rand(1)
                  ];
            recurse_quad_pts(CPW_my_format{ii},bdbx,p,0,lvls,1e-8);
end