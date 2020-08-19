function pts = recurse_quad_pts(CPWmat,bdbx,p,initlvls,lvls,tol)
%recurse_quad Performs recursive quadtree integration
%   Inputs:
%   CPWmat      One bezier object representing one boundary of domain
%   bdbx        bdbx a quad which intersects the boundary of the domain
%   funct       Function to be integrated f = @(x,y)...
%   p           Order of tensor-product Gaussian quad rule to be used
%   lvls        Levels of adaptive quadtree subdivision to perform
%   
%   Outputs:
%   Nquads      Totol number of quadrature points used
%   value       Estimate of value of integral over domain using quadtree
pts=zeros(3,0);
if initlvls>0
    daughtbxs{1} = [bdbx(1) (bdbx(1)+bdbx(2))/2 bdbx(3) (bdbx(3) + bdbx(4))/2];
    daughtbxs{2} = [(bdbx(1)+bdbx(2))/2 bdbx(2) bdbx(3) (bdbx(3) + bdbx(4))/2];
    daughtbxs{3} = [(bdbx(1)+bdbx(2))/2 bdbx(2) (bdbx(3) + bdbx(4))/2 bdbx(4)];
    daughtbxs{4} = [bdbx(1) (bdbx(1)+bdbx(2))/2 (bdbx(3) + bdbx(4))/2 bdbx(4)];
    for i=1:4
        pts2 = recurse_quad_pts(CPWmat,daughtbxs{i},p,initlvls-1,lvls,tol);
        pts = [pts; pts2];
    end
else
    if lvls==0
        quadpts = gauss2Dbdspts(bdbx,p);
%         pts2=shapeIEvalpts(quadpts,CPWmat,1,tol);
        pts = [pts; quadpts];
%         rectangle('position',[bdbx(1) bdbx(3) bdbx(2)-bdbx(1) bdbx(4)-bdbx(3)])
%         hold on
    else
        daughtbxs{1} = [bdbx(1) (bdbx(1)+bdbx(2))/2 bdbx(3) (bdbx(3) + bdbx(4))/2];
        daughtbxs{2} = [(bdbx(1)+bdbx(2))/2 bdbx(2) bdbx(3) (bdbx(3) + bdbx(4))/2];
        daughtbxs{3} = [(bdbx(1)+bdbx(2))/2 bdbx(2) (bdbx(3) + bdbx(4))/2 bdbx(4)];
        daughtbxs{4} = [bdbx(1) (bdbx(1)+bdbx(2))/2 (bdbx(3) + bdbx(4))/2 bdbx(4)];
        for i=1:4
%             if lvls==3
%                 bdpts = [daughtbxs{i}(1) daughtbxs{i}(3);
%                         daughtbxs{i}(1) daughtbxs{i}(4);
%                         daughtbxs{i}(2) daughtbxs{i}(3);
%                         daughtbxs{i}(2) daughtbxs{i}(4);];
% %                 scatter(bdpts(:,1),bdpts(:,2),100,'r.');
%             end
            isContID=isContainedBx(daughtbxs{i},CPWmat,1,tol,p);
            isContIDpert=isContainedBx(daughtbxs{i}+[.00000001 -.00000001 .00000001 -.00000001],CPWmat,1,tol,p);
            if isContIDpert == 0
                isContID = 0;
            end
            if isContID == 0
                quadpts = gauss2Dbdspts(daughtbxs{i},p);
%                 pts2=shapIEvalpts(quadpts,CPWmat,1,tol);
%                 rectangle('Position',[daughtbxs{i}(1) daughtbxs{i}(3) daughtbxs{i}(2)-daughtbxs{i}(1) daughtbxs{i}(4)-daughtbxs{i}(3)])
%                 hold on
                pts = [pts; quadpts];
            elseif isContID == 1
                pts2 = recurse_quad_pts(CPWmat,daughtbxs{i},p,0,lvls-1,tol);
                pts = [pts; pts2];
            end
        end
    end
end
end

