function [Nquads,value] = recurse_quad(CPWmat,bdbx,funct,p,initlvls,lvls,tol)
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
Nquads = 0;
value = 0;
if initlvls>0
    daughtbxs{1} = [bdbx(1) (bdbx(1)+bdbx(2))/2 bdbx(3) (bdbx(3) + bdbx(4))/2];
    daughtbxs{2} = [(bdbx(1)+bdbx(2))/2 bdbx(2) bdbx(3) (bdbx(3) + bdbx(4))/2];
    daughtbxs{3} = [(bdbx(1)+bdbx(2))/2 bdbx(2) (bdbx(3) + bdbx(4))/2 bdbx(4)];
    daughtbxs{4} = [bdbx(1) (bdbx(1)+bdbx(2))/2 (bdbx(3) + bdbx(4))/2 bdbx(4)];
    for i=1:4
        [nqt,vt] = recurse_quad(CPWmat,daughtbxs{i},funct,p,initlvls-1,lvls,tol);
        value = value + vt; Nquads = Nquads + nqt;
    end
else
    if lvls==0
        quadpts = gauss2Dbdspts(bdbx,p);
        quadpts(:,1:2)=shapIEvalpts(quadpts(:,1:2),CPWmat,
        shapeIEvalpts([x y],funct,CPWmat,1,tol);
        quadpts = gauss2Dbds(Ifun, bdbx, p);
        Nquads = p^2;
    else
        daughtbxs{1} = [bdbx(1) (bdbx(1)+bdbx(2))/2 bdbx(3) (bdbx(3) + bdbx(4))/2];
        daughtbxs{2} = [(bdbx(1)+bdbx(2))/2 bdbx(2) bdbx(3) (bdbx(3) + bdbx(4))/2];
        daughtbxs{3} = [(bdbx(1)+bdbx(2))/2 bdbx(2) (bdbx(3) + bdbx(4))/2 bdbx(4)];
        daughtbxs{4} = [bdbx(1) (bdbx(1)+bdbx(2))/2 (bdbx(3) + bdbx(4))/2 bdbx(4)];
        for i=1:4
            if lvls==3
                bdpts = [daughtbxs{i}(1) daughtbxs{i}(3);
                        daughtbxs{i}(1) daughtbxs{i}(4);
                        daughtbxs{i}(2) daughtbxs{i}(3);
                        daughtbxs{i}(2) daughtbxs{i}(4);];
%                 scatter(bdpts(:,1),bdpts(:,2),100,'r.');
            end
            isContID=isContainedBx(daughtbxs{i},CPWmat,1,tol,p);
            if isContID == 0
                value = value + gauss2Dbds(funct,daughtbxs{i},p);
                Nquads = Nquads + p^2;
            elseif isContID == 1
                [nqt,vt] = recurse_quad(CPWmat,daughtbxs{i},funct,p,0,lvls-1,tol);
                value = value + vt; Nquads = Nquads + nqt;
            end
        end
    end
end
end

