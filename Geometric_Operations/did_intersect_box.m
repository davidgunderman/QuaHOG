function bool = did_intersect_box(X,Y)
% Checks if two bounding boxes intersect
% X The points in the first point set which will give the first bounding
% box
% Y The points in the second point set
% Output parameters
% true if an intersection was found
% false if not
    if length(X(:,1))==2
        intersection=[max(min(X(1,:)),min(Y(1,:))),min(max(X(1,:)),max(Y(1,:))) ;
        max(min(X(2,:)),min(Y(2,:))), min(max(X(2,:)),max(Y(2,:)))];
    else
        intersection=[max(min(X(1,:)),min(Y(1,:))), max(min(X(2,:)),min(Y(2,:))) , max(min(X(3,:)),min(Y(3,:)));
        min(max(X(1,:)),max(Y(1,:))), min(max(X(2,:)),max(Y(2,:))), min(max(X(3,:)),max(Y(3,:)));];
    end
    % intersection stores a new bounding box with
    % 
bool = all((intersection(:,1)<=intersection(:,2)));    
end