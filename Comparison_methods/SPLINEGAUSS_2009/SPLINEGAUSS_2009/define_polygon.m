
function polygon_sides=define_polygon(polygon_type)

%-------------------------------------------------------------------------------
% POLYGON DEFINITION. 
%-------------------------------------------------------------------------------
%
% INPUT:
%
% [polygon_type]: PARAMETER THAT CHOOSES THE POLYGON. 
%
% OUTPUT:
% [polygon_sides]: DEFINES THE COORDINATES OF THE VERTICES (REPEATING AT THE END 
%                  THE FIRST POINT). THE VARIABLE "polygon_type" IS DEFINED ABOVE;
%                  "boundary_pts" IS A VECTOR HAVING AS COMPONENTS THE
%                  VERTICES OF THE COMPONENT, DESCRIBED COUNTERCLOCKWISE. 
%                  POLYGON IS DEFINED INSIDE THE SQUARE [0,1] x [0,1].
%                  
%-------------------------------------------------------------------------------
% OBSERVE THAT THE FIRST VERTEX IS REPEATED (TO CLOSE THE POLYGON). 
%-------------------------------------------------------------------------------

switch polygon_type
   
case 0
    k_arc=64;
    k1=k_arc; % NUMBER OF SUBDIVISIONS 1ST ARC ECLIPSE.
    k2=k_arc; % NUMBER OF SUBDIVISIONS 2ND ARC ECLIPSE.
    
    a(1)=-pi/2; a(2)=pi;
    h1=(a(2)-a(1))/k1; v1=(a(1):h1:a(2))';
    polygon_sides_1=[0.5*(cos(v1)+1) 0.5*(sin(v1)+1)];
    
    a(3)=pi/2; a(4)=0;
    h2=(a(3)-a(4))/k1; v2=(a(3)-h2:-h2:a(4))';
    polygon_sides_2=[0.5*cos(v2) 0.5*sin(v2)];
    
    polygon_sides=[polygon_sides_1; polygon_sides_2];
    
    fprintf('\n \t [POLYGON]: ECLYPSE 2, %4.0f, %4.0f',...
        size(polygon_sides_1,1),size(polygon_sides,1));
    
    plot(polygon_sides(:,1),polygon_sides(:,2),'r-');
    
case 1
    fprintf('\n \t [POLYGON]: UNIT SQUARE [0,1]^2'); 
    k=1;
    polygon_sides=k*[0 0; 1 0; 1 1; 0 1; 0 0]; 
    fprintf(' : %4.0f VERTICES',size(polygon_sides,1));
case 2
    fprintf('\n \t [POLYGON]: CONVEX POLYGON'); 
    polygon_sides=[0.1 0; 0.7 0.2; 1 0.5; 0.75 0.85; 0.5 1; 0 0.25; 0.1 0];
    fprintf(' : %4.0f VERTICES',size(polygon_sides,1));
case 3    
    fprintf('\n \t [POLYGON]: NON CONVEX POLYGON'); 
    polygon_sides=(1/4)*[1 2; 1 0; 3 2; 3 0; 4 2; 3 3; 3 0.85*4; 2 4; 0 3; 1 2];
    fprintf(' : %4.0f VERTICES',size(polygon_sides,1));
case 4
    fprintf('\n \t [POLYGON]: POLYGON'); 
    polygon_sides=(1/4)*[1 0; 3 2; 3 0; 4 2; 3 3; 3 4; 2 4; 0 3; 1 2; 1 0]; 
    fprintf(' : %4.0f VERTICES',size(polygon_sides,1));
case 5
    fprintf('\n \t [POLYGON]: POLYGON'); 
    polygon_sides=(1/4)*[0 0; 1 2; 2 0; 3 2; 4 0; 4 4; 3 2; 2 4; 1 2; 0 4; 0 0];
    fprintf(' : %4.0f VERTICES',size(polygon_sides,1));
case 6
    fprintf('\n \t [POLYGON]: UNIT SQUARE [-1,1]^2'); 
    polygon_sides=[-1 -1; 1 -1; 1 1; -1 1; -1 -1]; 
    fprintf(' : %4.0f VERTICES',size(polygon_sides,1));
case 7
    fprintf('\n \t [POLYGON]: UNIT TRIANGLE'); 
    polygon_sides=[0 0; 1 0; 1 1; 0 0];   
    fprintf(' : %4.0f VERTICES',size(polygon_sides,1));
case 8
    k_arc=64;
    k1=k_arc; % NUMBER OF SUBDIVISIONS 1ST ARC ECLIPSE.
    k2=k_arc; % NUMBER OF SUBDIVISIONS 2ND ARC ECLIPSE.
    
    a(1)=-pi/2; a(2)=pi; a(3)=3*pi/2;
    h1=(a(2)-a(1))/k1; v1=(a(1):h1:a(2))';
    h2=(a(3)-a(2))/k2; v2=(a(2):h2:a(3))';
    polygon_sides_1=[cos(v1) sin(v1)];
    
    % BUILDING OPPOSITE CONCAVITY.
    C=[cos(v2)  sin(v2)];  % ROW VECTORS!
    L=size(C,1);
    A=[C(1,:)];            % ROW VECTORS!
    B=[C(L,:)];            % ROW VECTORS!
    v=B-A; v=v/norm(v,2);
    polygon_sides_2=[];
    for index=2:L
        C_LOC=C(index,:);
        u=C_LOC-A; 
        u1=-u+2*(u*v')*v;
        polygon_sides_2=[polygon_sides_2; A+u1];
    end
    
    polygon_sides=[polygon_sides_1; polygon_sides_2];
    
    polygon_sides=0.5*ones(size(polygon_sides))+0.5*polygon_sides;
    fprintf('\n \t [POLYGON]: ECLYPSE, %4.0f, %4.0f',...
        size(polygon_sides_1,1),size(polygon_sides,1));
    
case 9  
    figure('position', get(0,'screensize'));
    axes('position', [0 0 1 1]);
    [x,y]=ginput(50);
    polygon_sides=[x y; x(1) y(1)];
    fprintf('\n \t [POLYGON]: MANUAL ENTRIES, %4.0f, %4.0f',...
        size(polygon_sides_1,1),size(polygon_sides,1));
    
case 10  
    
    k1=25; % NUMBER OF SUBDIVISIONS 1ST ARC CIRCLE.
    
    a(1)=0; a(2)=2*pi;
    h1=(a(2)-a(1))/k1; v1=(a(1):h1:a(2))';
    polygon_sides=[cos(v1) sin(v1)];
    
    polygon_sides=0.5*ones(size(polygon_sides))+0.5*polygon_sides;
    fprintf('\n \t [POLYGON]: DISK, %4.0f',size(polygon_sides,1));
case 11
    polygon_sides=[1 0; 0.75 0.5; 0 1; 0 0.5; 1 0];
    fprintf('\n \t [POLYGON]: BOOMERANG, %4.0f',size(polygon_sides,1));
case 12
    polygon_sides=[1 0; 0 1; -1 0; 0 -1; 1 0];
   fprintf('\n \t [POLYGON]: QUADRANGLE, %4.0f',size(polygon_sides,1));
case 13    
    fprintf('\n \t [POLYGON]: NON CONVEX POLYGON, HORSE LIKE'); 
    polygon_sides=(1/4)*[1 2; 1 0; 3 2; 3 0; 4 2; 3 3; 3 0.85*4; 4*0.6 4*0.95; 2 4; 0 3; 1 2];
    fprintf(' : %4.0f VERTICES',size(polygon_sides,1));
case 14
    fprintf('\n \t [POLYGON]: NON CONVEX POLYGON, HAND'); 
   polygon_sides=[   0.49804687500000   0.23832417582418
   0.52148437500000   0.29601648351648
   0.55507812500000   0.35233516483516
   0.58398437500000   0.40315934065934
   0.58242187500000   0.42513736263736
   0.56757812500000   0.43475274725275
   0.54414062500000   0.41826923076923
   0.52226562500000   0.38942307692308
   0.49570312500000   0.38392857142857
   0.48710937500000   0.41552197802198
   0.48945312500000   0.47184065934066
   0.49882812500000   0.54601648351648
   0.50195312500000   0.62019230769231
   0.49257812500000   0.64766483516484
   0.47773437500000   0.63667582417582
   0.46601562500000   0.55975274725275
   0.45039062500000   0.47733516483516
   0.45273437500000   0.57760989010989
   0.45117187500000   0.65728021978022
   0.43085937500000   0.67788461538462
   0.41757812500000   0.65178571428571
   0.41679687500000   0.57074175824176
   0.41445312500000   0.48008241758242
   0.40195312500000   0.55700549450549
   0.38632812500000   0.63667582417582
   0.37304687500000   0.64766483516484
   0.36210937500000   0.62568681318681
   0.36992187500000   0.53640109890110
   0.37695312500000   0.44848901098901
   0.34492187500000   0.50068681318681
   0.32070312500000   0.54876373626374
   0.30195312500000   0.55151098901099
   0.29882812500000   0.53090659340659
   0.32460937500000   0.45398351648352
   0.34570312500000   0.40041208791209
   0.34960937500000   0.32074175824176
   0.35976562500000   0.23832417582418
   0.49804687500000   0.23832417582418];
fprintf(' : %4.0f VERTICES',size(polygon_sides,1));
case 15
    fprintf('\n \t [POLYGON]: CUSTOM POLYGON'); 
    polygon_sides=(1/4)*[0 0; 1 1; 2 0; 3 1; 4 0; 4 4; 3 3; 2 4; 0 0];
    fprintf(' : %4.0f VERTICES',size(polygon_sides,1));
end

