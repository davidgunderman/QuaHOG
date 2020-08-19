
function [nodes_x,nodes_y,weights]=splinegauss_2009b(N,control_points,...
    rotation,P,Q,spline_order_vett,cumulative,SPLtypestring,cubature_type)

%--------------------------------------------------------------------------
% INPUTS.
%--------
%
% [N]      : DEGREE OF PRECISION OF THE ALGEBRAIC CUBATURE RULE.
%
% [control_points]: IF THE POLYGON HAS "L" SIDES, "control_points" IS A 
%           VARIABLE CONTAINING ITS VERTICES, ORDERED COUNTERCLOCKWISE. 
%           AS LAST ROW MUST HAVE THE COMPONENTS OF THE FIRST VERTEX. 
%           IN OTHER WORDS, THE FIRST ROW AND LAST ROW ARE EQUAL. 
%           "control_points" IS A "L+1 x 2" MATRIX.
%           IN CASE OF CURVILINEAR POLYGONS, "control_points" ARE
%           VERTICES OF THE CURVILINEAR POLYGON.
%
% [rotation]:SUPPOSE x_min, x_max, y_min, y_max ARE THE MINIMUM AND MAXIMUM
%           VALUES IN x AND y REACHED BY THE POLYGON, I.E.
%           "R=[x_min,x_max] x [y_min,y_max]" IS THE SMALLEST RECTANGLE 
%           WITH SIDES PARALLEL TO THE AXIS x AND y, CONTAINING THE POLYGON.
%
%           [rotation=0]: GENERAL CASE, BUT OBSERVE THAT THE FUNCTION MUST 
%                         BE DEFINED IN THE RECTANGLE "R" DESCRIBED ABOVE.
%           [rotation=1]: GOOD FOR CONVEX POLYGONS. WITH THIS CHOICE, FOR  
%                         CONVEX POLYGONS IT SUFFICES THAT THE FUNCTION IS 
%                         DEFINED IN THE POLYGON. FOR CURVILINEAR POLYGONS
%                         IT PROVIDES A GOOD REFERENCE SEGMENT BUT DOES NOT
%                         INSURE THAT ALL THE CUBATURE NODES ARE INSIDE THE
%                         DOMAIN.
%           [rotation=2]: THE USER CHOOSES A SPECIAL REFERENCE SEGMENT "PQ".
%                         CHECK [1] FOR FURTHER DETAILS. SEE THE VARIABLES 
%                         "P", "Q" BELOW.
%
% [P, Q]   : IF [rotation=2] THEN THE ALGORITHM CHOOSES A PREFERRED SEGMENT
%            "PQ" HAVING "P", "Q" AS EXTREMA, AS SET BY THE USER.
%            "P" AND "Q" ARE "1 x 2" ROW VECTORS. (SEE [1] FOR DETAILS).
%
% [spline_order_vett]: SPLINE ORDER AND BLOCK (VECTOR).
%             
%            [spline_order_vett(i,1)=2] : PIECEWISE LINEAR SPLINES.
%            [spline_order_vett(i,1)=4] : CUBIC SPLINES 
%                         (DEPENDING ON "SPLtypestring" IN INPUT).
%            [spline_order_vett(i,1)=k] : IN THE CASE k IS NOT 2 OR 4 IT 
%                          CHOOSES THE k-TH ORDER SPLINES (FOR ADDITIONAL
%                          HELP DIGIT "help spapi" IN MATLAB SHELL).                                  
%
%            [spline_order_vett(:,2)]: VECTOR OF FINAL COMPONENTS OF A BLOCK. 
%
%            EXAMPLE:
%
%            "spline_order_vett=[2 31; 4 47; 8 67]" MEANS THAT FROM THE 1st 
%             VERTEX TO THE 31th VERTEX WE HAVE AN ORDER 2 SPLINE (piecewise 
%             linear), FROM THE 32th VERTEX TO THE 47th WE USE A 4th ORDER 
%             SPLINE (i.e. A CUBIC AND PERIODIC SPLINE BY DEFAULT), FROM  
%             THE 48th TO THE 67th (AND FINAL!) WE USE AN 8th ORDER SPLINE.
%
% [cumulative]: IT CHOOSES THE PARAMETRIZATION.
%            [cumulative=0]: 1:N.
%            [cumulative=1]: CUMULATIVE.
%
% [SPLtypestring]: IF [spline_order_vett=4] IT DECIDES THE TYPE OF END 
%             CONDITIONS OF THE CUBIC SPLINE. IT IS A STRING. IT CAN BE:
%
%             'complete'   : match endslopes (as given in VALCONDS, with 
%                     default as under *default*).
%             'not-a-knot' : make spline C^3 across first and last interior  
%                     break (ignoring VALCONDS if given).
%             'periodic'   : match first and second derivatives at first 
%                     data point with those at last data point (ignoring
%                     VALCONDS if given).
%             'second'     : match end second derivatives (as given in
%                    VALCONDS, with default [0 0], i.e., as in variational).
%             'variational': set end second derivatives equal to zero
%                     (ignoring VALCONDS if given).
%
% [cubature_type]: 
%      
%            [cubature_type=0]: PADUA POINTS.
%            [cubature_type=1]: FEJER 1 TYPE (TENSORIAL).
%            [cubature_type=2]: FEJER 2 TYPE (TENSORIAL).
%            [cubature_type=3]: CLENSHAW-CURTIS (TENSORIAL).
%            [cubature_type=4]: GAUSS-LEGENDRE (TENSORIAL).
%            [cubature_type=5]: GAUSS-LEGENDRE-LOBATTO (TENSORIAL).
%
%----------
% OUTPUTS.
%----------
%
% [nodes_x, nodes_y]: THE CUBATURE RULE PRODUCES THE NODES           
%                          "(nodes_x,nodes_y)".
%            "nodes_x" AND "nodes_y" ARE COLUMN VECTORS.
%            
%
% [weights]: CUBATURE WEIGHT. COLUMN VECTOR.
%
%--------------------------------------------------------------------------
% REFERENCE PAPERS.
%-------------------
%
% [1]. A. SOMMARIVA and M. VIANELLO "Gauss-like and triangulation-free 
%      cubature over polygons". BIT Numerical Methematics 47 (2007),
%      441-453.
%
% [2]. A. SOMMARIVA and M. VIANELLO "Gauss-Green cubature and moment 
%      computation over arbitrary geometries". JCAM.
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Copyright (C) 2008 Alvise Sommariva, Marco Vianello.
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%%
%% Author:  Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>       
%% Date: June 15, 2007 - Nov. 18, 2009.
%--------------------------------------------------------------------------

%----------------------------------------------------------------------
% BOUNDARY PTS.
%----------------------------------------------------------------------
x_bd=control_points(:,1);
y_bd=control_points(:,2);

%----------------------------------------------------------------------
% "MINIMUM" RECTANGLE CONTAINING POLYGON.
%----------------------------------------------------------------------
x_min=min(x_bd); x_max=max(x_bd);
y_min=min(y_bd); y_max=max(y_bd);

%--------------------------------------------------------------------------
% POLYGON ROTATION (IF NECESSARY).
%--------------------------------------------------------------------------

switch rotation
case 0
    
    rot_matrix=eye(2);
    axis_abscissa=[x_min 0]; % NON REALLY AN AXIS. 
                             % IT IS IMPORTANT TO HAVE THE COORDINATE OF 
                             % THE MINIMUM VALUE OF THE POLYGON.
case 1
    [control_points,rot_matrix,rot_angle,axis_abscissa,P,Q]=...
        auto_rotation(control_points,[],[]);
case 2
    nrm_vect=norm(Q-P);
    if (nrm_vect > 0)
        direction_axis=(Q-P)/nrm_vect;
        [control_points,rot_matrix,rot_angle,axis_abscissa,P,Q]=...
            auto_rotation(control_points,P,Q);
    else
        [control_points,rot_matrix,rot_angle,axis_abscissa,P,Q]=...
            auto_rotation(control_points,P,Q);
    end
end

% Notation: "xi" as in the paper.
xi=axis_abscissa(1);

%--------------------------------------------------------------------------
% COMPUTE NODES.
%--------------------------------------------------------------------------

nodes_x=[];
nodes_y=[];
weights=[];

% Number of Blocks.
L=size(spline_order_vett,1);

for block_index=1:L
    
    % Initial and final indices of "control points" in the block.    
    if (block_index ==1)
        initial_point_index=1;
    else
        initial_point_index=spline_order_vett(block_index-1,2);
    end
    final_point_index=spline_order_vett(block_index,2);
    
    % Spline order in the block.
    spline_block_order=spline_order_vett(block_index,1);
    spline_block_degree=spline_block_order-1;
    
    % Control points (x_loc,y_loc) in the block.
    x_loc=control_points(initial_point_index:final_point_index,1);
    y_loc=control_points(initial_point_index:final_point_index,2);
    
    % Parametrical description of the block.
    s_loc=[];
    
    if cumulative == 0
        s_loc=(1:length(x_loc))';
    else
        s_loc(1,1)=0;
        for index=1:(length(x_loc)-1)
            s_loc(index+1,1)=s_loc(index,1)+...
                sqrt( (x_loc(index+1)-x_loc(index)).^2 + ...
                (y_loc(index+1)-y_loc(index)).^2 );
        end
    end
    
    % Computing the spline parametrical description of the block.  
    % "ppx", "ppy" describe S_i1, S_i2 in the paper, while "ppy1" 
    % describe S'_i2.
    
    switch spline_block_order 
    case 2
    case 4
        
        % CUBIC SPLINES BY CSAPE. AS DEFAULT WE USE PERIODIC CUBIC SPLINES.
        % Derivatives parameters are computed as well.
        ppx=csape(s_loc,x_loc,SPLtypestring);
        ppy=csape(s_loc,y_loc,SPLtypestring);
        [breaks_y,coeffs_y]=unmkpp(ppy);
        N_y=size(coeffs_y,1);
        dcoeffs_y=[zeros(N_y,1) 3*coeffs_y(:,1) 2*coeffs_y(:,2) ...
            coeffs_y(:,3)];
        ppy1=mkpp(breaks_y,dcoeffs_y);
        
        otherwise   
        
        ppx=spapi(spline_block_order,s_loc,x_loc);
        ppy=spapi(spline_block_order,s_loc,y_loc);
        ppy1=fnder(ppy,1);
        
    end
    
    
    
    % Every block is subdivided in "number_of_subblocks" curves determined 
    % by successive control points.
    number_of_subblocks=final_point_index-initial_point_index;
    
    % Cubature rule on the square [-1,1] x [-1,1].
    % Padua Points: 0, Gauss-Legendre: 4.
    [x_pts, y_pts, wpd]=cubature_manager(N,spline_block_degree,...
        cubature_type);
    
    
    % Computing quadrature points from a general sub-block. The cases in 
    % which the order is 2 is a little different from other spline orders. 
    % Consequently, we distinguish between them.
    for index_control_point=1:number_of_subblocks
        
        if (spline_block_order == 2)
            
            x1=x_loc(index_control_point);
            x2=x_loc(index_control_point+1);
            
            y1=y_loc(index_control_point);
            y2=y_loc(index_control_point+1);
            
            if ~(x2 == xi & x1 == xi)
                if ~( (y2-y1) == 0)
                    
                    % Computing nodes.
                    
                    qij_u= (y_pts +1)/2;
                    
                    si1_ij=x1+(x2-x1)*qij_u;
                    si2_ij=y1+(y2-y1)*qij_u;
                    
                    x=((si1_ij-xi)/2).*x_pts + (si1_ij+xi)/2;
                    y=si2_ij;
                    
                    % Applying rotation, if required.
                    rot_pts=[x y]; 
                    
                    rot_pts=rot_pts*rot_matrix;
                    
                    x_rot=rot_pts(:,1);
                    y_rot=rot_pts(:,2);
                    
                    nodes_x=[nodes_x; x_rot];
                    nodes_y=[nodes_y; y_rot];
                    
                    % Computing weights via formula (18).
                    diff_s=s_loc(index_control_point+1)-...
                        s_loc(index_control_point);
                    diff_y=y2-y1;
                    S1_i2_q_ij_u=diff_y/diff_s;
                    
                    initial_scaling_factor=(s_loc(index_control_point+1)...
                        -s_loc(index_control_point))/4; 
                    scaling_fact_minus=si1_ij-xi;
                    local_weights=initial_scaling_factor.*...
                        scaling_fact_minus.*S1_i2_q_ij_u.*wpd;
                    
                    weights=[weights; local_weights];
                    
                end
            end
            
        else
            % spline_block_order ~= 2.
            
            %-----------------------------------------
            % COMPUTING NODES.
            %-----------------------------------------
            
            % Computing "q_ij_u" (see formula (19)).
            
            half_pt_s=(s_loc(index_control_point+1)+...
                s_loc(index_control_point))/2;
            half_length_s=(s_loc(index_control_point+1)-...
                s_loc(index_control_point))/2;
            q_ij_u=half_pt_s+half_length_s*y_pts;
            
            % Evaluating S_i1(q_ij(u)), S_i2(q_ij(u)).
            switch spline_block_order 
            case 4
                S_i1_q_ij_u=ppval(ppx,q_ij_u);
                S_i2_q_ij_u=ppval(ppy,q_ij_u);
                
            otherwise
                S_i1_q_ij_u=fnval(ppx,q_ij_u);
                S_i2_q_ij_u=fnval(ppy,q_ij_u);
            end
 
            
            % See formula (18). Terms involving the spline in the arguments
            % of "f".
            scaling_fact_plus=(S_i1_q_ij_u+xi)/2; 
            scaling_fact_minus=(S_i1_q_ij_u-xi)/2;
            
            x=scaling_fact_minus.*x_pts+scaling_fact_plus;
            y=S_i2_q_ij_u;
            
            % Applying rotation, if required.
            rot_pts=[x y]; 
            rot_pts=rot_pts*rot_matrix;
         
            x_rot=rot_pts(:,1); 
            y_rot=rot_pts(:,2);
            
            nodes_x=[nodes_x; x_rot];
            nodes_y=[nodes_y; y_rot];
            
            %-----------------------------------------
            % COMPUTING WEIGHTS.
            %-----------------------------------------
            
            switch spline_block_order
            case 4
                
                S1_i2_q_ij_u=ppval(ppy1,q_ij_u);
                % SCALING FACTOR: WE PASS FROM INDEXES [i,i+1].
                initial_scaling_factor=(s_loc(index_control_point+1)-...
                    s_loc(index_control_point))/4; 
                
                otherwise
                
                S1_i2_q_ij_u=fnval( fnder( ppy, 1), q_ij_u );
                % SCALING FACTOR: WE PASS FROM INDEXES [i,i+1].
                initial_scaling_factor=(s_loc(index_control_point+1)-...
                    s_loc(index_control_point))/4; 
                
            end
            
            local_weights=initial_scaling_factor.*...
                (2*scaling_fact_minus).*S1_i2_q_ij_u.*wpd;
            
            weights=[weights; local_weights];
        end
    end
    
end





%--------------------------------------------------------------------------
% FUNCTIONS USED IN THE ALGORITHM.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% auto_rotation.
%--------------------------------------------------------------------------

function [control_pts_rot,rot_matrix,rot_angle,axis_abscissa,...
    vertex_1,vertex_2]=auto_rotation(control_pts,vertex_1,vertex_2)

%--------------------------------------------------------------------------
% OBJECT:
%-----------
%
% AUTOMATIC ROTATION OF A CONVEX POLYGON SO THAT "GAUSSIAN POINTS" 
% ARE ALL (OR "MOST OF THEM" IF THE DEGREE OF THE BORDER IS BIGGER THAN 2) 
% CONTAINED IN THE CONVEX POLYGON. SEE THE PAPER FOR DETAILS.
%
%--------------------------------------------------------------------------
% INPUTS:
%----------
%
% control_pts : CONTROL POINTS OF THE CURVILEAR POLYGON.
%
% vertex_1, vertex_2: ROTATION VERTICES (IF ANY ARE CHOOSEN).
%
%----------
% OUTPUTS:
%----------
%
% control_pts_rot: CONTROL POINTS OF THE CURVILEAR POLYGON, 
%       AFTER THE ROTATION.
%
% rot_matrix: ROTATION MATRIX.
%
% rot_angle: ROTATION ANGLE.
%
% axis_abscissa: ABSCISSA OF THE REFERENCE AXIS.
%
% vertex_1, vertex_2: ROTATION VERTICES (IF ANY ARE CHOOSEN OR COMPUTED BY 
%       THE ROUTINE).
%
%--------------------------------------------------------------------------
% ADDITIONAL ROUTINES :
%----------------------
%
% 1. points2distances
%
%--------------------------------------------------------------------------


% FIND DIRECTION AND ROTATION ANGLE.
if length(vertex_1) == 0
    % COMPUTING ALL THE DISTANCES BETWEEN POINTS.A LITTLE TIME CONSUMING AS
    % PROCEDURE.
    distances = points2distances(control_pts);  
    [max_distances,max_col_comp]=max(distances,[],2);
    [max_distance,max_row_comp]=max(max_distances,[],1);
    vertex_1=control_pts(max_col_comp(max_row_comp),:);
    vertex_2=control_pts(max_row_comp,:);
    direction_axis=(vertex_2-vertex_1)/max_distance;
else
    direction_axis=(vertex_2-vertex_1)/norm(vertex_2-vertex_1);
end

rot_angle_x=acos(direction_axis(1));
rot_angle_y=acos(direction_axis(2));

if rot_angle_y <= pi/2
    if rot_angle_x <= pi/2
        rot_angle=-rot_angle_y;
    else
        rot_angle=rot_angle_y;
    end
else
    if rot_angle_x <= pi/2
        rot_angle=pi-rot_angle_y;
    else
        rot_angle=rot_angle_y;
    end
end


% CLOCKWISE ROTATION.
rot_matrix=[cos(rot_angle) sin(rot_angle); 
    -sin(rot_angle) cos(rot_angle)];

number_sides=size(control_pts,1)-1;

control_pts_rot=(rot_matrix*control_pts')';

axis_abscissa=rot_matrix*vertex_1';





%--------------------------------------------------------------------------
% points2distances.
%--------------------------------------------------------------------------

function distances = points2distances(points)  

% Create euclidean distance matrix from point matrix.

% Get dimensions.
[numpoints,dim]=size(points);                    

% All inner products between points.
distances=points*points';                     

% Vector of squares of norms of points.
lsq=diag(distances);                            

% Distance matrix.
distances=sqrt(repmat(lsq,1,numpoints)+repmat(lsq,1,numpoints)'-2*distances);






%--------------------------------------------------------------------------
% cubature_manager.
%--------------------------------------------------------------------------

function [nodes_x, nodes_y, weights]=cubature_manager(M,p_i,cubature_type)


%--------------------------------------------------------------------------
% OBJECT:
%-----------
%
% THE TARGET OF THIS ROUTINE IS TO PROVIDE NODES AND WEIGHTS OF A
% QUADRATURE ROUTINE ON [-1,1]. THE CODES ARE DESCRIBED BY L.N. TREFETHEN 
% IN HIS CLENSHAW-CURTIS PAPER AND BY WALDVOGEL (PUBLISHED BY BIT).
%
%--------------------------------------------------------------------------
% INPUTS:
%----------
%
% M : DEGREE OF PRECISION OF 2D CUBATURE RULE.
%
% p_i: SPLINE DEGREE
%
% cubature_type: IT POINTS A CUBATURE ON THE SQUARE [-1,1] x [-1,1].
%
%     [cubature_type=0]: PADUA POINTS.
%     [cubature_type=1]: FEJER 1 (TENSORIAL).
%     [cubature_type=2]: FEJER 2 (TENSORIAL).
%     [cubature_type=3]: CLENSHAW CURTIS (VIA WEIGHTS) (TENSORIAL).
%     [cubature_type=4]: GAUSS-LEGENDRE (TENSORIAL).
%     [cubature_type=5]: GAUSS-LEGENDRE-LOBATTO (TENSORIAL).
%
%----------
% OUTPUTS:
%----------
%
% nodes : M x 2 MATRIX OF NODES, IN THE INTERVAL [-1,1].
% 
% weights: M x 1 COLUMN VECTOR OF WEIGHTS.
%
%--------------------------------------------------------------------------
% ADDITIONAL ROUTINES :
%----------------------
%
% 1. quadrature_rules_1D
% 2. pdcub
%
%--------------------------------------------------------------------------

switch cubature_type
    
case 4
    
    % lower degree of a Gauss-Legendre rule with ADE equal to M.
    Ntau=floor(M/2)+1; 
    % In the quadrature rule, for M as input we have a rule with M+1 pts.
    Ntau_cub=Ntau-1; 
    [nodes_tau,weights_tau]=quadrature_rules_1D(Ntau_cub,cubature_type);
    
    Nu_deg=(M+2)*p_i-1;
    % lower degree of a Gauss-Legendre rule with ADE equal to Nu_deg.
    Nu=floor(Nu_deg/2)+1; 
    Nu_cub=Nu-1;
    [nodes_u,weights_u]=quadrature_rules_1D(Nu_cub,cubature_type);
    
    [nodes_x,nodes_y]=meshgrid(nodes_tau,nodes_u);
    nodes_x=nodes_x(:);
    nodes_y=nodes_y(:);
    
    [weights1,weights2]=meshgrid(weights_tau,weights_u);
    weights=weights1(:).*weights2(:);
    
case 0
    
    % Degree of Padua points for achieving degree of precision "N".
    Ntau=(M+2)*p_i-1;
    Nu=Ntau;
    % Computing Padua points and weights in the square [-1,1]x[-1,1]
    PadL = pdcub(Ntau,[-1 1 -1 1]);
    % Working with vectors.
    nodes_x=PadL(:,1); nodes_y=PadL(:,2); weights=PadL(:,3); 
    
    otherwise
    
    Ntau=M; % lower degree of a rule with ADE equal to M.
    [nodes_tau,weights_tau]=quadrature_rules_1D(Ntau,cubature_type);
    % lower degree of a rule with ADE equal to Nu_deg=(M+2)*p_i-1.
    Nu=(M+2)*p_i-1; 
    
    [nodes_u,weights_u]=quadrature_rules_1D(Nu,cubature_type);
    
    [nodes_x,nodes_y]=meshgrid(nodes_tau,nodes_u);
    nodes_x=nodes_x(:);
    nodes_y=nodes_y(:);
    
    [weights1,weights2]=meshgrid(weights_tau,weights_u);
    weights=weights1(:).*weights2(:);
    
end





%--------------------------------------------------------------------------
% pdcub.
%--------------------------------------------------------------------------

function [varargout] = pdcub(n,xyrange,varargin)

%--------------------------------------------------------------------------
% OBJECT.
%---------
%
% PadL = pdcub(n,xyrange)
% [cubature,PadL] = pdcub(n,xyrange,funct,opt1,opt2,...)
%
% This function computes the Padua points defined in the rectangle 
% [xyrange(1),xyrange(2)] x [xyrange(3),xyrange(4)] and the 
% corresponding cubature weights as a matrix of abscissas (first
% column), ordinates (second column) and weights (third column).
% Optionally, it computes the cubature of the function f (optional 
% input argument).
%
%--------------------------------------------------------------------------
% INPUTS.
%---------
%
% n        : interpolation degree.
% xyrange  : a vector [a,b,c,d] defining the rectangle [a,b] x [c,d].
% funct    : function to be interpolated in the form 
%            funct(x,y,opt1,opt2,...), where opt1, opt2, ... are
%            optional arguments for funct. 
%
% OUTPUT.
%
% PadL     : a matrix with the abscissas of Padua points in the 
%            first column, the ordinates in the second and the 
%            cubature weights in the third.
% cubature : cubature of the integrand funct.
%
%--------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS CODE.
%--------------------------------
%
% 1. pdpts
% 2. pdwtsMM
%
%--------------------------------------------------------------------------
% EXAMPLES.
%-----------
%
% 1)
% Compute the Padua points and the cubature weights of degree 50
%
% xyrange = [0,1,0,1];
% PadL = pdcub(50,xyrange);
%
% 2)
% Compute the cubature of the Franke's function defined in funct.m
%
% xyrange = [0,1,0,1];
% [cubature,PadL] = pdcub(50,xyrange,@f);
%
%--------------------------------------------------------------------------
% PROCEDURES USED BY THE CODE:
%------------------------------
%
% 1. pdpts
% 2. pdwtsMM
%
% THESE FUNCTIONS ARE ATTACHED BELOW.
% 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Copyright (C) 2008-2009 
% Marco Caliari, Stefano De Marchi, Alvise Sommariva, Marco Vianello.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%
% Author:  
%          Marco Caliari     <marco.caliari@univr.it>
%          Stefano De Marchi <demarchi@euler.math.unipd.it>   
%          Alvise Sommariva  <alvise@euler.math.unipd.it>
%          Marco Vianello    <marcov@euler.math.unipd.it>   
%
% Date: November 19, 2009.
%--------------------------------------------------------------------------

% Compute the Padua points in the rectangle defined by xyrange
PadL = pdpts(n,xyrange);
% Compute the cubature weights
PadL(:,3) = pdwtsMM(n,xyrange);
if (nargin < 2)
    error('Too few input arguments')
elseif (nargin >= 3)
    funct = varargin{1};
    varargout{1} = PadL(:,3)'*feval(funct,PadL(:,1),PadL(:,2),varargin{2:end});
    varargout{2} = PadL;
else
    varargout{1} = PadL;
end





%--------------------------------------------------------------------------
% pdwtsMM.
%--------------------------------------------------------------------------

function [varargout] = pdwtsMM(n,varargin)

%--------------------------------------------------------------------------
% USAGE of "pdwtsMM".
%
% L = pdwtsMM(n)
% L = pdwtsMM(n,xyrange)
% [L1,L2,L] = pdwtsMM(n)
% [L1,L2,L] = pdwtsMM(n,xyrange)
%
% Compute the cubature weights L so that, if Pad is the
% matrix of Padua points computed through the call
% Pad = pdpts(n) or Pad = pdpts(n,xyrange), then the cubature of 
% the function f is given by L'*f(Pad(:,1),Pad(:,2)). 
% Otherwise, one can compute the cubature weights L1 and L2, so 
% that, if X1,Y1 and X2,Y2 are the subgrids of Padua points 
% computed through the call [X1,Y1,X2,Y2] = pdpts(n) or 
% [X1,Y1,X2,Y2] = pdpts(n,xyrange), then the cubature of the 
% function f is given by 
% sum(sum(L1.*f(X1,Y1)))+sum(sum(L1.*f(X2,Y2))).
%--------------------------------------------------------------------------
% INPUT.    
%
% n       : interpolation degree
% xyrange : an optional vector [a,b,c,d] defining the rectangle 
%           [a,b] x [c,d]. By default, xyrange = [-1,1,-1,1]
%
% OUTPUT.   
%
% L       : cubature weights associated to the matrix Pad of 
%           Padua points
% L1,L2   : cubature weights associated to the subgrids 
%           X1,Y1 and X2,Y2 of Padua poins.
%--------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS CODE:
% no external function is used.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Copyright (C) 2008-2009 
% Marco Caliari, Stefano De Marchi, Alvise Sommariva, Marco Vianello.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%
% Author:  
%          Marco Caliari     <marco.caliari@univr.it>
%          Stefano De Marchi <demarchi@euler.math.unipd.it>   
%          Alvise Sommariva  <alvise@euler.math.unipd.it>
%          Marco Vianello    <marcov@euler.math.unipd.it>   
%
% Date: November 19, 2009.
%--------------------------------------------------------------------------

if (nargin == 1)
    xyrange = [-1,1,-1,1];
else
    xyrange = varargin{1};
end
if (n == 0)
    % degree 0
    L1 = (xyrange(2)-xyrange(1))*(xyrange(4)-xyrange(3));
    L2 = zeros(1,0);
    L = L1;
else
    argn = linspace(0,pi,n+1);
    argn1 = linspace(0,pi,n+2);
    k = [0:2:n]';
    l = (n-mod(n,2))/2+1;
    % even-degree Chebyshev polynomials on the subgrids
    TE1 = cos(k*argn(1:2:n+1));
    TE1(2:l,:) = TE1(2:l,:)*sqrt(2);
    TO1 = cos(k*argn(2:2:n+1));
    TO1(2:l,:) = TO1(2:l,:)*sqrt(2);
    TE2 = cos(k*argn1(1:2:n+2));
    TE2(2:l,:) = TE2(2:l,:)*sqrt(2);
    TO2 = cos(k*argn1(2:2:n+2));
    TO2(2:l,:) = TO2(2:l,:)*sqrt(2);
    % even,even moments matrix
    mom = 2*sqrt(2)./(1-k.^2);
    mom(1) = 2;
    [M1,M2] = meshgrid(mom);
    M = M1.*M2;
    Mmom = fliplr(triu(fliplr(M)));
    % interpolation weights matrices
    W1 = 2*ones(l)/(n*(n+1));
    W2 = 2*ones((n+mod(n,2))/2+1,(n+mod(n,2))/2)/(n*(n+1));
    W1(:,1) = W1(:,1)/2;
    W2(1,:) = W2(1,:)/2;
    if (mod(n,2) == 0)
        %  Mmom(n/2+1,1) = Mmom(n/2+1,1)/2;
        Mmom(1,n/2+1) = Mmom(1,n/2+1)/2;
        W1(:,n/2+1) = W1(:,n/2+1)/2;
        W1(n/2+1,:) = W1(n/2+1,:)/2;
    else
        W2((n+1)/2+1,:) = W2((n+1)/2+1,:)/2;
        W2(:,(n+1)/2) = W2(:,(n+1)/2)/2;
    end
    % cubature weights as matrices on the subgrids.
    L1 = W1.*(TO2'*Mmom*TE1);
    L2 = W2.*(TE2'*Mmom*TO1);
    if (mod(n,2) == 0)
        L = zeros(n/2+1,n+1);
        L(:,1:2:n+1) = L1;
        L(:,2:2:n+1) = L2;
        L = L(:);
    else
        L = zeros((n+1)/2,(n+2));
        L = [L1',L2']';
        L = L(:);
    end
    L = L*(xyrange(2)-xyrange(1))*(xyrange(4)-xyrange(3))/4;
end
if (nargout == 0 | nargout == 1)
    varargout{1} = L;
else
    varargout{1} = L1;
    varargout{2} = L2;
    varargout{3} = L;
end





%--------------------------------------------------------------------------
% pdpts.
%--------------------------------------------------------------------------

function [varargout] = pdpts(n,varargin)

%--------------------------------------------------------------------------
% OBJECT.
%---------
%
% Pad = pdpts(n)
% Pad = pdpts(n,xyrange)
% [X1,Y1,X2,Y2] = pdpts(n)
% [X1,Y1,X2,Y2] = pdpts(n,xyrange)
%
% Compute the (first family of) Padua points, either as a matrix Pad 
% with their abscissas in the first column and their ordinates in the
% second, or as two (mesh)grids X1,Y1 and X2,Y2, respectively.
%--------------------------------------------------------------------------
% INPUTS.   
%--------
%
% n           : interpolation degree
% xyrange     : an optional vector [a,b,c,d] defining the rectangle 
%               [a,b] x [c,d]. Otherwise, xyrange = [-1,1,-1,1]
%
%----------
% OUTPUTS.  
%----------
%
% Pad         : matrix of size ((n+1)*(n+2)/2) x 2 such
%               that (Pad(:,1),Pad(:,2)) defines the Padua points in the
%               rectangle [xyrange(1),xyrange(2)] x [xyrange(3),xyrange(4)].  
% X1,Y1,X2,Y2 : the two (mesh)grids X1,Y1 and X2,Y2 defining the Padua 
%               points
%--------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS CODE.
%-------------------------------
%
% No external function is used.
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Copyright (C) 2008-2009 
% Marco Caliari, Stefano De Marchi, Alvise Sommariva, Marco Vianello.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%
% Author:  
%          Marco Caliari     <marco.caliari@univr.it>
%          Stefano De Marchi <demarchi@euler.math.unipd.it>   
%          Alvise Sommariva  <alvise@euler.math.unipd.it>
%          Marco Vianello    <marcov@euler.math.unipd.it>   
%
% Date: November 19, 2009.
%--------------------------------------------------------------------------

if (nargin == 1)
    % standard square [-1,1] x [-1,1]
    xyrange = [-1,1,-1,1];
else
    % rectangle [xyrange(1),xyrange(2)] x [xyrange(3),xyrange(4)]
    xyrange = varargin{1};
end
if (n == 0)
    % degree 0
    if (nargout ~= 4)
        % points as a single matrix
        Pad = [xyrange(1),xyrange(3)];
        varargout = {Pad};
    else
        % points as two (mesh)grids
        X1 = xyrange(1);
        Y1 = xyrange(3);
        X2 = zeros(1,0);
        Y2 = zeros(1,0);
        varargout = {X1,Y1,X2,Y2};
    end  
else
    % degree > 0
    zn = (xyrange(1)+xyrange(2)+(xyrange(2)-xyrange(1))*...
        cos(linspace(0,1,n+1)*pi))/2;
    zn1 = (xyrange(3)+xyrange(4)+(xyrange(4)-xyrange(3))*...
        cos(linspace(0,1,n+2)*pi))/2;
    if (nargout ~= 4)
        % points as a single matrix
        [Pad1,Pad2] = meshgrid(zn,zn1);
        [M1,M2] = meshgrid([0:n],[0:n+1]);
        findM = find(mod(M1+M2,2));
        Pad = [Pad1(findM),Pad2(findM)];
        varargout = {Pad};
    else
        % points as two (mesh)grids  
        En = zn(1:2:n+1);
        On = zn(2:2:n+1);
        En1 = zn1(1:2:n+2);
        On1 = zn1(2:2:n+2);
        [X1,Y1] = meshgrid(En,On1);
        [X2,Y2] = meshgrid(On,En1);
        varargout = {X1,Y1,X2,Y2};
    end  
end 





%--------------------------------------------------------------------------
% quadrature_rules_1D.
%--------------------------------------------------------------------------

function [nodes,weights]=quadrature_rules_1D(n,quadrature_type)

%--------------------------------------------------------------------------
% OBJECT:
%-----------
% THE TARGET OF THIS ROUTINE IS TO PROVIDE NODES AND WEIGHTS OF A
% QUADRATURE ROUTINE ON [-1,1]. THE CODES ARE DESCRIBED BY L.N. TREFETHEN 
% IN HIS CLENSHAW-CURTIS PAPER AND BY WALDVOGEL (PUBLISHED BY BIT).
%
%--------------------------------------------------------------------------
% INPUTS:
%----------
%
% n : NUMBER OF NODES OF THE QUADRATURE RULE (NOT THE DEGREE!!).
%
% quadrature_type: IT POINTS A QUADRATURE RULE
%           [quadrature_type=1]: FEJER 1.
%           [quadrature_type=2]: FEJER 2.
%           [quadrature_type=3]: CLENSHAW CURTIS (VIA WEIGHTS).
%           [quadrature_type=4]: GAUSS-LEGENDRE.
%           [quadrature_type=5]: GAUSS-LEGENDRE-LOBATTO.
%           [quadrature_type=6]: COMPOSITE TRAPEZOIDAL RULE.
%
%----------
% OUTPUTS:
%----------
%
% nodes : M x 2 MATRIX OF NODES, IN THE INTERVAL [-1,1].
% 
% weights: M x 1 COLUMN VECTOR OF WEIGHTS.
%
%--------------------------------------------------------------------------
% ADDITIONAL ROUTINES :
%----------------------
%
% 1. r_jacobi
% 2. gauss
% 3. lobatto_jacobi
%
% THESE ROUTINES ARE WRITTEN BY D. LAURIE AND W. GAUTSCHI, AND CAN BE FOUND
% IN W. GAUTSCHI HOMEPAGE. THEY ARE ATTACHED IN THIS FILE (SEE THE BOTTOM
% OF THE FILE). NO EXTERNAL FILE IS NEEDED.
%
%--------------------------------------------------------------------------
% EXAMPLE:
%----------
%
% >> [nodes,weights]=quadrature_rules_1D(5,5)
% 
% nodes =
%
%   -1.0000
%   -0.8302
%   -0.4688
%    0.0000
%    0.4688
%    0.8302
%    1.0000
%
%
% weights =
%
%    0.0476
%    0.2768
%    0.4317
%    0.4876
%    0.4317
%    0.2768
%    0.0476
%
%--------------------------------------------------------------------------
% RELATED PAPERS:
%-----------------
%
% [1] W. GAUTSCHI, ORTHOGONAL POLYNOMIALS AND QUADRATURE.
%
% [2] LLOYD N. TREFETHEN, IS GAUSS QUADRATURE BETTER THAN CLENSHAW?CURTIS?
%
% [3] J. WALDVOGEL, FAST CONSTRUCTION OF THE FEJER AND CLENSHAW-CURTIS
%     QUADRATURE RULES.
%
%--------------------------------------------------------------------------

switch quadrature_type

    case 1 % FEJER 1.
        N=[1:2:n-1]'; l=length(N); m=n-l; K=[0:m-1]';
        v0=[2*exp(i*pi*K/n)./(1-4*K.^2); zeros(l+1,1)];
        v1=v0(1:end-1)+conj(v0(end:-1:2));
        weights=ifft(v1);
        k=(1/2):(n-(1/2)); nodes=(cos(k*pi/n))';

    case 2 % FEJER 2.
        N=[1:2:n-1]'; l=length(N); m=n-l; K=[0:m-1]';
        v0=[2./N./(N-2); 1/N(end); zeros(m,1)];
        v2=-v0(1:end-1)-v0(end:-1:2);
        wf2=ifft(v2); weights=[wf2;0];
        k=0:n; nodes=(cos(k*pi/n))';

    case 3 % CLENSHAW CURTIS.
        n=n-1;
        N=[1:2:n-1]'; l=length(N); m=n-l; K=[0:m-1]';
        g0=-ones(n,1); g0(1+l)=g0(1+l)+n; g0(1+m)=g0(1+m)+n;
        g=g0/(n^2-1+mod(n,2));
        v0=[2./N./(N-2); 1/N(end); zeros(m,1)];
        v2=-v0(1:end-1)-v0(end:-1:2);
        wcc=ifft(v2+g); weights=[wcc;wcc(1,1)];
        k=0:n; nodes=(cos(k*pi/n))';

    case 4 % GAUSS LEGENDRE.
        beta=0.5./sqrt(1-(2*(1:n)).^(-2));
        T=diag(beta,1)+diag(beta,-1);
        [V,D]=eig(T);
        x=diag(D); [x,index]=sort(x); x=x';
        w=2*V(1,index).^2;
        nodes=x';
        weights=w';

    case 5 % GAUSS LEGENDRE LOBATTO.
        xw=lobatto_jacobi(n);
        nodes=xw(:,1);
        weights=xw(:,2);

    case 6 % TRAPEZOIDAL RULE.
        h=2/n;
        x=-1:h:1;
        w=h*[0.5 ones(1,length(x)-2) 0.5];
        nodes=x';
        weights=w';

end





%--------------------------------------------------------------------------
% ADDITIONAL ROUTINES BY W. GAUTSCHI.
%--------------------------------------------------------------------------

%---------------------------
% r_jacobi.
%---------------------------

function ab=r_jacobi(N,a,b)

nu=(b-a)/(a+b+2);
mu=2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);
if N==1
 ab=[nu mu]; return
end

N=N-1;
n=1:N;
nab=2*n+a+b;
nuadd=(b^2-a^2)*ones(1,N)./(nab.*(nab+2));
A=[nu nuadd];
n=2:N;
nab=nab(n);
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
abadd=[mu; B1; B'];
ab=[A' abadd];





%---------------------------
% gauss.
%---------------------------

function xw=gauss(N,ab)
N0=size(ab,1); if N0<N, error('input array ab too short'), end
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
  J(n,n-1)=sqrt(ab(n,2));
  J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
[D,I]=sort(diag(D));
V=V(:,I);
xw=[D ab(1,2)*V(1,:)'.^2];





%---------------------------
% lobatto_jacobi.
%---------------------------

function xw=lobatto_jacobi(N,a,b)
if nargin<2, a=0; end; if nargin<3, b=a; end
ab=r_jacobi(N+2,a,b);
ab(N+2,1)=(a-b)/(2*N+a+b+2);
ab(N+2,2)=4*(N+a+1)*(N+b+1)*(N+a+b+1)/((2*N+a+b+1)* ...
(2*N+a+b+2)^2);
xw=gauss(N+2,ab);


