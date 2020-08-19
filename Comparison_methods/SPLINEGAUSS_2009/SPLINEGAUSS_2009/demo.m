
clearvars;

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
% THIS DEMO SHOWS HOW TO USE THE PROCEDURE "splinegauss" FOR COMPUTING CU-
% BATURE OF CONTINUOUS FUNCTIONS OVER CERTAIN CURVILINEAR POLYGONS .
%--------------------------------------------------------------------------
% FUNCTIONS USED.
%
% 1. define_polygon
% 2. splinegauss_2009b
% 3. fct2D
%--------------------------------------------------------------------------
%% Copyright (C) 2007-2009. Alvise Sommariva and Marco Vianello 
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
%--------------------------------------------------------------------------
%% Authors:  
%% Marco Vianello    <marcov@euler.math.unipd.it>
%% Alvise Sommariva  <alvise@euler.math.unipd.it>
%% Date: November 19, 2009.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% SEVERAL POLYGONS ARE DEFINED IN THE m-FILE "define_polygon.m". 
%
% [polygon_type=1]: UNIT SQUARE [0,1]^2.
% [polygon_type=2]: A 6 SIDES CONVEX POLYGON USED IN PAPER.
% [polygon_type=3]: A 9 SIDES NONCONVEX POLYGON USED IN THE PAPER.
% [polygon_type=4]: POLYGON.
% [polygon_type=5]: POLYGON.
% [polygon_type=6]: UNIT SQUARE [-1,1]^2.
% [polygon_type=7]: UNIT TRIANGLE.
% [polygon_type=8]: LUNE LIKE DOMAIN.
% [polygon_type=9]: POINTS ON MONITOR.
% [polygon_type=10]: CIRCLE.
% [polygon_type=11]: PENTANGLE.
% [polygon_type=12]: QUADRANGLE.
% [polygon_type=13]: HORSE LIKE DOMAIN.
% [polygon_type=14]: HAND.
% [polygon_type=15]: CUSTOM POLYGON 1.
%
% TO DEFINE YOUR OWN POLYGON, PLEASE OPEN THE m-FILE "define_polygon.m" 
% AND ADD IT TO THE LIST. OBSERVE THAT THE FIRST VERTEX IS REPEATED AS 
% FINAL VERTEX. 
% FOR EXAMPLE, THE UNIT SQUARE IS DEFINED AS "[0 0; 1 0; 1 1; 0 1; 0 0]"
% AND NOT "[0 0; 1 0; 1 1; 0 1]".
%--------------------------------------------------------------------------
polygon_type=10;

%--------------------------------------------------------------------------
% SEVERAL FUNCTIONS ARE DEFINED IN THE m-FILE "fct2D.m". MORE PRECISELY
%
% [function_type=1]: z=franke(x,y);
% [function_type=2]: z=( (x-0.5).^2 +(y-0.5).^2 ).^(1/2);
% [function_type=3]: z=(x+y).^k;
% [function_type=4]: z=exp(loc_arg); loc_arg=-( (x-0.5).^2+(y-0.5).^2 ); 
% [function_type=5]: z=exp(loc_arg); loc_arg=-100*((x-0.5).^2+(y-0.5).^2); 
% [function_type=6]: z=cos(30*(x+y));
% [function_type=7]: z=ones(size(x));
% [function_type=8]: z=exp(x+y);
% [function_type=9]: z=1./den; den=1+16*(x.^2+y.^2); 
% [function_type=10]:z=(x^2+y^2)^(3/2)
%
% TO DEFINE YOUR OWN FUNCTION, PLEASE OPEN THE m-FILE "fct2D.m" AND ADD IT 
% TO THE LIST. 
%--------------------------------------------------------------------------
function_type=3;

%--------------------------------------------------------------------------
% [spline_order_vett]: SPLINE ORDER AND BLOCK (VECTOR).
%
%            [spline_order_vett(:,1)]: VECTOR OF SPLINES ORDERS OF EACH
%                                      CURVILINEAR SIDE (aka "BLOCK").
%            PARTICULAR CASES:
%            [spline_order_vett(i,1)=2] : PIECEWISE LINEAR SPLINES.
%            [spline_order_vett(i,1)=4] : PERIODIC CUBIC SPLINES (CAN  
%                                         CHANGE DEFAULT, MODIFYING THE 
%                                         CODE).
%            [spline_order_vett(i,1)=k] : IN THE CASE k IS NOT 2 OR 4 IT 
%                                         CHOOSES THE k-TH ORDER SPLINES 
%                                        (FOR ADDITIONAL HELP, DIGIT  
%                                         "help spapi"IN MATLAB SHELL).
%
%            [spline_order_vett(:,2)]   : VECTOR OF FINAL COMPONENTS OF 
%                                         A BLOCK. 
%
%             EXAMPLE:
%
%            "spline_order_vett=[2 31; 4 47; 8 68]" MEANS THAT FROM THE 1st 
%             VERTEX TO THE 31th VERTEX WE HAVE AN ORDER 2 SPLINE 
%             (piecewise linear), FROM THE 32th VERTEX TO THE 47th WE USE 
%             A 4th ORDER SPLINE (i.e. A CUBIC AND PERIODIC SPLINE BY 
%             DEFAULT), FROM THE 48th TO THE 67th (AND FINAL!) WE USE AN
%             8th ORDER SPLINE.
%--------------------------------------------------------------------------
% spline_order_vett=[repmat(3 ,64,1) (3:2:129)'];
% spline_order_vett=[repmat(3 ,4,1) (3:2:9)'];
spline_order_vett=[4 26];
% [4 65; 4 129]

%--------------------------------------------------------------------------
% [spline_type]: IF [spline_order_vett=4] IT DECIDES THE TYPE OF END 
%                CONDITIONS OF THE CUBIC SPLINE. IT IS A STRING. IT CAN BE:
%
%             1. 'complete'    : match endslopes (as given in VALCONDS, 
%                                with default as under *default*).
%             2. 'not-a-knot'  : make spline C^3 across first and last 
%                                interior break (ignoring VALCONDS if 
%                                given).
%             3. 'periodic'    : match first and second derivatives at 
%                                first data point with those at last data 
%                                point (ignoring VALCONDS if given).
%             4. 'second'      : match end second derivatives (as given in 
%                                VALCONDS,  with default [0 0], i.e., as in
%                                variational).
%             5. 'variational' : set end second derivatives equal to zero
%                               (ignoring VALCONDS if given).
%
%--------------------------------------------------------------------------
spline_type=3;

%--------------------------------------------------------------------------
% [cumulative]: THIS VARIABLE CHOOSES THE PARAMETRIZATION. 
%              
%              EXAMPLES:
%              "cumulative=0": 1:N.
%              "cumulative=1": CUMULATIVE.
%--------------------------------------------------------------------------
cumulative=1;

%--------------------------------------------------------------------------
% WRITE THE NUMBER OF POINTS OF THE 1-DIMENSIONAL QUADRATURE RULE USED BY 
% "splinegauss". IT CORRESPONDS TO COMPUTE A QUADRATURE RULE ON THE DOMAIN
% HAVING ALGEBRAIC DEGREE OF PRECISION "2*N-1".
%--------------------------------------------------------------------------
N=25; 

%--------------------------------------------------------------------------
% FOR THE DESCRIPTION OF HOW A ROTATION WORKS SEE PAPER [1].
%
% SUPPOSE x_min, x_max, y_min, y_max ARE THE MINIMUM AND MAXIMUM VALUES IN
% x AND y REACHED BY THE POLYGON, I.E. "R=[x_min,x_max] x [y_min,y_max]" 
% IS THE SMALLEST RECTANGLE WITH SIDES PARALLEL TO THE AXIS x AND y, 
% CONTAINING THE POLYGON.
%
% [rotation=0]: GENERAL CASE, BUT OBSERVE THAT THE INTEGRAND MUST BE 
%               DEFINED IN THE RECTANGLE "R" DESCRIBED ABOVE.
% [rotation=1]: GOOD FOR CONVEX POLYGONS. FOR CONVEX POLYGONS IT SUFFICES 
%               THAT THE INTEGRAND BE DEFINED IN THE POLYGON.
% [rotation=2]: THE USER CHOOSES A SPECIAL REFERENCE SEGMENT "PQ". CHECK   
%               [1] FOR FURTHER DETAILS. SEE THE VARIABLES "P", "Q" BELOW.
%--------------------------------------------------------------------------
rotation=0;

%--------------------------------------------------------------------------
% IF [rotation = 2] THEN THE ALGORITHM CHOOSES A PREFERRED SEGMENT "PQ" 
% HAVING "P", "Q" AS EXTREMA.
%--------------------------------------------------------------------------
P=[0 0.5]; Q=[1 0.5];

%--------------------------------------------------------------------------
% [plot_polygon = 0] NO PLOTS.
% [plot_polygon = 1] IT PLOTS THE DOMAIN AND THE CUBATURE NODES.                   
%--------------------------------------------------------------------------
plot_polygon = 1;

%--------------------------------------------------------------------------
% [cubature_type]: IT CHOOSES THE 1D QUADRATURE RULE NEEDED TO GENERATE A 
%                  QUADRATURE RULE ON THE CURVILINEAR POLYGON.
%
%                  EXAMPLES:
%                  "cubature_type=1": FEJER 1.
%                  "cubature_type=2": FEJER 2.
%                  "cubature_type=3": CLENSHAW CURTIS (VIA WEIGHTS).
%                  "cubature_type=4": GAUSS-LEGENDRE.
%--------------------------------------------------------------------------
cubature_type=4;

%--------------------------------------------------------------------------
% [is_octave]: CHOOSING WORKSPACE.
%           [is_octave=0]: MATLAB.
%           [is_octave=1]: OCTAVE.
%--------------------------------------------------------------------------
is_octave=0;

%--------------------------------------------------------------------------
% [exact_result]: IF THE EXACT RESULT IS KNOWN IT IS INSERTED HERE,
%                 OTHERWISE SET "exact_result=[];"
%                 IT WILL ALLOW TO TEST ABSOLUTE AND RELATIVE ERRORS.
%--------------------------------------------------------------------------
% exact_result=[1.883541031230328];
%0.203076269853
%--------------------------------------------------------------------------
%                    THE MAIN PROGRAM STARTS HERE.
%--------------------------------------------------------------------------

switch spline_type
case 1
    SPLtypestring='complete';
case 2
    SPLtypestring='not-a-knot';
case 3
    SPLtypestring='periodic';
case 4
    SPLtypestring='second';
otherwise
    SPLtypestring='variational';
end

%--------------------------------------------------------------------------
% DEFINE POLYGON.
%--------------------------------------------------------------------------
% "define_polygon.m" DEFINES THE COORDINATES OF THE VERTICES (REPEATING AT
% THE END THE FIRST POINT). THE VARIABLE "polygon_type" IS DEFINED ABOVE. 
% "vertices" IS A VECTOR HAVING AS COMPONENTS THE VERTICES OF THE 
% COMPONENT, DESCRIBED COUNTERCLOCKWISE. 
%--------------------------------------------------------------------------

vertices=define_polygon(polygon_type);

%--------------------------------------------------------------------------
% CORRECTION ON THE NUMBER OF POLYGON SIDES OR ORDER, IF NECESSARY.
%--------------------------------------------------------------------------
L=size(spline_order_vett,1);
Lps=size(vertices,1);
if spline_order_vett(L,2) ~= Lps
    fprintf('\n \t [WARNING]: spline_order_vett IS BADLY DEFINED');
    fprintf('\n \t            CHANGING NUMBER OF VERTICES;');
    fprintf('\n \t            USING LINEAR PARAMETRIC SPLINES');
    spline_order_vett=[2 Lps]; 
end

%--------------------------------------------------------------------------
% COMPUTING INTEGRAL 
%--------------------------------------------------------------------------
% "cubature_result", NODES "(x_nodes,y_nodes)" AND WEIGHTS "weights".
% "fct2D.m" IS AN m-file DEFINED BY THE USER, HAVING "function_type" AS
% VARIABLE THAT SWITCHES BETWEEN FUNCTIONS. SEE THE FILE "fct2D.m" FOR 
% DETAILS.
%--------------------------------------------------------------------------

process_cputime_2009(1)=cputime;

if (cubature_type == 4)
    N_2009=2*N-1;
end

[x_nodes_2009,y_nodes_2009,weights_2009]=splinegauss_2009b(N_2009,... 
    vertices,rotation,P,Q,spline_order_vett,cumulative,SPLtypestring,...
    cubature_type);

f_nodes_2009=fct2D(x_nodes_2009,y_nodes_2009,function_type);
cubature_result_2009=weights_2009'*f_nodes_2009;

   process_cputime_2009(2)=cputime;

%--------------------------------------------------------------------------
% SOME STATS.
%--------------------------------------------------------------------------

% GENERAL PARAMETERS.
fprintf('\n \t ---------------------------------------------------------');
fprintf('\n \t [SIDES POLYGON]: %5.0f [CUBATURE TYPE]: %1.0f',...
    length(vertices(:,1))-1,cubature_type);

fprintf('\n \t [PTS 1D RULE]: %5.0f [PTS]: %5.0f',...
    N,size(x_nodes_2009,1)*size(x_nodes_2009,2));

fprintf('\n \t [POLYGON TYPE ]: %5.0f [FUNCTION]: %2.0f',...
    polygon_type,function_type);

fprintf('\n \t [CUMULATIVE]: %1.0f', cumulative);

% DESCRIBING CURVILINEAR SIDES.

fprintf('\n \t ---------------------------------------------------------');
for index=1:size(spline_order_vett,1)
    if spline_order_vett(index,1) == 4
        fprintf('\n \t [FINAL PT]: %5.0f [ORDER]: %5.0f ', ...
            spline_order_vett(index,2),spline_order_vett(index,1) ); 
        SPLtypestringext=upper(strcat('[END CONDS] [',SPLtypestring,']'));
        fprintf(SPLtypestringext);
    else
        fprintf('\n \t [FINAL PT]: %5.0f [ORDER]: %5.0f',...
            spline_order_vett(index,2),spline_order_vett(index,1) ); 
    end
end

% DESCRIBING DOMAIN ROTATION.
fprintf('\n \t ---------------------------------------------------------');

fprintf('\n \n \t >> SPLINEGAUSS 2009.');
fprintf('\n \t ---------------------------------------------------------');

if rotation == 2
    fprintf('\n \t [ROTATION] [P]:[%1.1e,%1.1e] [Q]:[%1.1e,%1.1e]',...
        rotation,P(1),P(2),Q(1),Q(2));
else
    fprintf('\n \t [ROTATION TYPE]: %1.0f',rotation);
end

number_pts_2009=length(weights_2009); 

fprintf('\n \t [PTS]: %1.0f',number_pts_2009);

% DESCRIBING WEIGHTS.
sw_2009(1)=sum(abs(weights_2009)); 
sw_2009(2)=sum(weights_2009); 
sw_2009(3)=sum(abs(weights_2009))-sum(weights_2009);
neg_w_2009=length(find(weights_2009 < 0));
pos_w_2009=length(find(weights_2009 > 0));
fprintf('\n \t [SUM WEIGHTS ABS. VALUES ]: %2.2e ',sw_2009(1)); 
fprintf('\n \t [SUM WEIGHTS]: %2.2e [SUM -]: %2.2e',sw_2009(2),...
    sw_2009(3)/2);
fprintf('\n \t [WEIGHTS] [NEG.]: %2.0f [POS.]: %2.0f [ALL]: %2.0f',...
    neg_w_2009,pos_w_2009,length(weights_2009)); 
fprintf('\n \t [WEIGHTS] [MAX.]: %2.2e [MIN.]: %2.2e',...
    max(weights_2009),min(weights_2009)); 
fprintf('\n \t ---------------------------------------------------------');
fprintf('\n \t [POLYGSS RES.]: %5.25f',cubature_result_2009);
if length(exact_result) > 0
    abs_err_2009=abs(cubature_result_2009-exact_result);
    fprintf('\n \t [ABS. ERR.]: %2.2e',abs_err_2009); 
    abs_exr_2009=abs(exact_result);
    if abs_exr_2009 > 0
        rel_err_2009=abs_err_2009/abs_exr_2009;
        fprintf(' [REL. ERR.]: %2.2e',rel_err_2009);  
    end
end
fprintf('\n \t ---------------------------------------------------------');


% CHECKING FORMULA (8), (12).

sov2=spline_order_vett(:,2); 
sov2=[1; sov2];
for index=1:size(spline_order_vett,1)
   p_i=spline_order_vett(index,1)-1;
   n_loc=N;
   ni_loc=n_loc*p_i+floor((p_i+1)/2);
   mi=sov2(index+1)-sov2(index)+1;
   loc_points(index)=n_loc*(mi-1)*ni_loc;
end

tot_pts=sum(loc_points);
fprintf('\n \n \t >> OTHER STATS.');
fprintf('\n \t ---------------------------------------------------------');
fprintf('\n \t [EXACT NUMBER OF POINTS]: %7.0f',tot_pts);

fprintf('\n \t [CPUTIME-2009]: %2.2e',diff(process_cputime_2009));
fprintf('\n \t ---------------------------------------------------------');

fprintf('\n');

%--------------------------------------------------------------------------
% PLOT CURVILINEAR POLYGON.
%--------------------------------------------------------------------------
if plot_polygon == 1  
    close all; clf;
    
    spline_blocks_fin=spline_order_vett(:,2);
    if length(spline_blocks_fin) > 1
        sbf=spline_blocks_fin(1:length(spline_blocks_fin)-1);
        add_spline_blocks_fin=sbf;
    else
        add_spline_blocks_fin=[];
    end
    spline_blocks_init=[1; add_spline_blocks_fin];
    
    % PLOTTING EACH CURVILINEAR SIDE.
    
    for index_block=1:length(spline_blocks_fin)
        
        spline_order=spline_order_vett(index_block,1);
        index_start=spline_blocks_init(index_block);
        index_end=spline_blocks_fin(index_block);
        
        local_sides=(index_start:index_end)';
        vertices_loc=vertices(local_sides,:);
        
        %------------------------------------
        % PARAMETRIC SPLINE: PARAMETRIZATION. 
        %------------------------------------
        x_loc=vertices_loc(:,1);
        y_loc=vertices_loc(:,2);
        
        s_loc=[];
        if cumulative == 0
            s_loc=(1:length(vertices_loc(:,1)))';
        else
            s_loc(1,1)=0;
            for index=1:(length(x_loc)-1)
                s_loc(index+1,1)=s_loc(index,1)+...
                    sqrt( (x_loc(index+1)-x_loc(index)).^2 + ...
                    (y_loc(index+1)-y_loc(index)).^2 );
            end
        end
             
        %-----------------------------------
        % SPLINE CURVILINEAR BOUNDARY PTS.
        %-----------------------------------
        s_loc_min=min(s_loc); s_loc_max=max(s_loc);
        Nss=1000;
        ss_loc=s_loc_min:(s_loc_max-s_loc_min)/Nss:s_loc_max;
        
        switch spline_order 
        case 2
            x_border=x_loc;
            y_border=y_loc;
        case 4 
 
            %--------------------------------------------------------------
            % CUBIC SPLINES BY CSAPE. ONE CAN ADD ADDITIONAL REQUESTS 
            % AS THE TYPE OF CUBIC SPLINE ('complete','not-a-knot',
            % 'periodic','second','variational').
            % AS DEFAULT WE USED PERIODIC CUBIC SPLINES.
            %--------------------------------------------------------------
            ppx=csape(s_loc,x_loc,SPLtypestring);
            ppy=csape(s_loc,y_loc,SPLtypestring);
            x_border=ppval(ppx,ss_loc);
            y_border=ppval(ppy,ss_loc);
        otherwise   
            %-----------------------------
            % K-th ORDER SPLINES BY SPAPI.
            %-----------------------------
            ppx=spapi(spline_order,s_loc,x_loc);
            ppy=spapi(spline_order,s_loc,y_loc);
            x_border=fnval(ppx,ss_loc);
            y_border=fnval(ppy,ss_loc);
        end
        
        plot(x_loc,y_loc,'go'); hold on;  % POLYGON VERTICES.
        plot(x_border,y_border,'g-'); hold on;  % BOUNDARY POINTS.
        plot(x_nodes_2009,y_nodes_2009,'m.'); hold on; % POLYGAUSS NODES.
        
    end 
end

axis tight;
