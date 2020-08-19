
%--------------------------------------------------------------------------
% REFERENCE PAPERS.
%-------------------
%
% [1]. A. SOMMARIVA and M. VIANELLO "Gauss-like and triangulation-free 
%      cubature over polygons". BIT Numerical Methematics 47 (2007),
%      441-453.
%
% [2]. A. SOMMARIVA and M. VIANELLO "Gauss-Green cubature an moment 
%      computation over arbitrary geometries". JCAM.
%
%--------------------------------------------------------------------------
%
% THIS DEMO SHOWS HOW TO USE THE PROCEDURE "splinegauss" FOR COMPUTING
% CUBATURE OF CONTINUOUS FUNCTIONS OVER CERTAIN CURVILINEAR POLYGONS.
% IN PARTICULAR IT COMPUTES THE MOMENTS OF LEGENDRE OR MONOMIAL OR 
% CHEBYSHEV BASES (OF A FIXED DEGREE basis_degree) ON SOME SPLINE 
% CURVILINEAR POLYGONS.
%--------------------------------------------------------------------------
% FUNCTIONS USED.
%
% 1. define_polygon
% 2. splinegauss_2009b
% 3. legendre_poly
% 4. vander_poly
% 5. chebyshev_poly
%--------------------------------------------------------------------------
%% Copyright (C) 2007 Marco Vianello and Alvise Sommariva
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
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%--------------------------------------------------------------------------
%% Authors:  
%% Marco Vianello    <marcov@euler.math.unipd.it>
%% Alvise Sommariva  <alvise@euler.math.unipd.it>
%% Date: NOVEMBER 19, 2009
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
%
% TO DEFINE YOUR OWN POLYGON, PLEASE OPEN THE m-FILE "define_polygon.m" 
% AND ADD IT TO THE LIST. OBSERVE THAT THE FIRST VERTEX IS REPEATED. 
% FOR EXAMPLE, THE UNIT SQUARE IS DEFINED AS "[0 0; 1 0; 1 1; 0 1; 0 0]"
% AND NOT "[0 0; 1 0; 1 1; 0 1]".
%--------------------------------------------------------------------------
polygon_type=14;

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
spline_order_vett=[4 37; 2 38]; 

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
spline_type=2;

%--------------------------------------------------------------------------
% [cumulative]: THIS VARIABLE CHOOSES THE PARAMETRIZATION. 
%              
%              EXAMPLES:
%              "cumulative=0": 1:N.
%              "cumulative=1": CUMULATIVE.
%--------------------------------------------------------------------------
cumulative=1;

%--------------------------------------------------------------------------
% MOMENTS ARE COMPUTED UP TO DEGREE "basis_degree".
%--------------------------------------------------------------------------
basis_degree=16;                                  

%--------------------------------------------------------------------------
% WRITE A VECTOR OF 2 COMPONENTS CONTAINING THE NUMBER OF POINTS OF TWO 1D
% QUADRATURE RULES USED BY "splinegauss".
% 
% FIRST COMPONENT.  IT IS THE 1D RULE USED AS REFERENCE RESULT FOR MAKING
%                   THE COMPARISON.
%
% SECOND COMPONENT. IT IS THE 1D RULE THAT THE USER INTENDS TO USE FOR 
%                   COMPUTING THE MOMENTS.
%. 
% THEY CORRESPOND TO COMPUTE A CUBATURE RULE ON THE DOMAIN HAVING 
% ALGEBRAIC DEGREES OF EXACTNESS (ADE)  "2*N-1".
% 
% EXAMPLE:
%
% SUPPOSE basis_degree=16
% N_vett=[ceil((basis_degree+1)/2) 7] WILL COMPUTE TWO CUBATURE RULES.
% THE FIRST ONE, USED AS REFERENCE, HAS ADE=17.
% THE SECOND ONE HAS ADE=13 AND IT IS THE CUBATURE RULES CHOOSEN BY THE 
% USER.
%--------------------------------------------------------------------------
N_vett=[ceil((basis_degree+1)/2) 10];

%--------------------------------------------------------------------------
% [basis_type]: BASIS TYPE IN MOMENTS CALCULATIONS. 
%               1. LEGENDRE.
%               2. MONOMIAL.
%               3. CHEBYSHEV.
%--------------------------------------------------------------------------
basis_type=1;                   

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
rotation=2; 

%--------------------------------------------------------------------------
% IF [rotation = 2] THEN THE ALGORITHM CHOOSES A PREFERRED SEGMENT "PQ" 
% HAVING "P", "Q" AS EXTREMA.
% example: 
% P=[0.43 0.68]; Q=[0.43 0.236];
%--------------------------------------------------------------------------
% P=[0 0.5]; Q=[1 0.5];
P=[0.43 0.68]; Q=[0.43 0.236];

%--------------------------------------------------------------------------
% [plot_polygon = 0] NO PLOTS.
% [plot_polygon = 1] IT PLOTS THE DOMAIN AND THE CUBATURE NODES.                   
%--------------------------------------------------------------------------
plot_polygon = 1;

%--------------------------------------------------------------------------
% [cubature_type]: IT CHOOSES THE 1D QUADRATURE RULE .
%
%           [cubature_type=1]: FEJER 1.
%           [cubature_type=2]: FEJER 2.
%           [cubature_type=3]: CLENSHAW CURTIS (VIA WEIGHTS).
%           [cubature_type=4]: GAUSS-LEGENDRE.
%
% PS: WE SUGGEST TO USE "cubature_type=4".
%--------------------------------------------------------------------------
cubature_type=4;

%--------------------------------------------------------------------------
% [is_octave]: CHOOSING WORKSPACE.
%           [is_octave=0]: MATLAB.
%           [is_octave=1]: OCTAVE.
%--------------------------------------------------------------------------
is_octave=0;

%--------------------------------------------------------------------------
% [moments_method]: METHODS USED FOR MOMENTS COMPUTATION.
%           [moments_method=0]: VIA WEIGHTED VANDERMONDE.
%           [moments_method=1]: VIA VANDERMONDE.
%--------------------------------------------------------------------------
moments_method=0;

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

if rotation < 2
    P_vett=[0 0]; Q_vett=[0 1];
end



%--------------------------------------------------------------------------
% DEFINE POLYGON.
%--------------------------------------------------------------------------
% "define_polygon.m" DEFINES THE COORDINATES OF THE VERTICES (REPEATING AT
% THE END THE FIRST POINT). THE VARIABLE "polygon_type" IS DEFINED ABOVE. 
% "polygon_sides" IS A VECTOR HAVING AS COMPONENTS THE VERTICES OF THE 
% COMPONENT, DESCRIBED COUNTERCLOCKWISE. 
%--------------------------------------------------------------------------

if is_octave == 1 % OCTAVE DOES NOT IMPLEMENT SPAPI AND GINPUT. PATCH.
    
    if polygon_type == 9
        err_str='[ERROR]: GINPUT NOT IMPLEMENTED IN OCTAVE.';
        fprintf(err_str);
        fprintf('\n \t TRY ANOTHER POLYGON, BY DIGITING ANOTHER NUMBER');
        polygon_type=input('\n \t INTEGER BETWEEN 1 AND 14 (BUT NOT 9).');
        vertices=define_polygon(polygon_type);
    else
        vertices=define_polygon(polygon_type);
    end
    
    L=size(spline_order_vett,1);
    Lps=size(vertices,1);
    if spline_order_vett(L,2) ~= Lps
        fprintf('\n \t [WARNING]: spline_order_vett IS BADLY DEFINED. ');
        err_str1='CHANGING NUMBER OF VERTICES. USING LINEAR PARAM.SPLINES';
        fprintf(err_str1);
        spline_order_vett=[2 Lps]; 
    end
    
    for index=1:size(spline_order_vett,1)
        sov=spline_order_vett(index,1);
        if (sov ~= 2) & (sov ~= 4)
            errstr='OCTAVE DOES NOT IMPLEMENT SPAPI. USING CSAPE.';
            fprintf(errstr);
            spline_order_vett(1,1)=4;
        end
    end  
    
else
    
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
    
end




%--------------------------------------------------------------------------
% COMPUTING INTEGRAL 
%--------------------------------------------------------------------------
% "cubature_result", NODES "(x_nodes,y_nodes)" AND WEIGHTS "weights".
% "fct2D.m" IS AN m-file DEFINED BY THE USER, HAVING "function_type" AS
% VARIABLE THAT SWITCHES BETWEEN FUNCTIONS. SEE THE FILE "fct2D.m" FOR 
% DETAILS.
%--------------------------------------------------------------------------
for index=1:2
    
    N=N_vett(index);
    
    [x_nodes,y_nodes,weights]=splinegauss_2009b(N_2009,... 
    vertices,rotation,P,Q,spline_order_vett,cumulative,SPLtypestring,...
    cubature_type);
    
    sqrt_weights=sqrt(weights);
    
    xmax=max(vertices(:,1)); 
    xmin=min(vertices(:,1));
    
    ymax=max(vertices(:,2)); 
    ymin=min(vertices(:,2));
    
    %-------------------------------------------------------------------------
    % EVALUATION OF LEGENDRE OR MONOMIAL OR CHEBYSHEV POLYNOMIALS. 
    % IF REQUIRED, IT MULTIPLIES THEIR VALUES FOR sqrt(weights) SO THAT THE  
    % (WEIGHTED) BIVARIATE VANDERMONDE MATRIX IS EASILY COMPUTED.
    %-------------------------------------------------------------------------
    switch basis_type
    case 1
        % LEGENDRE.
        xt_nodes=(x_nodes-(xmax+xmin)/2)*(2/(xmax-xmin));
        yt_nodes=(y_nodes-(ymax+ymin)/2)*(2/(ymax-ymin));
        [lx_A_LS, lx_DA_LS]=legendre_poly(basis_degree,xt_nodes',sqrt_weights);
        [ly_A_LS, ly_DA_LS]=legendre_poly(basis_degree,yt_nodes',sqrt_weights);
        
    case 2
        
        % CANONICAL.
        lx_DA_LS = vander_poly(basis_degree,x_nodes,sqrt_weights);
        ly_DA_LS = vander_poly(basis_degree,y_nodes,sqrt_weights);
        lx_A_LS = vander_poly(basis_degree,x_nodes,ones(size(sqrt_weights)));
        ly_A_LS = vander_poly(basis_degree,y_nodes,ones(size(sqrt_weights)));
        
    case 3
        
        % CHEBYSHEV.
        xt_nodes=(x_nodes-(xmax+xmin)/2)*(2/(xmax-xmin));
        yt_nodes=(y_nodes-(ymax+ymin)/2)*(2/(ymax-ymin));
        [lx_A_LS, lx_DA_LS]=chebyshev_poly(basis_degree,xt_nodes',sqrt_weights);
        [ly_A_LS, ly_DA_LS]=chebyshev_poly(basis_degree,yt_nodes',sqrt_weights);
        
    end
    
    
    %------------------------------------------------------------------
    % EVERY BIVARIATE POLYNOMIAL IS PRODUCT OF UNIVARIATE POLYNOMIALS:
    % P(x,y)=L_i(x) L_j(y). 
    % HERE WE CHOOSE THE COUPLES (i,j).
    %------------------------------------------------------------------
    A=[];
    for k=1:(basis_degree+1)
        for j=1:k
            Aloc=[j k-j+1];
            A=[A; Aloc];
        end
    end
    
    %------------------------
    % COMPUTING THE MOMENTS.
    %------------------------
    if moments_method == 1
        LLX_DA_LS=lx_DA_LS(A(:,1),:);
        LLY_DA_LS=ly_DA_LS(A(:,2),:);
        DA_LS=(LLX_DA_LS.*LLY_DA_LS)'; 
        b_moms= (sum(DA_LS))';
    else
        LLX_A_LS=lx_A_LS(A(:,1),:);
        LLY_A_LS=ly_A_LS(A(:,2),:);
        A_LS=(LLX_A_LS.*LLY_A_LS)'; 
        b_moms=A_LS'*weights;
    end
    
    %------------------------
    % REFERENCE MOMENTS.
    %------------------------
    if index == 1
        exact_moms=b_moms;
    end
    
end

%--------------------------------------------------------------------------
% STATS.
%--------------------------------------------------------------------------
fprintf('\n \t ----------------------------------------------------------');

fprintf('\n \t [POLYGON]: %3.0f [BASIS DEGREE]: %3.0f',...
    polygon_type,basis_degree);
fprintf('\n \t [POL.DIM]: %4.0f [BASIS TYPE]: %1.0f',...
    (basis_degree+1)*(basis_degree+2)/2,basis_type);
fprintf('\n \t [MOMENTS] INTEGRALS COMPUTED: %6.0f ',length(b_moms));
fprintf('\n \t -----------------------------------------------------------');
index_weights_positive=find(weights > 0);
fprintf('\n \t [WEIGHTS-SPLINEGAUSS] PTS: %10.0f POSITIVE: %10.0f',...
    length(weights),length(index_weights_positive));
fprintf('\n \t [WEIGHTS-SPLINEGAUSS] SUM: %2.2e  ABS.SUM : %2.2e',...
    sum(weights),sum(abs(weights)));
fprintf('\n \t -----------------------------------------------------------');
fprintf('\n \t [MOMENTS] ERROR INF.TY NORM: %2.2e ',...
    norm(b_moms-exact_moms,inf));
fprintf('\n \t -----------------------------------------------------------');
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
        
        plot(x_loc,y_loc,'ko'); hold on;  % POLYGON VERTICES.
        plot(x_border,y_border,'k-'); hold on;  % BOUNDARY POINTS.
        plot(x_nodes,y_nodes,'m.'); hold on; % POLYGAUSS NODES.
        
    end 
end

axis tight;

