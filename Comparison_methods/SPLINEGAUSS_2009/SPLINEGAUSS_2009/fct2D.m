function z=fct2D(x,y,function_type)

%-------------------------------------------------------------------------------
% FUNCTION DEFINITION. 
%-------------------------------------------------------------------------------
%
% INPUT:
%
% [x,y]: EVALUATE THE FUNCTION IN THE POINTS (x,y). "x", "y" ARE COLUMN VECTORS.
%
% [function_type]: PARAMETER THAT CHOOSES THE FUNCTION.
%
% OUTPUT:
%
% [z]: VALUE OF THE FUNCTION IN THE POINTS (x,y). IT IS A COLUMN VECTOR.
%
%-------------------------------------------------------------------------------

switch function_type
case 1
    z= .75*exp(-((9*x-2).^2 + (9*y-2).^2)/4) + ...
         .75*exp(-((9*x+1).^2)/49 - (9*y+1)/10) + ...
         .5*exp(-((9*x-7).^2 + (9*y-3).^2)/4) - ...
         .2*exp(-(9*x-4).^2 - (9*y-7).^2);
case 2
    z=( (x-0.5).^2 +(y-0.5).^2 ).^(1/2);
case 3
    k=3;
    z=(x+y).^k;
case 4
    loc_arg=(x-0.5).^2+(y-0.5).^2;
    z=exp(-loc_arg);
case 5
    loc_arg=(x-0.5).^2+(y-0.5).^2;
    z=exp(-100*loc_arg);;
case 6
    degree_loc=20;
    z=cos(degree_loc*(x+y));
case 7
    z=ones(size(x));
case 8 
    z=exp(x+y);
case 9
    den=1+16*(x.^2+y.^2);
    z=1./den;
case 10
    exponent=3/2;
    z=(x.^2+y.^2).^(exponent);
end

