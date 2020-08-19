% Example 3: Numerical approximation of an integral on the interval [0,1]  
%            where the integrand has one multiple complex conjugate pair
%            of singularities
addpath("../Src","../Tests","../ThirdPartySupportingCode")
clear; 
close; 
clc; 

disp('One multiple complex conjugate pair of poles \pm i*B ');
disp('int_{0}^{1} sinh( 1/(x^2+B^2) ) dx');
fun = @(x) ((1-sqrt(2))*x.^2 + (sqrt(2)-2)*x+1)./((2-sqrt(2))*x.^2+(sqrt(2)-2)*x+1);
% B = [1:1:10]; % ten different cases

% Exact = [0.87974890970966687009990261325026,0.23393310806757554569295582399079, ...
%  0.10745651578818868161254072898437,0.061282997467641006451117098989414, ... 
%  0.03948937231856718098266526520331,0.027528255918594052999339126939024, ...
%  0.020272396252557931274876024669302,0.015545000357541637593303195735998, ...
%  0.012295556596944788759891908194807,0.0099670302697005460949357802401303];
for jj=1:10000
sgl=repmat(roots([(2-sqrt(2)),(sqrt(2)-2),1]),2,1);
%     disp(sprintf('\nPole is +/- %1.3fi',B(j)));
    % Whenever the complex conjugate counterpart of a complex pole is 
    % missing in the sequence of poles, it will automatically be added 
    % to the sequence by the program rfejer itself.    
%     n = 16;
%     disp(sprintf('The maximal number of iterations is %1.3f',n+1));
%     sgl = 1i*B(j)*ones(1,n);
    % First we need to map the interval [0,1] onto the interval [-1,1]
    [fout,sglout ] = transf( @(x)fun(x) , sgl , [0,1] );
    x = rfejer(sglout);
    NumInt = fout(x(1,:))*x(2,:)';
%     disp(sprintf('Computed value: %1.16e',NumInt));
%     disp(sprintf('Estimated relative error : %1.16e',Err));
    ErrExact = abs(1-NumInt/0.6336517494806966);
end
%     disp(sprintf('Exact relative error : %1.16e',ErrExact));
