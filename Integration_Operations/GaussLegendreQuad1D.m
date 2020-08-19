function [ Xq,Wq,pExact ] = GaussLegendreQuad1D( N,a,b )
%function [ Xq,Wq,pExact ] = GaussLegendreQuad1D( N,a,b ) generates the
%Gaussian quadrature scheme based on the Legendre nodes containing N points
%and weights which is accurate up to degree p = 2*N-1. The interval of 
%integration is [a,b]. If a and be are not specified by the user, the 
%default interval is [0,1].
%
%   Inputs: N - Number of quadrature points and weights desired
%           a - left endpoint of integration interval (default a = 0)
%           b - right endpoint of integration interval (default b = 1)
%
%   Outputs: Xq - list of N quadrature nodes
%            Wq - list of N quadrature weights
%            pExact - polynomial exactness of quadrature scheme (p = 2*N-1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing code
%
% N = 3;
%
% [ Xq,Wq,pExact ] = GaussLegendreQuad1D( N );
%
% plot(Xq,0,'ro')
%
% f =@(x) x.^pExact;
%
% x = linspace(0,1);
% ExactInt = 1/(pExact+1);
% QuadInt = Wq*f(Xq);
%
% fprintf('function: x^%d \n', pExact)
% fprintf('Exact Integral: %e \n',ExactInt);
% fprintf('Approx. Integral: %e \n',QuadInt);
% fprintf('Error: %e \n',abs(ExactInt-QuadInt));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Based on the algorithm introduced in:
%   Golub, Gene H., and John H. Welsch. "Calculation of Gauss quadrature rules."
%   Mathematics of computation 23.106 (1969): 221-230.
%
% Written by: Joseph Benzaken 
% Last Updated: 5/4/16
% 
% CMGLab

% if no specified interval of integration, set to default of [0,1]
if nargin == 1
    a = 0;
    b = 1;
end

%initialize arrays to store reccurence coefficients
alpha = zeros(N,1);
beta = zeros(N,1);

%Legendre 3-term recursion relationship
for i=1:N
    alpha(i) = 0;
    beta(i) = (i-1)^2/(4*(i-1)^2-1);
end

%assemble Jacobi matrix
Tjacobi = diag(alpha) + sqrt(diag(beta(2:end),1)) + sqrt(diag(beta(2:end),-1));

[U,V]=eig(Tjacobi); %compute eigenvalues and eigenvectors of Jacobi matrix

Xq = toStandardInterval(diag(V),-1,1,a,b); %quadrature nodes are eigenvalues
Wq = (b-a)*U(1,:).^2; %quadrature weights are first component of eigenvectors squared and scaled by appropriate interval
pExact = 2*N-1; %polynomial exactness is 2*N-1
end

function [ stdX ] = toStandardInterval( X,oldMin,oldMax,newMin,newMax )
%function [ stdX ] = toStandardInterval( X,oldMin,oldMax,newMin,newMax )
%   Linearly maps X in [oldMin,oldMax] to stdX in [newMin,newMax]
stdX = (newMax-newMin)/(oldMax-oldMin)*(X-oldMin)+newMin;

end