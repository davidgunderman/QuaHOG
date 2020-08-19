function shape= generateTestFigures(figNum)
% This file outputs one of the 8 test figures that are used in "Spectral 
% mesh-free quadrature for planar regions bounded by rational parametric 
% curves"
%
% figNum ==1 outputs a square plate with a circular hole
% figNum ==2 outputs an L-bracket with 3 holes
% figNum ==3 outputs a wrench figure
% figNum ==4 outputs a guitar-shaped object
% figNum ==5 outputs a treble clef
% figNum ==6 outputs a difficult-to-mesh rotor
% figNum ==7 outputs a gear shape
% figNum ==8 outputs a circle
% Each of the outputs have a corresponding NURBS file
%
% The format of the output is a cellArray. Each cell contains a 
% connected, oriented boundary component.
% 
% Each cell is a rational bezier curve with dimensions
% 3n X (m+1) where n is the number of curve components and m is the degree
% of the curve. Every three rows comprises one rational bezier curve, where
% the first row is the x component, second row is the y component, and
% third curve is the weights

% Circle primitive for creating objects 1, 2, 7 and 8
Circletemp=load("circle_rational.mat"); seedFigure1= Circletemp.Circle1;
Squaretemp=load("square_rational.mat"); square=Squaretemp.square;

if figNum==1
    plate= offset(rorient(square),[1;1]);
    hole = offset(stretch(rorient(seedFigure1),[.5;.5]),[1;1]);
    shape{1}=plate;
    shape{2}=hole;
end

if figNum==2
    plate= offset(stretch(rorient(square),[1.5;2]),[.5;0]);
    minusplate= offset(stretch(square,[1.5;2]),[1.5;-1]);
    middleplate= offset(stretch(rorient(square),[.5;.5]),[0;1]);
    middlecircle= offset(stretch(rorient(rotate(seedFigure1,.5)),[.5;.5]),[.5;.5]);
    [~,ss]= RatboolEls(middleplate,middlecircle,1);
    [~,lbracket]=RatboolEls(plate,minusplate,1);
    lbracket{1}(2,:)=[-2 -.25 .5];
    lbracket{1}(4,:)=[.5 1.25 2];
    middlecirclecut=ss{1}(1:6,:);
    culbracket=[lbracket{1}(1:3,:); middlecirclecut; lbracket{1}(4:end,:)];
    hole = stretch(rorient(seedFigure1),[.25;.25]);
    hole1 = offset(hole,[-.5;1.5]);
    hole2 = offset(hole,[1.5;1.5]);
    hole3 = offset(hole,[-.5;-1.5]);
    shape={culbracket,hole1,hole2,hole3};
end

if figNum==3
    load guitar_rational
      shape=guitar;
end

if figNum==4
    load treble_clef_rational;
    shape=clef;
end

if figNum==5
    load rotor_rational;
end
      
if figNum==6
    load wrench_rational;
end

if figNum==7
    lroffset=[.3957106819596820 -1.5728593603867]; Circle= seedFigure1;
    Circle(1:3:end,:)=Circle(1:3:end,:)+lroffset(1).*Circle(3:3:end,:);
    Circle(2:3:end,:)=Circle(2:3:end,:)+lroffset(2).*Circle(3:3:end,:);
    shape{1}=Circle;
end

if figNum==8
    Octagon=[1          (2+sqrt(2))/4   sqrt(2)/2;
            0           sqrt(2)/4       sqrt(2)/2;
            1           1               1;
            sqrt(2)/2   sqrt(2)/4       0;
            sqrt(2)/2   (2+sqrt(2))/4   1;
            1           1               1;
            0           -sqrt(2)/4      -sqrt(2)/2;
            1           (2+sqrt(2))/4   sqrt(2)/2;
            1           1               1;
            -sqrt(2)/2  -(2+sqrt(2))/4  -1;
            sqrt(2)/2   sqrt(2)/4       0;
            1           1               1;
            -1          -(2+sqrt(2))/4  -sqrt(2)/2;
            0           -sqrt(2)/4      -sqrt(2)/2;
            1           1               1;
            -sqrt(2)/2  -sqrt(2)/4      -0;
            -sqrt(2)/2  -(2+sqrt(2))/4  -1;
            1           1               1;
            0           sqrt(2)/4      sqrt(2)/2;
            -1          -(2+sqrt(2))/4  -sqrt(2)/2;
            1           1               1;
            sqrt(2)/2   (2+sqrt(2))/4   1;
            -sqrt(2)/2  -sqrt(2)/4       0;
            1           1               1;
            ];
            Circle= seedFigure1;
            Notch=stretch(Circle,[.25,.25]);
            Centers=Octagon(:,1); Centers(3:3:end)=[];
            Centers=reshape(Centers,2,8)';
            ss{1}=Octagon;
            for ii=1:8
                [~,sstemp]=RatboolEls(ss{ii},rorient(offset(Notch,Centers(ii,:)')),1);
                ss{ii+1}= sstemp{1};
            end
            clear sstemp
            shape{1}=stretch(square,[2 2]);
            shape{2}=rorient(ss{9});
end
end


% Geometric operations (rotate, shift, stretch, reverseorient)
function rotatemat = rotate(Element,theta)
RotateMat= [cos(theta), -sin(theta); sin(theta), cos(theta)];
rotatemat=Element;
for i=1:3:size(Element,1)
    rotatemat(i:(i+1),:)=RotateMat*Element(i:(i+1),:);
end
end

function offsetmat = offset(Element,xy)
for i=1:3:size(Element,1)
        Element(i:(i+1),:)=Element(i:(i+1),:)+xy.*Element(i+2,:);
end
offsetmat=Element;
end

function stretchmat = stretch(Element,xy)
for i=1:3:size(Element,1)
        Element(i,:)=xy(1).*Element(i,:);
        Element(i+1,:)=xy(2).*Element(i+1,:);
end
stretchmat=Element;
end

function reversemat = rorient(Element)
numRows=size(Element,1);
for i=1:3:numRows
    Element2(i:(i+2),:)=fliplr(Element((numRows-(i+1)):(numRows-i+1),:));
end
reversemat=Element2;
end