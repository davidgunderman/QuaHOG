Circletemp=load("Circle.mat"); seedFigure1= Circletemp.Circle1;
square=[-1 0 1;
             1 1 1;
              1 1 1;
              1 1 1;
              1 0 -1;
              1 1 1;
              1 0 -1
              -1 -1 -1;
              1 1 1;
              -1 -1 -1;
              -1 0 1;
              1 1 1;
              ]
% wrenchhead= stretch(seedFigure1,[1;1.3])
% wrenchnotch= offset(rotate(stretch(square,[.75;.5]),.2),[-.6; 0])
% wrenchhandle= offset(stretch(rorient(square),[2.5; .5]),[1.75;0])
% plot_rat_bern_poly(wrench,2,.001,"b.")
seedFigure2=[ 0 .2381692783 .8167239668 1;
            -.3 .129368 -.41672386 0;
            1 .7671823 .12692386 1;
            1 1.67123846719234 1.283657146 1.3;
            0 .1672361023 .63486123406 1;
            1 .19620376491 .76838194 1;
            1.3 .6386712306  .168235061 -.2;
            1 1.236817234869 .7681234869173 1.1;
            1 .67123856712384 .368129348 1;
            -.2 .1623857691 -.168234769 0;
            1.1 .65182376891324 .12635781234 -.3;
            1 .7682349876 .571846893467 1;
            ]
        
lll=stretch(rotate(offset(seedFigure2,[.3; .3]),.5),[2 1]);

seedFigure3=seedFigure2(:,[1,2,4]);
nowFigs={seedFigure1 seedFigure1};
clear newFig;
for jj=1:10
newFig{1}= seedFigure1; 
for figNum=1:6
    close
    offset1=(.5*(rand(2,1)-1));
    offset2=(1*(rand(2,1)-1));
    figs{1}=stretch(rotate(offset(nowFigs{1},(.5*(rand(2,1)-1))),.3+rand()),((rand(2,1)+.5)));
    figs{2}=stretch(rotate(offset(nowFigs{2},(.5*(rand(2,1)-1))),.3+rand()),((rand(2,1)+.5)));
    choose=floor((rand())+1.5);
    newFig(figNum+1) = RatboolEls(figs{choose},newFig{figNum},0);
    plot_rat_bern_poly({newFig{figNum}},2,.0001,"k",1)
end
testFig{jj}=newFig{4};
end

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
