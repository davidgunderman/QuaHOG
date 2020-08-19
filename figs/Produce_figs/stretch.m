function stretchmat = stretch(Element,xy)
for i=1:3:size(Element,1)
        Element(i,:)=xy(1).*Element(i,:);
        Element(i+1,:)=xy(2).*Element(i+1,:);
end
stretchmat=Element;
end