function offsetmat = offset(Element,xy)
for i=1:3:size(Element,1)
        Element(i:(i+1),:)=Element(i:(i+1),:)+xy.*Element(i+2,:);
end
offsetmat=Element;
end