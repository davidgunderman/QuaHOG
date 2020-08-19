function reversemat = rorient(Element)
numRows=size(Element,1);
for i=1:3:numRows
    Element2(i:(i+2),:)=fliplr(Element((numRows-(i+1)):(numRows-i+1),:));
end
reversemat=Element2;
end