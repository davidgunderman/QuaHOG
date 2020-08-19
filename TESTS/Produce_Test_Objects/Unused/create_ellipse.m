function ellipseCPWS = create_ellipse(c0,c1,theta,r0,r1)
RotateMat= [cos(theta), -sin(theta); sin(theta), cos(theta)];
Circletemp=load("Circle.mat"); Circle= Circletemp.Circle1;
for i=1:3:size(Circle,1)
    Circle(i,:)=Circle(i,:).*r0;
    Circle(i+1,:)=Circle(i+1,:).*r1;
    Circle(i:(i+1),:)=RotateMat*Circle(i:(i+1),:);
    Circle(i,:)=Circle(i,:)+c0.*Circle(i+2,:);
    Circle(i+1,:)=Circle(i+1,:)+c1.*Circle(i+2,:);
end
ellipseCPWS=Circle;
end