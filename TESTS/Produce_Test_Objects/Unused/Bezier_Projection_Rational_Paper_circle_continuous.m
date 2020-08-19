% clear all
% close all
% clearvars
% circlebreakpoints{1}=[1,11];
%{2×105 double}    {2×51 double}    {2×19 double}    {2×25 double}
% circlebreakpoints{1}=[1,2,3,4,5];
% circlebreakpoints{2}=[1,5];
circlebreakpoints{1}=[1,15,20,28,41,43,50,53];
circlebreakpoints{2}=[1,2,5,10,17,19,26];
circlebreakpoints{3}=[1,8,11];
circlebreakpoints{4}=[1,13];
% circle=generateTestFigures(5);
circle=generateTestFigures(5);
for CURVE=3:length(circle)
    newcirclecurve=circle{CURVE};
    circle{CURVE}(1:3:end,:)=circle{CURVE}(1:3:end,:)./circle{CURVE}(3:3:end,:);
    circle{CURVE}(2:3:end,:)=circle{CURVE}(2:3:end,:)./circle{CURVE}(3:3:end,:);
    P3=[]; W3=[]; Xi3=[]; W2=[]; P=[]; P2=[]; W=[]; counter=1; ctr=1;
        starter=1;
        ender=size(circle{CURVE},1)/3;
        testShape=circle{CURVE};
        nb=size(testShape,1)/3;
        p=size(testShape,2)-1;
        Xitemp=[zeros(1,p) 0:nb nb*ones(1,p)];
            Xitemp=zeros(1,p);
            for ii=1:(length(circlebreakpoints{CURVE})-1)
                Xitemp=[Xitemp (circlebreakpoints{CURVE}(ii):(circlebreakpoints{CURVE}(ii+1))-1)-1 repmat((circlebreakpoints{CURVE}(ii+1)),1,p-1)-1];
            end
%             Xitemp(end-(p-1):end)=[];
            Xitemp=[Xitemp nb*ones(1,p)];
            Xi2=Xitemp;
%             Xi2=[zeros(1,p) Xitemp(:)' nb*ones(1,p)];
        if starter==ender
            P=testShape(1:2,:);
            W2=testShape(3,:);
            figure(1);
            plotNURBScurve(2,p,Xi2,P(1:2,:),W2)
            hold on
        else
            n=length(Xi2)-p-1;
            [n_el,C,IEN2] = Extract_Basis_1D(p,n,Xi2);
            P2{nb}=0;
            W{nb}=0;
            P2(1:nb)={zeros(2,n)};
            W(1:nb)={zeros(1,n)};
            P=zeros(2,n);
            W2=zeros(1,n);
            figure(1);
            for i=1:nb
                plotNURBScurve(2,p,[zeros(1, p+1) ones(1, p+1)],testShape(((i-1)*3+1):(3*i-1),:),testShape((3*i),:))
                hold on
                P2{i}(:,IEN2(:,i))=testShape(((i-1)*3+1):(3*i-1),:)/C(:,((i-1)*(p+1))+1:(i*(p+1)));
                W{i}(:,IEN2(:,i))=testShape(3*i,:)/C(:,((i-1)*(p+1))+1:(i*(p+1)));
                scatter(P2{i}(1,IEN2(:,i)),P2{i}(2,IEN2(:,i)));
                hold on
            end
            for j=1:n
                counter=0;
                for k=1:nb
                    if j==22
                        j;
                    end
                    if any(P2{k}(:,j)~=0)
                        P(:,j)=P(:,j)+(P2{k}(:,j));
                        W2(j)=W2(j)+(W{k}(j));
                        counter=counter+1;
                    end
                end
                if counter==0
                    counter=1;
                end
               P(:,j)=P(:,j)/counter;
               W2(j)=W2(j)/counter;
               if W2(j)==0
                   W2(j)=1;
               end
            end
        end
%         P=P./W2;
        perturb=3*rand(1,length(W2)-2);
        perturb=[0];
        for ii=1:(length(circlebreakpoints{CURVE})-1)
            perturb=[perturb 3*rand(1,circlebreakpoints{CURVE}(ii+1)-(circlebreakpoints{CURVE}(ii))) 0];
        end
%         W2(2:(end-1))=W2(2:(end-1))+3*rand(1,length(W2)-2);
%         P=P.*W2;
        P3(:,(ctr):(ctr+length(W2)-1))=P(:,1:(end));
        W3(:,(ctr):(ctr+length(W2)-1))=W2(1:(end));
        ctr=ctr+length(W2);
    newcircle{CURVE}=testShape;
    CP{CURVE}=P3;
    W{CURVE}=W3;
    Xi{CURVE}=Xi2;
end
% bkpts{1}= [5];
for i=1:length(circle)
    bkpts{i}=circlebreakpoints{i}(2:end);
    bkpts{i}(end)=bkpts{i}(end)-1;
end
    