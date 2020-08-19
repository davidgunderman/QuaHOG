clear all
close all
clearvars
rotorbreakpoints{1}=[1,16];
rotorbreakpoints{2}=[1,61];
rotorbreakpoints{3}=[1,54];
rotor=generateTestFigures(6);
for CURVE=1:length(rotor)
    newrotorcurve=rotor{CURVE};
    P3=[]; W3=[]; Xi3=[]; W2=[]; P=[]; P2=[]; W=[]; counter=1; ctr=1;
    for segment=1:(length(rotorbreakpoints{CURVE})-1)
        starter=rotorbreakpoints{CURVE}(segment);
        ender=rotorbreakpoints{CURVE}(segment+1)-1;
        testShape=rotor{CURVE}(((starter-1)*3+1):(ender*3),:);
        nb=size(testShape,1)/3;
        p=size(testShape,2)-1;
    %     Xi=[zeros(1,p) 0:nb nb*ones(1,p)];
            Xitemp=repmat(0:1:nb,p,1);
            Xi2=[0 Xitemp(:)' nb];
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
            end
            for i=1:2
                for j=1:n
                    counter=0;
                    for k=1:nb
                        if P2{k}(i,j)~=0
                            P(i,j)=P(i,j)+(P2{k}(i,j));
                            counter=counter+1;
                        end
                    end
                   P(i,j)=P(i,j)/counter;
                end
            end
            for j=1:n
                counter=0;
                for k=1:nb
                    if P2{k}(i,j)~=0
                        W2(j)=W2(j)+(W{k}(j));
                        counter=counter+1;
                    end
                end
               W2(j)=W2(j)/counter;
            end
        end
        W2(2:(end-1))=W2(2:(end-1))+1*rand(1,length(W2)-2);
        P3(:,(ctr):(ctr+length(W2)-1))=P(:,1:(end));
        W3(:,(ctr):(ctr+length(W2)-1))=W2(1:(end));
        Xitemp=repmat(0:1:((length(W3))/(p+1)),p+1,1);
        Xi3=Xi2;
        ctr=ctr+length(W2);
%         figure(2);
%         plotNURBScurve(2,p,Xi,P(1:2,:),W2)
%         hold on
        newtestShape{segment}=testShape;
%         for i=1:nb
%     %         P=P.*W2;
%             newtestShape{segment}(((i-1)*3+1):(3*i-1),:)=P(:,IEN2(:,i))*C(:,((i-1)*(p+1))+1:(i*(p+1)));
%             newtestShape{segment}(3*i,:)=W2(:,IEN2(:,i))*C(:,((i-1)*(p+1))+1:(i*(p+1)));
%             figure(3)
%             plotNURBScurve(2,p,[zeros(1, p+1) ones(1, p+1)],newtestShape{segment}(((i-1)*3+1):(3*i-1),:),newtestShape{segment}((3*i),:))
%             hold on
%         end
        newrotorcurve(((starter-1)*3+1):(ender*3),:)=newtestShape{segment};
    end
    newrotor{CURVE}=newrotorcurve;
    newcurveCP{CURVE}=P3;
    newcurveW{CURVE}=W3;
    newcurveXi{CURVE}=Xi3;
end
rotorbreakpoints{1}=[5,6,13,14,15];
rotorbreakpoints{2}=[7,10,14,21,41,60];
rotorbreakpoints{3}=[53];
bkpts= rotorbreakpoints;
clear CP W Xi
CP=newcurveCP; W=newcurveW; Xi=newcurveXi;
    