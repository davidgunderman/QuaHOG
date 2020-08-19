function [int2val,numEvals] = integral2_bools(CPWs,funct,orientations)
global evalCounter;
jj=1;
timing =0;
bds=zeros(length(CPWs),4);
for ii=1:length(CPWs)
    bds(ii,:)=[min(CPWs{ii}(1:3:end,:),[],'all'),...
                max(CPWs{ii}(1:3:end,:),[],'all'),...
                min(CPWs{ii}(2:3:end,:),[],'all'),...
                max(CPWs{ii}(2:3:end,:),[],'all')]...
                +1e-14*[-rand(1), rand(1), -rand(1), rand(1)];
end
while jj<4
    evalCounter=0;
    tic;
    for ii=1:length(CPWs)
        Ifun = @(x,y) shapeIEval([x y],funct,CPWs{ii},1);
        arrayIfun = @(x,y) arrayfun(Ifun,x,y);
        int2valjj(ii)=integral2(arrayIfun,bds(ii,1),bds(ii,2),bds(ii,3),bds(ii,4),'RelTol',2^(-jj));
    end
    int2val(jj)=sum(orientations.*int2valjj);
    numEvals(jj)=evalCounter;
    timing=toc;
    jj=jj+1;
end
    