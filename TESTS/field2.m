function val = field2(x,y,funct)
global evalCounter;
evalCounter=evalCounter+length(x(:));
val=funct(x,y);
global scatterEvals
if scatterEvals==1
    if size(x,2)>1
        figure(7)
    else
        figure(1)
    end
    scatter(x(:),y(:),50,[0 .7 0],'s','filled');
    hold on
end
% scatter(x,y,200,[0 .7 0],'.');
% scatter(x,y,200,[0 0 .7],'.');
% scatter(x,y,200,[.7 .35 0],'.');
% hold on
end