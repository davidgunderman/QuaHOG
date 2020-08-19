function toSplineFile(filename,ControlPoints,Weights,KV)

format long;
filename = [filename,'.spline'];
fileID = fopen(filename,'w+t');

NCURVES=length(ControlPoints);
NSD=length(ControlPoints);

fprintf(fileID,'%s \n','HEADER');
fprintf(fileID,'%s %s \n','NCURVES','NSD');
fprintf(fileID,'%u %u \n',NCURVES,NSD);
for i=1:NCURVES
%     ControlPoints{i}=ControlPoints{i}.*Weights{i};
    NCP=length(ControlPoints{i});
    LKV=length(KV{i});
    p=LKV-NCP-1;
    fprintf(fileID,'%s %s %s \n','NCP','LKV','p');
    fprintf(fileID,'%u %u %u \n',NCP,LKV,p);
    fprintf(fileID,'%s \n','CONTROL POINTS');
    for j=1:(length(Weights{i}))
        fprintf(fileID,'%.16f %.16f %.16f \n',ControlPoints{i}(1,j),ControlPoints{i}(2,j),Weights{i}(j));
    end
    fprintf(fileID,'%s \n','KNOT VECTOR');
    fprintf(fileID,['' repmat('%g  ', 1, numel(KV{i})-1), '%g \n'],KV{i});
    fprintf(fileID,'%s \n','BOUNDARY FLAGS');
    for j=1:(KV{i}(end))
        fprintf(fileID,'%u %u \n',j,1);
    end
end
fprintf(fileID,'%s \n','BOUNDARY CONDITIONS');
fprintf(fileID,'%s \n','NBOUNDS');
fprintf(fileID,'%u \n', 1);
fprintf(fileID,'%u %u \n',1, 6);
fclose(fileID);
end