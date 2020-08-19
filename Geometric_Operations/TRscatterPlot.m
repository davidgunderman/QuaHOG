function f1 = TRscatterPlot(ctrlpoints,dim,res,color)
    ptcntr = 1;
    ord = size(ctrlpoints,1)-1;
    for i=0:res
        for j=0:(res-i)
            k=res-i-j;
            pts2plot(ptcntr,:)=TdCReval(ctrlpoints,[i/res,j/res,k/res],dim);
            ptcntr=ptcntr+1;
        end
    end
    if dim==3
        f1=scatter3(pts2plot(:,1),pts2plot(:,2),pts2plot(:,3),20,-sqrt(sum(pts2plot.^2,2)),color,'filled');
    else
        f1=scatter(pts2plot(:,1),pts2plot(:,2),20,-sqrt(sum(pts2plot.^2,2)),color,'filled');
end