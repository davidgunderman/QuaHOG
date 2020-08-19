% This file's purpose is to compare quadtree, gauss-green, and szego-green...
% quadrature methods using a small suite of test cases defined by two
% domains and various integrands:
% Domains:
close all;
clear all;
shape=3;
testIntegrands=0;

counter=1;
for i=0:5
    for j=0:i
        a=i-j; b=j;
        monfuncts(counter) = {@(x,y) x.^a.*y.^(b)};
        bla(counter,:)=[a,b];
        counter=counter+1;
    end
end
% 2. Three polynomials of degree 2 (bilinear), 4 (biquadratic), and 6
% (bicubic)
polyfuncts={@(x,y) (2*x.^2 +x.*y - y +2);
            @(x,y) (2*x.^2.*y.^2 +.3*x.^2.*y - y.^4 + 3*x +2);
            @(x,y) (x.^5 - 5*y.^3.*x.^3 + .2*x.^2 + 2*y.*x.^2 +3);};
% 3. A rational function of degree 4 and an exponential function.
otherfuncts={@(x,y) (y.^3 - (x.^3.*y.^2) - (x.*y) -3)./((x.^2).*(y.^2) +10);
             @(x,y) 10*(exp( - x.^2 ) + 2*y);
             @(x,y) sqrt((x+10).^2+(x+10).*(y+10).^2 +x)};
addpath("../Rational_Quadrature/Matlab/Src",...
"../Rational_Quadrature/Matlab/Tests",...
"../Rational_Quadrature/Matlab/ThirdPartySupportingCode")
d=2;
gaussOrders=[2:25];

% figNum ==1 yields a square plate with a circular hole
% figNum ==2 yields an L-bracket with 3 holes
% figNum ==3 yields a wrench figure
% figNum ==4 yields a guitar-shaped object
% figNum ==5 yields a treble clef
filenames={'plate_with_hole','l_bracket','guitar','treble_clef','two_rotors'};
f{1}=figure('Position',[0 0 900 700],'Color','w');
f{2}=figure('Position',[0 0 900 700],'Color','w');
for jjj=1:5
    for mondegs=1:3%length(integrands)
        functs=monfuncts(((mondegs-1)*(mondegs)/2 +1):((mondegs)*(mondegs+1)/2));
        kg=invTri(mondegs-1)+1;
        orientations=[1 -1 -1 -1 -1];
        [avgDAT,DAT]=convergenceAnalysis(filenames{jjj},...
                                            orientations,...
                                            functs,...
                                            kg);
        figure(f{1})
        subplot(5,4,(jjj-1)*4+mondegs+1)
        loglog(avgDAT{1}(1,:),abs(avgDAT{1}(2,:)),'mo-','MarkerSize',5);
        hold on
        loglog(avgDAT{2}(1,:),abs(avgDAT{2}(2,:)),'co-','MarkerSize',5);
        loglog(avgDAT{3}(1,:),abs(avgDAT{3}(2,:)),'ko-','MarkerSize',5);
        loglog(avgDAT{4}(1,:),abs(avgDAT{4}(2,:)),'bx-.','MarkerSize',5);
        loglog(avgDAT{5}(1,:),abs(avgDAT{5}(2,:)),'gx-.','MarkerSize',5);
        loglog(avgDAT{6}(1,:),abs(avgDAT{6}(2,:)),'rx-.','MarkerSize',5);
%             eqn=func2str(funct); eqn=replace(eqn,'.',''); eqn=replace(eqn,'*','');
%             eqn=eqn(7:end);
%             title({'Convergence Analysis of various quadrature schemes',['for the  ', replace(filename,'_',' '),sprintf(' example and $f(x,y) = %s$',eqn)]},'Interpreter','Latex')
%             xlabel('Number of Quadrature Points')
%             legend({'Exact mesh + 3rd order Gauss',...
%                 'Linear mesh + 1st order Gauss',...
%                 'Quadtree + 3rd order Gauss',...
%                 'Cubic Spline appr. + Gauss-Green',...
%                 'Parametric Gauss-Green',...
%                 'Exact Rational-Green'},'Location','southeast')
        pl=get(gcf);
        pl.FontName='times';
        pl.FontSize=12;
        set(gca,'yscale','log')
        set(gca,'xscale','log')
        yyaxis('right')
        set(gca,'yscale','log')
        xlim([1e1 1e7])
        xticks([ 1e2 1e4 1e6])
        ylim([1e-17 1])
        yticks([1e-15 1e-10 1e-5 1])
        set(gca,'YColor','black','FontName','times');
        if jjj==3 && mondegs==3
            ylabel('Integration Error','FontSize',14,'Interpreter','Latex')
        end
        if jjj==5 && mondegs==2
            xlabel('# of Quad points','FontSize',14)
        end
        yyaxis('left')
        ylim([1e-17 1])
        set(gca,'YTickLabel',[]);
        if jjj==1
            title(sprintf('%dth order',mondegs-1),'Interpreter','Latex','FontSize',12);
        end
        
        
        
        figure(f{2})
        subplot(5,4,(jjj-1)*4+mondegs+1)
                    
            loglog(avgDAT{1}(3,:),abs(avgDAT{1}(2,:)),'mo-','MarkerSize',5);
            hold on
            loglog(avgDAT{2}(3,:),abs(avgDAT{2}(2,:)),'co-','MarkerSize',5);
            loglog(avgDAT{3}(3,:),abs(avgDAT{3}(2,:)),'ko-','MarkerSize',5);
            loglog(avgDAT{4}(3,:),abs(avgDAT{4}(2,:)),'bx-.','MarkerSize',5);
            loglog(avgDAT{5}(3,:),abs(avgDAT{5}(2,:)),'gx-.','MarkerSize',5);
            loglog(avgDAT{6}(3,:),abs(avgDAT{6}(2,:)),'rx-.','MarkerSize',5);
%             eqn=func2str(funct); eqn=replace(eqn,'.',''); eqn=replace(eqn,'*','');
%             eqn=eqn(7:end);
%             title({'Convergence Analysis of various quadrature schemes',['for the  ', replace(filename,'_',' '),sprintf(' example and $f(x,y) = %s$',eqn)]},'Interpreter','Latex')
%             xlabel('Number of Quadrature Points')
%             legend({'Exact mesh + 3rd order Gauss',...
%                 'Linear mesh + 1st order Gauss',...
%                 'Quadtree + 3rd order Gauss',...
%                 'Cubic Spline appr. + Gauss-Green',...
%                 'Parametric Gauss-Green',...
%                 'Exact Rational-Green'},'Location','southeast')
            pl=get(gcf);
        pl.FontName='times';
        pl.FontSize=12;
        set(gca,'yscale','log')
        set(gca,'xscale','log')
        yyaxis('right')
        set(gca,'yscale','log')
        xticks([10^-2 10^-1 1 10 100]);
        xlim([10^-3 10^2])
        ylim([1e-17 1])
        yticks([1e-15 1e-10 1e-5 1])
        set(gca,'YColor','black','FontName','times');
        if jjj==3 && mondegs==3
            ylabel('Integration Error','FontSize',14,'Interpreter','Latex')
        end
        if jjj==5 && mondegs==2
            xlabel('Timing (s)','FontSize',14,'Interpreter','Latex')
        end
        yyaxis('left')
        ylim([1e-17 1])
        set(gca,'YTickLabel',[]);
        if jjj==1
            title(sprintf('%dth order',mondegs-1),'Interpreter','Latex','FontSize',12);
        end
    end
end
for iii=1:2
    figure(iii)
for jjj=1:5
    subplot(5,4,(jjj-1)*4 +1)
    delete(gca)
    subplot(5,4,(jjj-1)*4 +1)
    
    ss=generateTestFigures(jjj);
    vals = plot_rat_bern_poly(ss,2,.01,{},[1 0 0 0])
    ll=get(gca,'Children')
%     for i=1:length(ll)
%         set(ll(i),'SizeData',5)
%     end
    if jjj==5 || jjj==2
        axis equal
        if jjj==2
            ylim([-2.1,2.1]);
        end
        if jjj==5
            ylim([0,2.25]);
            xlim([-.1,4.2]);
        end
    elseif jjj==4
        axis ij
        xlim([-1 2]);
        ylim([1.1,4.3]);
    else
        axis ij
        axis square
    end
    if jjj==1
        title("Region",'Interpreter','Latex','FontSize',12)
        xlim([-.1,2.1])
        ylim([-.1,2.1])
    end
    msize = @(bla) size(bla,1);
    xlabel(sprintf("$M=%d,n_c=%d,p=%d$",length(ss),sum(cellfun(msize,ss))/3,size(ss{1},2)-1),'Interpreter','Latex','FontSize',10)
    axis on
    h=get(gca,'xlabel');
    set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
    set(h, 'Color', 'k')
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
end
% subplot(5,4,2+4)
% legend('Exact mesh','Linear mesh',...
%     'Quadtree','Cubic spline',...
%     'Present method','Mod. method','Interpreter','Latex')
% ax = get(gcf,'children');
% ind = find(isgraphics(ax,'Legend'));
% set(gcf,'children',ax([ind:end,1:ind-1]))
end
% 
set(gcf, 'PaperUnits', 'normalized')
set(gcf, 'PaperPosition', [0 0 1 1])
print -dpdf paperFig_convergence_analysis_monomial_comparisons
