% This file's purpose is to compare quadtree, gauss-green, and szego-green...
% quadrature methods using a small suite of test cases defined by two
% domains and various integrands:
% Domains:
% close all;
% clear all;
shape=3;
testIntegrands=0;

timings=0;
strgctr=1;
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

% figNum ==1 yields a square plate with a circular hole
% figNum ==2 yields an L-bracket with 3 holes
% figNum ==3 yields a wrench figure
% figNum ==4 yields a guitar-shaped object
% figNum ==5 yields a treble clef
filenames={'circle'};
f{1}=figure('PaperPositionMode','auto','Units','inches','Position',[-0.010416666666667  -0.010416666666667   5.5   3.25],'Color','w');
f{2}=figure('PaperPositionMode','auto','Units','inches','Position',[-0.010416666666667  -0.010416666666667   5.5   3.25],'Color','w');
for jjj=1:1
    for mondegs=1:6%length(integrands)
        functs=monfuncts(((mondegs-1)*(mondegs)/2 +1):((mondegs)*(mondegs+1)/2));
        kg=mondegs-1;
        orientations=[1 -1 -1 -1 -1];
        [avgDAT,DAT]=convergenceAnalysis(filenames{jjj},...
                                            orientations,...
                                            functs,...
                                            kg,...
                                            timings);
            figure(f{1})
            subplot(4,2,mondegs+2)
            loglog(avgDAT{1}(1,:),abs(avgDAT{1}(2,:)),'mo-','MarkerSize',3);
            hold on
            loglog(avgDAT{2}(1,:),abs(avgDAT{2}(2,:)),'co-','MarkerSize',3);
            loglog(avgDAT{3}(1,:),abs(avgDAT{3}(2,:)),'ko-','MarkerSize',3);
            loglog(avgDAT{4}(1,:),abs(avgDAT{4}(2,:)),'bx-.','MarkerSize',3);
            loglog(avgDAT{5}(1,:),abs(avgDAT{5}(2,:)),'gx-.','MarkerSize',3);
            loglog(avgDAT{6}(1,:),abs(avgDAT{6}(2,:)),'rx-.','MarkerSize',3);
            loglog(avgDAT{7}(1,:),abs(avgDAT{7}(2,:)),'x-.','Color',[.5 .5 0],'MarkerSize',3);
            %             eqn=func2str(funct); eqn=replace(eqn,'.',''); eqn=replace(eqn,'*','');
%             eqn=eqn(7:end);
            xlabel('Number of Quadrature Points')
            pl=get(gcf);
            pl.FontName='times';
            pl.FontSize=12;
            set(gca,'yscale','log')
            set(gca,'xscale','log')
            yyaxis('right')
            set(gca,'yscale','log')
            xlim([1e1 1e5])
            xticks([ 1e2 1e4 1e6])
            ylim([1e-17 1])
            yticks([1e-15 1e-10 1e-5 1])
            set(gca,'YColor','black','FontName','times');
            if mondegs==4
                ylabel('Integration Error','FontSize',12,'FontName','times','Interpreter','Latex')
            end
            if mondegs==6
                xlabel('# of Quad Points','FontSize',12,'FontName','times','Interpreter','Latex')
            end
            yyaxis('left')
            ylim([1e-17 1])
            set(gca,'YTickLabel',[]);
            if jjj==1
                title(sprintf('%dth degree',mondegs-1),'Interpreter','Latex','FontSize',12,'FontName','times');
            end
            
%             
%             figure(f{2})
%             subplot(4,2,mondegs+2)
%             loglog(avgDAT{1}(3,:),abs(avgDAT{1}(2,:)),'mo-','MarkerSize',3);
%             hold on
%             loglog(avgDAT{2}(3,:),abs(avgDAT{2}(2,:)),'co-','MarkerSize',3);
%             loglog(avgDAT{3}(3,:),abs(avgDAT{3}(2,:)),'ko-','MarkerSize',3);
%             loglog(avgDAT{4}(3,:),abs(avgDAT{4}(2,:)),'bx-.','MarkerSize',3);
%             loglog(avgDAT{5}(3,:),abs(avgDAT{5}(2,:)),'gx-.','MarkerSize',3);
%             loglog(avgDAT{6}(3,:),abs(avgDAT{6}(2,:)),'rx-.','MarkerSize',3);
% %             eqn=func2str(funct); eqn=replace(eqn,'.',''); eqn=replace(eqn,'*','');
% %             eqn=eqn(7:end);
% %             title({'Convergence Analysis of various quadrature schemes',['for the  ', replace(filename,'_',' '),sprintf(' example and $f(x,y) = %s$',eqn)]},'Interpreter','Latex')
% %             xlabel('Number of Quadrature Points')
% %             legend({'Exact mesh + 3rd order Gauss',...
% %                 'Linear mesh + 1st order Gauss',...
% %                 'Quadtree + 3rd order Gauss',...
% %                 'Cubic Spline appr. + Gauss-Green',...
% %                 'Parametric Gauss-Green',...
% %                 'Exact Rational-Green'},'Location','southeast')
%             pl=get(gcf);
%             pl.FontName='times';
%             pl.FontSize=12;
%             set(gca,'yscale','log')
%             set(gca,'xscale','log')
%             yyaxis('right')
%             set(gca,'yscale','log')
%             xlim([1e-3 10])
%             xticks([10^-3 10^-2 10^-1 1 10 100]);
%             ylim([1e-17 1])
%             yticks([1e-15 1e-10 1e-5 1])
%             set(gca,'YColor','black','FontName','times');
%             if mondegs==4
%                 ylabel('Integration Error','FontSize',12,'FontName','times','Interpreter','Latex')
%             end
%             if mondegs==6
%                 xlabel('Time (s)','FontSize',12,'FontName','times','Interpreter','Latex')
%             end
%             yyaxis('left')
%             ylim([1e-17 1])
%             set(gca,'YTickLabel',[]);
%             if jjj==1
%                 title(sprintf('%dth degree',mondegs-1),'Interpreter','Latex','FontSize',12,'FontName','times');
%             end
%             strg{strgctr}=avgDAT;
%             strgctr=strgctr+1;
    end
end

for i=1:2
    figure(f{i})
    ax = get(gcf,'children');
    ind = find(isgraphics(ax,'Legend'));
    legend('DD-Rational mesh','DD-Linear mesh',...
    'DD-Quadtree','GT-Cubic spline',...
    "GT-SPECTRAL",'GT-SPECTRAL PE','GT-Linear',...
    'FontName','times','Interpreter','Latex',...
    'FontSize',12,'NumColumns',2)
    set(gcf,'children',ax([ind:end,1:ind-1]))
    set(legend,'FontSize',12)
end
figure(f{i})
subplot(4,2,2)
hh1=loglog([.5*10^3 .5*10^5],[10^-1 10^-3],'k');
hold on
hh1t=text(.5*10^5,7*10^-3,'$\mathcal{O}(N^{-1})$','Interpreter','Latex','HorizontalAlignment','Left')
hh3=loglog([.5*10^3 .5*10^5],[10^-2 10^-8],'k');
hh3t=text(.5*10^5,5*10^-8,'$\mathcal{O}(N^{-3})$','Interpreter','Latex','HorizontalAlignment','Left')
hh6=loglog([.5*10^3 .5*10^5],[10^-3 10^-15],'k');
hh6t=text(.5*10^5,5*10^-15,'$\mathcal{O}(N^{-6})$','Interpreter','Latex','HorizontalAlignment','Left')
xlim([1e1 1e5])
ylim([1e-17 1])
axis off

for ii=3:8
    hh=subplot(4,2,ii);
    pos=get(hh,'Position');
    pos2=[pos(1)-.05 pos(2) pos(3) pos(4)];
    set(hh,'Position',pos2);
end

set(gcf, 'PaperUnits', 'normalized')
set(gcf, 'PaperPosition', [0 0 1 1])
print -dpdf paperFig_convergence_analysis_circle_comparisons
