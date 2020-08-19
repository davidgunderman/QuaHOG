% This file's purpose is to compare quadtree, gauss-green, and szego-green...
% quadrature methods using a small suite of test cases defined by two
% domains and various integrands:
% Domains:
% close all;
% clear all;

% Load all needed files.
mydir  = which(mfilename);
idcs   = strfind(mydir,filesep);
folder = mydir(1:idcs(end-1)-1);
addpath(genpath(folder));


shape=3;
testIntegrands=0;
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
d=2;
gaussOrders=[2:25];

% figNum ==1 yields a square plate with a circular hole
% figNum ==2 yields an L-bracket with 3 holes
% figNum ==3 yields a wrench figure
% figNum ==4 yields a guitar-shaped object
% figNum ==5 yields a treble clef
timings=1;
filenames={'plate_with_hole','l_bracket','guitar','treble_clef','two_rotors'};
% f{1}=figure('PaperPositionMode','auto','Units','inches','Position',[-0.010416666666667  -0.010416666666667   8.250000000000000   6.770833333333333],'Color','w');
% f{2}=figure('PaperPositionMode','auto','Units','inches','Position',[-0.010416666666667  -0.010416666666667   8.250000000000000   6.770833333333333],'Color','w');
for jjj=5:5
    for iifuncts=3:3%length(integrands)
        functs=polyfuncts(iifuncts);
        kg=2*iifuncts;
        orientations=[1 -1 -1 -1 -1];
        [avgDAT,DAT]=convergenceAnalysis(filenames{jjj},...
                                            orientations,...
                                            functs,...
                                            kg,...
                                            timings);
%         figure(f{1})
%         subplot(6,4,(jjj)*4+iifuncts+1)
%         loglog(avgDAT{1}(1,:),abs(avgDAT{1}(2,:)),'mo-','MarkerSize',3);
%         hold on
%         loglog(avgDAT{2}(1,:),abs(avgDAT{2}(2,:)),'co-','MarkerSize',3);
%         loglog(avgDAT{3}(1,:),abs(avgDAT{3}(2,:)),'ko-','MarkerSize',3);
%         loglog(avgDAT{4}(1,:),abs(avgDAT{4}(2,:)),'bx-.','MarkerSize',3);
%         loglog(avgDAT{5}(1,:),abs(avgDAT{5}(2,:)),'gx-.','MarkerSize',3);
%         loglog(avgDAT{6}(1,:),abs(avgDAT{6}(2,:)),'rx-.','MarkerSize',3);
        loglog(avgDAT{7}(1,:),abs(avgDAT{7}(2,:)),'x-.','Color',[.5 .5 0],'MarkerSize',3,'Parent',hh(iidx(5*(iifuncts-1)+jjj)));
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
%         pl=get(gcf);
%         pl.FontName='times';
%         pl.FontSize=12;
%         set(gca,'yscale','log')
%         set(gca,'xscale','log')
%         yyaxis('right')
%         set(gca,'yscale','log')
%         xlim([1e1 1e7])
%         xticks([ 1e2 1e4 1e6])
%         ylim([1e-17 1])
%         yticks([1e-15 1e-10 1e-5 1])
%         set(gca,'YColor','black','FontName','times');
%         if jjj==3 && iifuncts==3
%             ylabel('Integration Error','FontSize',12,'Interpreter','Latex')
%         end
%         if jjj==5 && iifuncts==2
%             xlabel('# of Quad points','FontSize',12,'Interpreter','Latex')
%         end
%         yyaxis('left')
%         ylim([1e-17 1])
%         set(gca,'YTickLabel',[]);
%         if jjj==1
%             title(sprintf('$p_%d(x,y)$',iifuncts),'Interpreter','Latex','FontSize',12);
%         end
%         
%         
%         
%         figure(f{2})
%         subplot(6,4,(jjj)*4+iifuncts+1)
%                     
%             loglog(avgDAT{1}(4,:),abs(avgDAT{1}(2,:)),'mo-','MarkerSize',3);
%             hold on
%             loglog(avgDAT{2}(4,:),abs(avgDAT{2}(2,:)),'co-','MarkerSize',3);
%             loglog(avgDAT{3}(4,:),abs(avgDAT{3}(2,:)),'ko-','MarkerSize',3);
%             loglog(avgDAT{4}(4,:),abs(avgDAT{4}(2,:)),'bx-.','MarkerSize',3);
%             loglog(avgDAT{5}(4,:),abs(avgDAT{5}(2,:)),'gx-.','MarkerSize',3);
%             loglog(avgDAT{6}(4,:),abs(avgDAT{6}(2,:)),'rx-.','MarkerSize',3);
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
%         pl.FontName='times';
%         pl.FontSize=12;
%         set(gca,'yscale','log')
%         set(gca,'xscale','log')
%         yyaxis('right')
%         set(gca,'yscale','log')
%         xticks([10^-3 10^-2 10^-1 1 10]);
%         xlim([10^-3 10^2])
%         ylim([1e-17 1])
%         yticks([1e-15 1e-10 1e-5 1])
%         set(gca,'YColor','black','FontName','times');
%         if jjj==3 && iifuncts==3
%             ylabel('Integration Error','FontSize',14,'Interpreter','Latex')
%         end
%         if jjj==5 && iifuncts==2
%             xlabel('Timing (s)','FontSize',14,'Interpreter','Latex')
%         end
%         yyaxis('left')
%         ylim([1e-17 1])
%         set(gca,'YTickLabel',[]);
%         if jjj==1
%             title(sprintf('$p_%d(x,y)$',iifuncts),'Interpreter','Latex','FontSize',12);
%         end
%         strg{strgctr}=avgDAT;
%         strgctr=strgctr+1;
        end
end
for iii=1:2
    figure(iii)
for jjj=1:5
    subplot(6,4,(jjj)*4 +1)
    delete(gca)
    subplot(6,4,(jjj)*4 +1)
    
    ss=generateTestFigures(jjj);
    vals = plot_rat_bern_poly(ss,2,.01,{},[1 0 0 0])
    ll=get(gca,'Children')
%     for i=1:length(ll)
%         set(ll(i),'SizeData',5)
%     end
    if jjj==5 || jjj==2
%         axis equal
        if jjj==2
            ylim([-2.1,2.1]);
            xlim([-4,4]);
        end
        if jjj==5
            ylim([0,2.25]);
            xlim([-.1,4.2]);
        end
    elseif jjj==4
        axis ij
        xlim([-2 3]);
        ylim([1.1,4.3]);
    else
        axis ij
        axis square
    end
    if jjj==1
        title("Region",'Interpreter','Latex','FontName','times','FontSize',12)
        xlim([-.1,2.1])
        ylim([-.1,2.1])
    end
    msize = @(bla) size(bla,1);
    xlabel(sprintf("$n_c=%d,m=%d$",sum(cellfun(msize,ss))/3,size(ss{1},2)-1),'Interpreter','Latex','FontSize',12)
    axis on
    h=get(gca,'xlabel');
    set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
    set(h, 'Color', 'k')
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
end
subplot(6,4,2+4)
legend('DD-Rational mesh','DD-Linear mesh',...
    'DD-Quadtree','GT-Cubic spline',...
    "GT-SPECTRAL",'GT-SPECTRAL PE',...
    'FontName','times','Interpreter','Latex','NumColumns',2,...
    'FontSize',12)
ax = get(gcf,'children');
ind = find(isgraphics(ax,'Legend'));
set(gcf,'children',ax([ind:end,1:ind-1]))
set(legend,'FontSize',10)
end

tt=subplot(6,4,4)
hh1=loglog([.5*10^3 .5*10^5],[10^-1 10^-3],'k');
hold on
hh1t=text(.5*10^5,7*10^-3,'$\mathcal{O}(n_q^{-1})$','Interpreter','Latex','HorizontalAlignment','Left')
hh3=loglog([.5*10^3 .5*10^5],[10^-2 10^-8],'k');
hh3t=text(.5*10^5,5*10^-8,'$\mathcal{O}(n_q^{-3})$','Interpreter','Latex','HorizontalAlignment','Left')
hh6=loglog([.5*10^3 .5*10^5],[10^-3 10^-15],'k');
hh6t=text(.5*10^5,5*10^-15,'$\mathcal{O}(n_q^{-6})$','Interpreter','Latex','HorizontalAlignment','Left')
xlim([1e1 1e7])
ylim([1e-17 1])
axis off
% 

for ii=5:24
    hh=subplot(6,4,ii);
    pos=get(gca,'Position');
    pos2=[pos(1) pos(2) pos(3)*1.2 pos(4)*1.2];
    set(gca,'Position',pos2);
end
% set(gcf, 'PaperUnits', 'normalized')
% set(gcf, 'PaperPosition', [0 0 1 1])
% print -dpdf paperFig_timing_results
