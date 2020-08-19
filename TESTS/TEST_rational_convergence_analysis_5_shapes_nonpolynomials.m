% Produces figure 9 in "spectral mesh-free quadrature for planar regions
% bounded by rational parametric curves.
% Tests the SPECTRAL and SPECTRALPE methods' convergence curves against
% four other methods.
% close all;
% clear all;

% Load all needed files.
mydir  = which(mfilename);
idcs   = strfind(mydir,filesep);
folder = mydir(1:idcs(end-1)-1);
addpath(genpath(folder));

strgctr=1;
counter=1;

% Definitions of three test functions:
%   a rational function of degree 4
%   an exponential function
%   a square roots function
functs={@(x,y) (y.^3 - (x.^3.*y.^2) - (x.*y) -3)./((x.^2).*(y.^2) +30);
             @(x,y) (exp( - x.^2 ) + 2*y);
             @(x,y) sqrt((x+10).^2+(x+10).*(y+10).^2 +x)};
         
% Test on five functions defined by NURBS:
filenames={'plate_with_hole','l_bracket','guitar','treble_clef','two_rotors'};

% % Create figures for plotting
% f{1}=figure('PaperPositionMode','auto',...
%     'Units','inches',...
%     'Position',[-0.010416666666667  -0.010416666666667   8.250000000000000   6.770833333333333],...
%     'Color','w');

% This sets the "k" parameter for each function, which determines the
% degree of polynomial exactness used in the SPECTRALPE for each function.
% Currently, they are all set to "2.5", this ensures exactness for
% polynomials up to degree 2.
kgpts=[2.5 2.5 2.5];

% Perform the analysis for each of the 5 shapes and each of the 3
% integrands:
for jjj=1:5
    for iii=1:3
        funct=functs(iii);
        
        % All of the figures have one counter-clockwise oriented loops,
        % then the rest are clock-wise oriented. This is only used by the
        % comparison methods, not by the SPECTRAL or SPECTRALPE methods.
        orientations=[1 -1 -1 -1 -1];
        
        % Perform convergence analysis (Note: this could take a long time,
        % since the comparison methods can take a long time. To turn off
        % the comparison methods, edit the code inside "convergence
        % analysys")
        [avgDAT,DAT]=convergenceAnalysis(filenames{jjj},...
                                            orientations,...
                                            funct,...
                                            kgpts(iii),...
                                            0);
                                        
            % Plot data on log log plots
%             figure(f{1})
%             subplot(6,4,(jjj)*4+1+iii)
%             loglog(avgDAT{1}(1,:),abs(avgDAT{1}(2,:)),'mo-','MarkerSize',3);
%             hold on
%             loglog(avgDAT{2}(1,:),abs(avgDAT{2}(2,:)),'co-','MarkerSize',3);
%             loglog(avgDAT{3}(1,:),abs(avgDAT{3}(2,:)),'ko-','MarkerSize',3);
%             loglog(avgDAT{4}(1,:),abs(avgDAT{4}(2,:)),'bx-.','MarkerSize',3);
%             loglog(avgDAT{5}(1,:),abs(avgDAT{5}(2,:)),'gx-.','MarkerSize',3);
%             loglog(avgDAT{6}(1,:),abs(avgDAT{6}(2,:)),'rx-.','MarkerSize',3);
            loglog(avgDAT{7}(1,:),abs(avgDAT{7}(2,:)),'x-.','Color',[.5 .5 0],'MarkerSize',3,'Parent',hh(iidx(3*(jjj-1)+iii)));
            % Format plot
%             pl=get(gcf);
%             pl.FontName='times';
%             pl.FontSize=12;
%             set(gca,'yscale','log')
%             set(gca,'xscale','log')
%             yyaxis('right')
%             set(gca,'yscale','log')
%             xlim([1e1 1e7])
%             xticks([ 1e2 1e4 1e6])
%             ylim([1e-17 1])
%             yticks([1e-15 1e-10 1e-5 1])
%             set(gca,'YColor','black','FontName','times');
%             if jjj==3 && iii==3
%                 ylabel('Integration Error','FontSize',14,'Interpreter','Latex')
%             end
%             if jjj==5 && iii==2
%                 xlabel('# of Quad points','FontSize',14,'Interpreter','Latex')
%             end
%             yyaxis('left')
%             ylim([1e-17 1])
%             set(gca,'YTickLabel',[]);
%             if jjj==1
%                 title(sprintf("$f_%d(x,y)$",iii),'Interpreter','Latex','FontSize',12);
%             end
    end
end

% Plot regions
for jjj=1:5
    subplot(6,4,(jjj)*4 +1)
    delete(gca)
    subplot(6,4,(jjj)*4 +1)

    ss=generateTestFigures(jjj);
    vals = plot_rat_bern_poly(ss,2,.01,{},[1 0 0 0]);
    ll=get(gca,'Children');
    if jjj==5 || jjj==2
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

% Create legend
subplot(6,4,2+4)
legend('DD-Rational mesh','DD-Linear mesh',...
    'DD-Quadtree','GT-Cubic spline',...
    "GT-SPECTRAL",'GT-SPECTRAL PE',...
    'FontName','times','Interpreter','Latex','NumColumns',2);
ax = get(gcf,'children');
ind = find(isgraphics(ax,'Legend'));
set(gcf,'children',ax([ind:end,1:ind-1]));
set(legend,'FontSize',10);

% Create convergence lines in upper-right hand corner
subplot(6,4,4)
hh1=loglog([.5*10^3 .5*10^5],[10^-1 10^-3],'k');
hold on
hh1t=text(.5*10^5,7*10^-3,'$\mathcal{O}(n_q^{-1})$','Interpreter','Latex','HorizontalAlignment','Left');
hh3=loglog([.5*10^3 .5*10^5],[10^-2 10^-8],'k');
hh3t=text(.5*10^5,5*10^-8,'$\mathcal{O}(n_q^{-3})$','Interpreter','Latex','HorizontalAlignment','Left');
hh6=loglog([.5*10^3 .5*10^5],[10^-3 10^-15],'k');
hh6t=text(.5*10^5,5*10^-15,'$\mathcal{O}(n_q^{-6})$','Interpreter','Latex','HorizontalAlignment','Left');
xlim([1e1 1e7]);
ylim([1e-17 1]);
axis off
% 

% Increase size of all plots
for ii=4:24
    hh=subplot(6,4,ii);
    pos=get(hh,'Position');
    pos2=[pos(1) pos(2) pos(3)*1.2 pos(4)*1.2];
    set(hh,'Position',pos2);
end
