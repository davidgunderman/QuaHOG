% This file's purpose is to compare our gauss-green, and szego-green...
% with various other quadrature methods using a small suite of test cases 
% defined by 5 domains and various integrands:
% Domains:
close all;
clear all;

% 2. Three polynomials of degree 2 (bilinear), 4 (biquadratic), and 6
% (bicubic)
polyfuncts={@(x,y) (2*x.^2 +x.*y - y +2);
            @(x,y) (2*x.^2.*y.^2 +.3*x.^2.*y - y.^4 + 3*x +2);
            @(x,y) (x.^5 - 5*y.^3.*x.^3 + .2*x.^2 + 2*y.*x.^2 +3);};
% 3. A rational function of degree 4 and an exponential function.
otherfuncts={@(x,y) (y.^3 - (x.^3.*y.^2) - (x.*y) -3)./((x.^2).*(y.^2) +30);
             @(x,y) (exp( - x.^2 ) + 2*y);
             @(x,y) sqrt((x+10).^2+(x+10).*(y+10).^2 +x)};
addpath("../Rational_Quadrature/Matlab/Src",...
"../Rational_Quadrature/Matlab/Tests",...
"../Rational_Quadrature/Matlab/ThirdPartySupportingCode")
d=2;

timedfuncts={polyfuncts{3}; otherfuncts{3}};
% figNum ==1 yields a square plate with a circular hole
% figNum ==2 yields an L-bracket with 3 holes
% figNum ==3 yields a wrench figure
% figNum ==4 yields a guitar-shaped object
% figNum ==5 yields a treble clef
filenames={'plate_with_hole','l_bracket','guitar','treble_clef','two_rotors'};
f{1}=figure('Position',[0 0 400 600],'Color','w');
    for iii=1:2
        functs=timedfuncts(iii);
        if iii==1
            kg=6;
        else
            kg=2.5;
        end
        orientations=[1 -1 -1 -1 -1];
        [avgDAT,DAT]=convergenceAnalysis(filenames{5},...
                                            orientations,...
                                            functs,...
                                            kg);
            figure(f{1})
            subplot(3,2,iii+2)
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
            xticks([10^-1 1 10 100 1000]);
            xlim([10^-1 10^3])
            ylim([1e-17 1])
            yticks([1e-15 1e-10 1e-5 1])
            set(gca,'YColor','black','FontName','times');
            xlabel('Pre-processing (s)','FontSize',10)
            if iii==2
                ylabel('Integration Error','FontSize',12)
            end
            yyaxis('left')
            ylim([1e-17 1])
            set(gca,'YTickLabel',[]);
            if iii==1
                title("$p_3(x,y)$",'Interpreter','Latex','FontSize',12);
            end
            if iii==2
                title("$f_3(x,y)$",'Interpreter','Latex','FontSize',12);
            end

            subplot(3,2,iii+4)

                loglog(avgDAT{1}(4,:),abs(avgDAT{1}(2,:)),'mo-','MarkerSize',5);
                hold on
                loglog(avgDAT{2}(4,:),abs(avgDAT{2}(2,:)),'co-','MarkerSize',5);
                loglog(avgDAT{3}(4,:),abs(avgDAT{3}(2,:)),'ko-','MarkerSize',5);
                loglog(avgDAT{4}(4,:),abs(avgDAT{4}(2,:)),'bx-.','MarkerSize',5);
                loglog(avgDAT{5}(4,:),abs(avgDAT{5}(2,:)),'gx-.','MarkerSize',5);
                loglog(avgDAT{6}(4,:),abs(avgDAT{6}(2,:)),'rx-.','MarkerSize',5);
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
            xticks([10^-5 10^-4 10^-3 10^-2 10^-1 1]);
            xlim([10^-5 10])
            ylim([1e-17 1])
            yticks([1e-15 1e-10 1e-5 1])
            set(gca,'YColor','black','FontName','times');
            xlabel('Evaluation (s)','FontSize',10)
            yyaxis('left')
            ylim([1e-17 1])
            set(gca,'YTickLabel',[]);
    end

subplot(3,2,4)
legend('DD-Exact mesh','DD-Linear mesh',...
    'DD-Quadtree','GT-Cubic spline',...
    "GT-SPECTRAL",'GT-SPECTRAL PE',...
    'FontName','times','Interpreter','Latex','NumColumns',2)
ax = get(gcf,'children');
ind = find(isgraphics(ax,'Legend'));
set(gcf,'children',ax([ind:end,1:ind-1]))
set(legend,'FontSize',10)

for ii=3:6
    hh=subplot(3,2,ii);
    pos=get(hh,'Position');
    pos2=[pos(1)-.05 pos(2) pos(3) pos(4)];
    set(hh,'Position',pos2);
end