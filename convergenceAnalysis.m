function [avgDAT,DAT] = convergenceAnalysis(filename,orientations,functs,poldeg,timings)
    % Performs convergence analysis of the integral of funct over the
    % interior of the region defined by the set of oriented NURBS cuves
    % given by CP, W, and Xi using six strategies:
    % 1. Quadtree integration
    %    Uses my own adaptive quadtree implementation
    % 2. Linear meshing and order 3 Gaussian quadrature under h-refinement
    %    using a modified version of mesh2d (xmesh provided in TRIGA)
    % 3. Exact rational meshing and order 3 Gaussian quadrature under
    %    p-refinement using xmesh2d from TRIGA
    % 4. Polynomial Spline approximation of the boundary followed by
    %    Green's theorem-based quadrature using splinegauss
    % 5. SPECTRAL quadrature using our method
    % 6. SPECTRALPE quadrature using our method
    %
    %   Input:
    %       filename    A string containing the filename where the boundary 
    %                   NURBS curves are stored, with the extension _NURBS 
    %                   removed. (e.g. two_rotors will load the file
    %                   two_rotors_NURBS). These NURBS curves bound the 
    %                   integraiton region.
    %       orientations A vectors of 1's and -1's representing the
    %                   orientations of each of the boundary loops in the
    %                   region in filename
    %       functs      A list of functions for which integrals will be
    %                   computed.
    %       poldegs     If the functions are polynomial, this contains an
    %                   integer for each functs which describes its degree.
    %                   If the functions are nonpolynomial, this contains a
    %                   non-integer representing what degree of polynomial
    %                   exactness should be guaranteed for the SPECTRALPE
    %                   method.
    %       timings     0 if timings are not desired. 1 if timings are
    %                   Note that timings==1 takes much longer.
    %   Output:
    %       avgDAT      Contains four lists of data for each
    %                   method:
    %                       1. Number of quad points
    %                       2. Average quadrature error for all functions in
    %                       functs for the associated # of quad points
    %                       3. Pre-processing time for the associated # of quad
    %                       points (if timings==1)
    %                       4. Average evaluation time for all functions in
    %                       functs for the associate # of quadrature points
    %                       (if timings ==1)
    %       DAT         Formatted the same as avgDAT, except contains one cell
    %                   for each function (i.e., is not an average over the
    %               functs).
    
    % Turn off irritating warnings
    warning('off','all'); 
    
    % This hard-codes how many iterations of each strategy should be
    % employed.
    maxmeshQuad=10;
    if strcmp(filename,'treble_clef') || strcmp(filename,'guitar')
        maxquadmeshlin=5;
        maxGGquadssq=35;
        numQuadquad=3;
    elseif strcmp(filename,'two_rotors')
        maxquadmeshlin=5;
        maxGGquadssq=55;
        numQuadquad=4;
    else
        maxquadmeshlin=6;
        maxGGquadssq=25;
        numQuadquad=5;
    end
    maxGGquads=ceil(sqrt(maxGGquadssq));
    if mod(poldeg,1)==0
        maxSGquads=1;
        maxVianquad=10;
        maxVianquadLin=10;
    else
        maxSGquads=maxGGquads;
        maxVianquad=9;
        maxVianquadLin=9;
    end
    
    %% On off switch 
    % If only some of the comparison methods are desired, the others can be
    % turned off using this onoff list (0= off, 1= on). The order is the
    % same as in the description of the function.
    onoff=[0 0 0 0 1 1 0];
    
    for jj=1:1
        if ~onoff(1)
            maxmeshQuad=0;
        end
        if ~onoff(2)
            maxquadmeshlin=0;
        end
        if ~onoff(3)
            numQuadquad=0;
        end
        if ~onoff(4)
            maxVianquad=0;
        end
        if ~onoff(5)
            maxGGquads=0;
        end
        if ~onoff(6)
            maxSGquads=0;
        end
        if ~onoff(7)
            maxVianquadLin=0;
        end
    end
    
    %% Initialize Output
    avgDAT=cell(1,7);
    avgDAT{1}=zeros(4,maxmeshQuad);
    avgDAT{2}=zeros(4,maxquadmeshlin);
    avgDAT{3}=zeros(4,numQuadquad);
    avgDAT{4}=zeros(4,maxVianquad);
    avgDAT{5}=zeros(4,maxGGquads);
    avgDAT{6}=zeros(4,maxSGquads);
    avgDAT{7}=zeros(4,maxVianquadLin);
    for funi=1:length(functs)
        DAT{funi}=avgDAT;
    end
    
    % Load the NURBS curves.
    load([filename,'_NURBS'],'CP','W','Xi');
    
    % Convert the NURBS format to my format
    CPW_my_format=extract_my_format_from_NURBS(CP,W,Xi);
    
    % Save the NURBS format as a spline file for TRIGA
    toSplineFile(filename,CP,W,Xi);
    [xtrue,ytrue,wtrue] = SPECTRAL_quads(CPW_my_format,55,55);

           
    for funi=1:length(functs)
        % Creates new function which can track number of evaluations
        funct=functs{funi};

        % Integrate using szego-gauss
        fprintf('Performing our Fejer-Gauss method \n');
        p=length(Xi{1})-length(W{1})-1;
        if mod(poldeg,1)==0
            timingstrue=zeros(1,1);
            timingsAtrue=zeros(1,1);
            f = @() SPECTRALPE_quads(CPW_my_format,0,p*(poldeg+3),ceil((poldeg+1)/2));
            [xq,yq,wq]=f();
            nQuadsSzego=length(xq);
            g= @() applyRule(xq,yq,wq,funct);
            truev= g();
            if timings
                timingsAtrue=timeit(g);
                timingstrue=timeit(f);
            end
        else
            truev=zeros(1,maxSGquads);
            nQuadsSzego=zeros(1,maxSGquads);
            timingstrue=zeros(1,maxSGquads);
            timingsAtrue=zeros(1,maxSGquads);
            for ii=1:maxSGquads
                f = @() SPECTRALPE_quads(CPW_my_format,ii^2,p*(floor(poldeg)+3),ii^2);
                [xq,yq,wq]=f();
                g= @() applyRule(xq,yq,wq,funct);
                truev(ii)= g();
                if timings
                    timingstrue(ii)=timeit(f);
                    timingsAtrue(ii)=timeit(g);
                end
                nQuadsSzego(ii)=length(xq);
            end
        end


        % Integrate using green-gauss
        nQuadsGG=zeros(1,maxGGquads); valueGG=zeros(1,maxGGquads); timingsGG=zeros(1,maxGGquads);
        timingsAGG=zeros(1,maxGGquads);
        for ii=1:maxGGquads
            fprintf('Performing our Gauss-Green method with %d orders of accuracy \n',ii^2);
            tic
            if mod(poldeg,1)==0
                f = @() SPECTRAL_quads(CPW_my_format,ii^2,ceil((poldeg+1)/2));
    %             valueGG(ii)=-SPECTRALquads(CPW_my_format,field,ii^2,ceil((poldeg+1)/2));
            else
                f = @() SPECTRAL_quads(CPW_my_format,ii^2,ii^2);
    %             valueGG(ii)=-SPECTRALquads(CPW_my_format,field,ii^2,ii^2);
            end
            [xq,yq,wq]=f();
            nQuadsGG(ii)=length(xq);
            g= @() applyRule(xq,yq,wq,funct);
            valueGG(ii)= g();
            if timings
                timingsAGG(ii)=timeit(g);
                timingsGG(ii)=timeit(f);
            end
        end
        
        fprintf("Computing true value using SPECTRAL quadrature with 55 quad points \n");
        trueI(funi) = applyRule(xtrue,ytrue,wtrue,funct);
        DAT{funi}([5,6])={...
            [nQuadsGG; valueGG; timingsGG; timingsAGG],...
            [nQuadsSzego; truev; timingstrue; timingsAtrue]};
    end
         
    
    % Integrate using quadtree decomposition using 2nd order
    % Gaussian quadrature
    warning('off','all');
    nQuadsquad = zeros(length(functs),numQuadquad); valuequad = zeros(length(functs),numQuadquad);
    timingsQuadquad = zeros(1,numQuadquad); timingsAQuadquad=zeros(length(functs),numQuadquad);
%     initlvls=0;
    for intlvls=1:numQuadquad
        tol=1e-8;
        fprintf('Performing quadtree quadrature with %d maximum levels \n',intlvls);
        f=@() quadtreeIntegratePts(CPW_my_format,intlvls,3,tol,orientations);
        [xq,yq,wq]=f();
        nQuadsquad(:,intlvls)=length(xq);
        if timings
            timingsQuadquad(intlvls)=timeit(f);
        end
        for funi=1:length(functs)
            funct=functs{funi};
            g=@() applyRule(xq,yq,wq,funct);
            valuequad(funi,intlvls)=g();
            if timings
                timingsAQuadquad(funi,intlvls)=timeit(g);
            end
        end
    end
    for funi=1:length(functs)
        DAT{funi}(3)= {[nQuadsquad(funi,:); valuequad(funi,:); timingsQuadquad; timingsAQuadquad(funi,:);]};
    end
    
    % Integrate using Sommariva/Vianello strategy (approximate boundary by
    % polynomial spline, then effectively use gauss-green)
%     integrateVian(CPWmat,numSplinepts,
    nQuadsVian=zeros(length(functs),maxVianquad); valueVian=zeros(length(functs),maxVianquad);
    timingsVian=zeros(1,maxVianquad); timingsAVian=zeros(length(functs),maxVianquad);
    for ll=3:maxVianquad
        kk=2^(ll-1);
        fprintf('Performing cubic spline approximation quadrature with %d sample points per curve \n',kk);
            if mod(poldeg,1)==0
                N_2009=10;
            else
                N_2009=22;
            end
            f=@() performVianQuadPts(CPW_my_format,kk,N_2009);
            [xq,yq,wq]=f();
            for funi=1:length(functs)
                funct=functs{funi};
                g=@()applyRule(xq,yq,wq,funct);
                valueVian(funi,ll)= g();
                if timings
                    timingsAVian(funi,ll)=timeit(g);
                end
            end
            nQuadsVian(:,ll)=length(xq);
            if timings
                timingsVian(ll)=timeit(f);
            end
    end
    for funi=1:length(functs)
        DAT{funi}(4)= {[nQuadsVian(funi,:); valueVian(funi,:); timingsVian; timingsAVian(funi,:);]};
    end
   
        % Integrate using Sommariva/Vianello strategy (approximate boundary by
    % polynomial spline, then effectively use gauss-green)
%     integrateVian(CPWmat,numSplinepts,
    nQuadsVianLin=zeros(length(functs),maxVianquadLin); valueVianLin=zeros(length(functs),maxVianquadLin);
    timingsVianLin=zeros(1,maxVianquadLin); timingsAVianLin=zeros(length(functs),maxVianquadLin);
    for ll=3:maxVianquadLin
        kk=2^(ll-1);
        fprintf('Performing linear approximation quadrature with %d sample points per curve \n',kk);
            if mod(poldeg,1)==0
                N_2009=5;
            else
                N_2009=15;
            end
            f=@() performVianQuadPtsLin(CPW_my_format,kk,N_2009);
            [xq,yq,wq]=f();
            for funi=1:length(functs)
                funct=functs{funi};
                g=@()applyRule(xq,yq,wq,funct);
                valueVianLin(funi,ll)= g();
                if timings
                    timingsAVianLin(funi,ll)=timeit(g);
                end
            end
            nQuadsVianLin(:,ll)=length(xq);
            if timings
                timingsVianLin(ll)=timeit(f);
            end
    end
    for funi=1:length(functs)
        DAT{funi}(7)= {[nQuadsVianLin(funi,:); valueVianLin(funi,:); timingsVianLin; timingsAVianLin(funi,:);]};
    end
    % Exact and linear mesh the geometry, saved in filename.neu and
    % filenamelin.neu
    % Integrate over the exact and linear meshes using pth order Gaussian
    % quadrature
    nQuadsmesh= zeros(length(functs),maxmeshQuad); valuemesh= zeros(length(functs),maxmeshQuad);
    timingsQuadmesh=zeros(1,maxmeshQuad); timingsAQuadmesh=zeros(length(functs),maxmeshQuad);
    options.thresh=2;
    if maxmeshQuad>1
    fprintf('Attempting to exactly mesh the geometry \n');
    tic
    xmesh2d(filename,options);
    meshtime=toc;
    for i=1:maxmeshQuad
        tic
        fprintf('Performing exact rational meshing quadrature with %d Gauss points per triangle \n',i^2);
%         p=length(Xi{1})-length(W{1});
        for funi=1:length(functs)
            funct=functs{funi};
            f=@() integrate2dPts(filename,i);
            [xq,yq,wq]=f();
            nQuadsmesh(funi,i)=length(xq);
            g=@() applyRule(xq,yq,wq,funct);
            valuemesh(funi,i)=g();
            if timings
                timingsAQuadmesh(funi,i)=timeit(g);
            end
        end
%         refineMeshRefLvl(filename);
        if timings
            timingsQuadmesh(i)=timeit(f)+meshtime;
        end
    end
    end
    for funi=1:length(functs)
        DAT{funi}(1)={[nQuadsmesh(funi,:); valuemesh(funi,:);timingsQuadmesh; timingsAQuadmesh(funi,:);]};
    end
    
    
        
    nQuadsmeshlin= zeros(length(functs),maxquadmeshlin); valuemeshlin= zeros(length(functs),maxquadmeshlin);
    timingsQuadmeshlin=zeros(1,maxquadmeshlin);
    timingsAQuadmeshlin=zeros(length(functs),maxquadmeshlin);
    for i=1:maxquadmeshlin
        fprintf('Performing linear meshing quadrature with %d refinements \n',i-1);
        options.thresh=1 + 1*(10^(1-i));
        options.hmax=.3*(.5)^(1-i);
        tic
        [pts,tri]=xmesh2dlin(filename,options);
        meshtime=toc;
        for funi=1:length(functs)
            funct=functs{funi};
            f=@() integrate2dlinPts(pts,tri,p);
            [xq,yq,wq] = f();
            nQuadsmeshlin(funi,i)= length(xq); 
            g=@() applyRule(xq,yq,wq,funct);
            valuemeshlin(funi,i)=g();
            if timings
                timingsAQuadmeshlin(funi,i)=timeit(g);
            end
        end
        if timings
            timingsQuadmeshlin(i)=timeit(f)+meshtime;
        end
    end
   for funi=1:length(functs)
       DAT{funi}(2)={[nQuadsmeshlin(funi,:); valuemeshlin(funi,:); timingsQuadmeshlin; timingsAQuadmeshlin]};
   end

    
    for lpl=1:7
        if onoff(lpl)
    	for funi=1:length(functs)
            avgDAT{lpl}(2:4,:)=avgDAT{lpl}(2:4,:)+(abs(DAT{funi}{lpl}(2:4,:))-[abs(trueI(funi));0; 0])./(length(functs)*[abs(trueI(funi));1; 1]);
        end
        avgDAT{lpl}(1,:)=DAT{funi}{lpl}(1,:);
        end
    end
    
end

function [xq,yq,wq]=performVianQuadPts(CPW_my_format,kk,N_2009)
    xq=zeros(0,1); yq=zeros(0,1); wq=zeros(0,1);
    rotation=0; P=[0 0]; Q=[0 0]; cumulative=1;
    cubature_type=4; SPLtypestring='complete';
    numEvalspcurve=kk;evalPts=linspace(0,1-(1/numEvalspcurve),numEvalspcurve);
%         spline_order_vett{1}=[2*ones(4,1) (1:numEvalspcurve:(numEvalspcurve*4))'];
    for ii=1:length(CPW_my_format)
%             fst_spline_order=(bkpts{ii}-[1 bkpts{ii}(1:end-1)]);
%             fst_spline_order(fst_spline_order>0)=3; fst_spline_order=fst_spline_order+1;
%             spline_order_vett{ii}=[fst_spline_order' (kk*bkpts{ii}-kk+1)'];
       spline_order_vett{ii}=[repmat(4,size(CPW_my_format{ii},1)/3,1) ...
            (kk:kk:(kk*size(CPW_my_format{ii},1)/3))'+1];
        clear vertices;
        for jj=1:3:(size(CPW_my_format{ii},1))
            verts=dCR_eval(CPW_my_format{ii}(jj:(jj+2),:),evalPts);
            vertices((((jj-1)/3)*numEvalspcurve+1):((jj-1)/3+1)*numEvalspcurve,:)=verts;
        end
%             spline_order_vett{ii}=[4 length(vertices)+1];
        vertices(end+1,:)=vertices(1,:);
        f=@() splinegauss_2009b(N_2009,... 
            vertices,rotation,P,Q,spline_order_vett{ii},...
            cumulative,SPLtypestring,cubature_type);
        [x_nodes_2009,y_nodes_2009,weights_2009]=f();
        xq=[xq; x_nodes_2009]; yq=[yq; y_nodes_2009];
        wq=[wq; weights_2009];
    end
end

function [xq,yq,wq]=performVianQuadPtsLin(CPW_my_format,kk,N_2009)
    xq=zeros(0,1); yq=zeros(0,1); wq=zeros(0,1);
    rotation=0; P=[0 0]; Q=[0 0]; cumulative=1;
    cubature_type=4; SPLtypestring='complete';
    numEvalspcurve=kk;evalPts=linspace(0,1-(1/numEvalspcurve),numEvalspcurve);
%         spline_order_vett{1}=[2*ones(4,1) (1:numEvalspcurve:(numEvalspcurve*4))'];
    for ii=1:length(CPW_my_format)
%             fst_spline_order=(bkpts{ii}-[1 bkpts{ii}(1:end-1)]);
%             fst_spline_order(fst_spline_order>0)=3; fst_spline_order=fst_spline_order+1;
%             spline_order_vett{ii}=[fst_spline_order' (kk*bkpts{ii}-kk+1)'];
       spline_order_vett{ii}=[repmat(2,size(CPW_my_format{ii},1)/3,1) ...
            (kk:kk:(kk*size(CPW_my_format{ii},1)/3))'+1];
        clear vertices;
        for jj=1:3:(size(CPW_my_format{ii},1))
            verts=dCR_eval(CPW_my_format{ii}(jj:(jj+2),:),evalPts);
            vertices((((jj-1)/3)*numEvalspcurve+1):((jj-1)/3+1)*numEvalspcurve,:)=verts;
        end
%             spline_order_vett{ii}=[4 length(vertices)+1];
        vertices(end+1,:)=vertices(1,:);
        f=@() splinegauss_2009b(N_2009,... 
            vertices,rotation,P,Q,spline_order_vett{ii},...
            cumulative,SPLtypestring,cubature_type);
        [x_nodes_2009,y_nodes_2009,weights_2009]=f();
        xq=[xq; x_nodes_2009]; yq=[yq; y_nodes_2009];
        wq=[wq; weights_2009];
    end
end

function I = applyRule(xq,yq,wq,funct)
    I=wq'*funct(xq,yq);
end