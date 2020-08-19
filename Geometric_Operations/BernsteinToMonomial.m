function monCoeffs = BernsteinToMonomial(bernCoeffs)
    n=length(bernCoeffs)-1;
    monCoeffs=bernCoeffs;
    for k=(n:-1:1)
        for j=k:n
            monCoeffs(j+1)= monCoeffs(j+1)-monCoeffs(j);
        end
    end
    monCoeffs=monCoeffs.*(arrayfun(@(kk)nchoosek(n,kk),0:n));
end