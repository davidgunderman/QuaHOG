function [l,Bas] = ComputeSplineBasis(xi,p,Xi)
    if xi>=0.0 && xi<1.0
        l=which_knot(xi,Xi);
%     offset=p-l+1;
        N=zeros(p,length(Xi)-1);
        N(1,l)=1;
        for k=2:(p+1)
            for i=1:(length(Xi)-k)
    %             if i>=(1)
    %                 Bik = (xi-Xi(i))/(Xi(i+k)-Xi(i));
    %                 N(i,k) = Bik*N(i,k-1);
    %             elseif i<=0
    %                 Bi1k = (xi-Xi(i+1))/(Xi(i+1+k)-Xi(i+1));
    %                 N(i,k) = (1-Bi1k)*N(i+1,k-1);
    %             else
                    Bik = (xi-Xi(i))/(Xi(i+k-1)-Xi(i));
                    Bi1k = (xi-Xi(i+1))/(Xi(i+1+(k-1))-Xi(i+1));
                    if isnan(Bik)||Bik==-Inf || Bik==Inf
                        Bik=0;
                    end
                    if isnan(Bi1k)||Bi1k==-Inf || Bi1k==Inf
                        Bi1k=0;
                    end
                    N(k,i) = (1-Bi1k)*N(k-1,i+1)+Bik*N(k-1,i);
    %             end
            end
        end
    %     Bas=zeros(1,length(Xi)-p-1)
        Bas=N(end,:);
    else
        l=0;
        Bas=zeros(1,length(Xi)-p-1);
    end
end