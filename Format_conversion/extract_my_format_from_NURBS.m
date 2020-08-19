function CPWmat = extract_my_format_from_NURBS(CP,W,Xi)
CPWmat=cell(1,length(CP));
for segment=1:length(CP)
    p=length(Xi{segment})-length(W{segment})-1;
    n=length(W{segment});
    [n_el,C,IEN2]=Extract_Basis_1D(p,n,Xi{segment});
    CPWmat{segment}=zeros(n_el*3,p+1);
    for i=1:n_el
        CPWmat{segment}(3*i,:)=W{segment}(:,IEN2(:,i))*C(:,((p+1)*(i-1)+1):((p+1)*i));
        CPWmat{segment}(((i-1)*3+1):(3*i-1),:)=CPWmat{segment}(3*i,:).*CP{segment}(:,IEN2(:,i))*C(:,((p+1)*(i-1)+1):((p+1)*i));
    end
end