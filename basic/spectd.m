function [spectrum,projectors,nihil,V,D]=spectd(A,OPT)
if nargin<2
    JORDAN=0;
elseif strcmpi(OPT,'jordan')
    JORDAN=1;
else
    JORDAN=0;
end
if JORDAN
    [V,J]=jordan(A);
    iV=inv(V);
    n=length(A);
    from=1;
    if issym(A(1,1))
        spectrum=sym([]);
        projectors={};
        nihil={};
    else
        spectrum=[];
        projectors={};
        nihil={};
    end
    for to=1:n
        if to==n || J(to,to+1)==0 % next is new block
            spectrum(end+1)=J(from,from);
            projectors{end+1}=V*subblock(sym(eye(n)),from,to)*iV;
            nihil{end+1}=V*subblock(diag(ones(1,n-1),1),from,to)*iV;
            from=to+1;
        end
    end
    D=J;
else
    if ~strcmpi(class(A(1,1)),'sym')
        if ~issparse(A)
            [V,D]=eig(A);
        else
            [V,D]=eigs(A);
        end
        spectrum=diag(D);
        projectors={};
    else
        [V,D]=eig(A);
        spectrum=diag(D);
        projectors={};
        i=1:length(A);
    end
    i=1:length(A);
    spectrum=spectrum(i);
    iV=inv(V);
    for k=i(:)'
        projectors{end+1}=V(:,k)*iV(k,:);
    end
    nihil={};
end
spectrum=reshape(spectrum,1,length(spectrum));

    function X=subblock(X,from,to)
        X(1:(from-1),:)=0;
        X(:,1:(from-1))=0;
        X((to+1):end,:)=0;
        X(:,(to+1):end)=0;
    end

end