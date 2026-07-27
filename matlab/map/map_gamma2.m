function G2=map_gamma2(MAP)
if issym(MAP{1}(1,1))
    g=simplify(eig(map_embedded(MAP)));
    if g(1)==1
        G2=g(2);
    else
        G2=g(1);
    end
else
    P = map_embedded(MAP);       
    if issparse(P)
        eigP = eigs(P);
    else
        eigP = eig(P);
    end
    [~,J]=sort(abs(eigP),'descend');    
    G2=eigP(J(2));
end
end