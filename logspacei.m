function iset=logspacei(a,b,points)
%LOGSPACEI(a,b,numpoints) - logarithmically spaced intergers in [10^a,10^b]
    iset=round(logspace(log10(a),log10(b),points));
    iset(iset<a)=a;
    iset(iset>b)=b;
end