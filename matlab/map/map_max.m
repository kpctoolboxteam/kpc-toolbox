function M=map_max(A,B)
na=size(A{1},1);
nb=size(B{1},1);
a=-A{1}*ones(na,1);
b=-B{1}*ones(nb,1);
M0=[
    krons(A{1},B{1}),kron(a,eye(na)),kron(b,eye(nb))
    zeros(nb,na*nb),B{1},zeros(nb,nb)
    zeros(na,na*nb),zeros(na,nb),A{1}
    ];
pie = [kron(map_pie(A),map_pie(B)),zeros(1,na),zeros(1,nb)];
M1=pie'*ones(1,size(M0,1));
M={M0,M1};
end