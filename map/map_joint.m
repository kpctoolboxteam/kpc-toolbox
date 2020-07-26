function JM=map_joint(MAP,a,i)
% map_joint computes the joint moments of a MAP
%
% JM=map_joint(MAP,a,i) returns the joint moment of a MAP
% 
%   Input:      
%           MAP = a map in the form of {D0,D1}
%           a   = a vector (a1, a2, a3, ... ) which are the subscripts 
%                 of each term in E[(X_a1)^i1*(X_a2)^i2*...]
%           i   = a vector (i1, i2, i3, ... ) which specifying the 
%                 power of each element in the joint moments 
%                 E((X_a1)^i1*(X_a2)^i2*(X_a3)^i3*(X_a4)^i4*...]
%
%   Output: 
%           JM = E[(X_a1)^i1*(X_{a1+a2})^i2*(X_{a1+a2+a3})^i3*... ]
%
%
%  
%
a=cumsum(a);                 % the cumulative vector of vector a
P=inv(-MAP{1})*MAP{2};
invD0=inv(-MAP{1});
JM=1;
K=length(a);
for k=1:(K-1)
    JM=JM*factorial(i(k))*invD0^(i(k))*(P^(a(k+1)-a(k)));
end
JM=map_pie(MAP)*JM*factorial(i(K))*invD0^(i(K))*ones(length(P),1);
end