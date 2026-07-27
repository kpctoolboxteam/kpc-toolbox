function E = mvph_joint(alpha,S,T,D,n1,n2)
% E = mvph_joint(alpha,X,T,Y,n1,n2)
% Joint moment E[X^n1*Y^n2] of a bivariate phase-type distribution with
% transition matrxi D between X and Y and initial vector alpha

E = factorial(n1) * factorial(n2) * alpha * inv(-S)^(n1+1) * D * inv(-T)^(n1+1) * (-T) * ones(size(T,1),1);

end