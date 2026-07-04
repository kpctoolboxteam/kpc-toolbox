function pa = kpcqbd_solve(ARV,H,npa)
% Calculating pi as pi(i+1) = pi(i)*R where R is calculated at every step
% using the Woodbury identity 

% ARV = arrival process
% H = service MAPs
% npa = number of probability terms to compute

if nargin<3
    error('The approximation requires at least 3 parameters');
end

pa = zeros(npa,1);

% create service process
K = length(H);
ns=1;
mu=1;
S0 = 1;
for k=1:K
    S0 = kron(S0,sparse(H{k}{1}));
    ns = ns * length(H{k}{1});
    mu = mu * map_lambda(H{k});
end
if S0(1,1)>0
    S0 = -S0;
end
lambda = map_lambda(ARV);
rho = lambda / mu;

k=1;
H{k} = map_scale(H{k},rho/map_lambda(ARV));
Fkpc{k} = kron(ARV{2},eye(length(H{1}{2})));
Lkpc{k} = krons(ARV{1},H{k}{1});
Bkpc{k} = kron(eye(length(ARV{2})),H{k}{2});
Gkpc{k} = qbd_gmatrix(Bkpc{k}, Lkpc{k}, Fkpc{k});
alpha{k} = 0;
p0 = dtmc_solve(Gkpc{k});


for k=2:K
    H{k}=map_scale(H{k},rho);
    alpha{k} = map_pie(H{k});
    p0 = kron(p0, alpha{k}); %calculating pi_0
end

pa(1) = 1-rho;


% Decompose F*G as U*S*V (this is not the SVD decomposition but S is still diag and U and V are rectangular)
% We are using here that the arrival and the services apart from H{1} are
% PH distributions
% [U, S, V] = svds(Fkpc{1}*Gkpc{1}, length(H{1}{1}));
[U, S, V] = svds(Fkpc{1}*Gkpc{1}, rank(Fkpc{1}*Gkpc{1}));
V=V';

for k=2:K
    U = kron(U,ones(length(H{k}{1}),1));
    V = kron(V,alpha{k});
end

iC = diag(sparse(1)./diag(S));

F = kron(sparse(ARV{2}),speye(ns));
L = krons(sparse(ARV{1}),S0);
if sum(sum(abs(L-diag(diag(L))))) < 1e-8 % is diagonal?
    iL = diag(sparse(1)./diag(L));
else
    iL = inv(L); 
end
ViL = V*iL;
iCViLU = iC+ViL*U;

palast = p0*(1-rho);
for i=1:(npa-1)
    u0 = palast*(-F);
    u0=u0*iL;
    u1 =u0*U;
    u1t = (iCViLU'\u1')';
    u1v = u1t*V;
    u1=u1v*iL;
    pacur = u0 - u1; % woodbury rank-m update
    pa(1+i) = sum(pacur,2);
    palast=pacur;
end

end