function [XN,QN,UN,pqueue] = kpcqbd_mapmap1(MAPa, MAPs)
% [XN,QN,UN,PQUEUE] = KPCQBD_MAPMAP1(MAPA,MAPS) - exact analysis of a
% MAP/MAP/1 queue via the matrix-geometric method.
%
% MAPA    - arrival process {D0,D1}
% MAPS    - service process {D0,D1}
% XN      - throughput (arrival rate)
% QN      - mean queue length (number in system)
% UN      - utilization
% PQUEUE  - queue-length distribution, one row per level, columns indexed
%           by the (arrival x service) phase; row n+1 is level n
%
% Self-contained replacement for the LINE routine qbd_mapmap1 (which relies
% on QBD_CR/QBD_pi); the QBD G-matrix is obtained from qbd_gmatrix.

na = length(MAPa{1});
ns = length(MAPs{1});
m  = na*ns;

A2 = kron(MAPa{2}, eye(ns));    % up   : an arrival occurs
A1 = krons(MAPa{1}, MAPs{1});   % local: repeating level
A0 = kron(eye(na), MAPs{2});    % down : a service completes
L0 = kron(MAPa{1}, eye(ns));    % local at boundary level 0 (no service)

lambda = map_lambda(MAPa);
mu     = map_lambda(MAPs);

% G and R matrices of the repeating portion
G = qbd_gmatrix(A0, A1, A2);          % 0 = A0 + A1 G + A2 G^2
U = A1 + A2*G;
R = A2 / (-U);                        % R = A2*(-U)^{-1}, minimal nonnegative

% boundary vector [pi0 pi1] from the level-0 and level-1 balance equations
%   pi0*L0 + pi1*A0            = 0
%   pi0*A2 + pi1*(A1 + R*A0)   = 0
M  = [L0, A2; A0, A1 + R*A0];
[~,~,V] = svd(full(M'));       % smallest right singular vector = null space of M'
xi = V(:,end)';
pi0 = xi(1:m);
pi1 = xi(m+1:end);

% normalize:  pi0*1 + pi1*(I-R)^{-1}*1 = 1
tot = sum(pi0) + sum(pi1 * ((eye(m)-R)\ones(m,1)));
pi0 = pi0/tot;
pi1 = pi1/tot;
if sum(pi0) < 0            % fix sign of the null vector
    pi0 = -pi0; pi1 = -pi1;
end

% expand the matrix-geometric tail pi_n = pi1*R^(n-1) until the mass is captured
maxLevels = 20000;
pqueue = zeros(maxLevels, m);
pqueue(1,:) = pi0;
row = pi1;
captured = sum(pi0);
n = 1;
while n < maxLevels
    pqueue(n+1,:) = row;
    captured = captured + sum(row);
    if captured > 1-1e-12 && n >= 2
        n = n+1; break
    end
    row = row * R;
    n = n+1;
end
pqueue = pqueue(1:n,:);

XN = lambda;
UN = 1 - sum(pqueue(1,:));
QN = (0:size(pqueue,1)-1) * sum(pqueue,2);
end
