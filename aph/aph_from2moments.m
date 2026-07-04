function MAP = aph_from2moments(MEAN, SCV)
% MAP = APH_FROM2MOMENTS(MEAN,SCV) - acyclic phase-type distribution, in
% renewal MAP form {D0,D1}, matching a target mean and squared coefficient
% of variation SCV (=variance/mean^2).
%
% Canonical acyclic construction of [Bobbio,Horvath,Telek].

cv2    = SCV;
lambda = 1/MEAN;
N      = max(ceil(1/cv2), 2);
p      = 1/(cv2 + 1 + (cv2-1)/(N-1));

T = -lambda*p*N*eye(N);
for i=1:N-1
    T(i,i+1) = -T(i,i);
end
T(N,N) = -lambda*N;

alpha = zeros(1,N);
alpha(1) = p;
alpha(N) = 1 - p;

d   = -T*ones(N,1);
MAP = {T, d*alpha};
end
