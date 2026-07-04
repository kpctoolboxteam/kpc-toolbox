function G = qbd_gmatrix(A0, A1, A2)
% G = QBD_GMATRIX(A0,A1,A2) - minimal nonnegative solution of the QBD
% matrix equation
%
%        0 = A0 + A1*G + A2*G^2
%
% for a continuous-time Quasi-Birth-Death process with generator blocks
% A0 (down), A1 (local), A2 (up), where A0+A1+A2 has zero row sums.
%
% Self-contained cyclic reduction [Bini,Meini], used in place of the LINE
% routine QBD_CR so that the KPC-QBD code carries no external dependency.

m = size(A1,1);
I = eye(m);

% uniformize the continuous-time blocks into an equivalent discrete-time
% QBD (G is invariant under uniformization)
if any(diag(A1) < 0)
    lambda = max(-diag(A1));
    A0 = A0/lambda;
    A1 = A1/lambda + I;
    A2 = A2/lambda;
end

% basic cyclic reduction
A = A1; B = A2; C = A0; Ahat = A1;
check = 1; numit = 0; maxit = 50;
while check > 1e-14 && numit < maxit
    T  = (I - A) \ I;
    BT = B * T;
    CT = C * T;
    Ahat = Ahat + BT * C;
    A    = A + BT * C + CT * B;
    B    = BT * B;
    C    = CT * C;
    numit = numit + 1;
    check = min(norm(B, inf), norm(C, inf));
end
if numit == maxit && check > 1e-14
    warning('qbd_gmatrix:maxIter', ...
        'cyclic reduction stopped after %d iterations (residual %g)', maxit, check);
end

G = (I - Ahat) \ A0;
end
