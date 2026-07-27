function[acv] = autocov(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% autocov computes the autocovariance between two column vectors X and Y
% with same length N using the Fast Fourier Transform algorithm from 0 to
% N-2.
% The resulting autocovariance column vector acv is given by the formula:
% acv(p,1) = 1/(N-p) * \sum_{i=1}^{N} (X_{i} - X_bar) * (Y_{i+p} - Y_bar)
% where X_bar and Y_bar are the mean estimates: 
% X_bar = 1/N * \sum_{i=1}^{N} X_{i}; Y_bar = 1/N * \sum_{i=1}^{N} Y_{i}
% It satisfies the following identities:
% 1. variance consistency: if acv = autocov(X,X), then acv(1,1) = var(X)
% 2. covariance consistence: if acv = autocov(X,Y), then acv(1,1) = cov(X,Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Jacques Burrus on March 29th 2012.
% Url: www.bfi.cl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y=X;

% Verify the input consistency
if nargin < 1, error('Missing input vector.'); end

[M N] = size(X);
[P Q]= size(Y);
if M < 2, error('X is too short.'); end
if M~=P || N~=Q, error('Input vectors do not have the same size.'); end
if N ~= 1, error('X must be a column vector.'); end

% Compute the autocovariance
X_pad = [X - mean(X); zeros(M,1)];%%paddle
Y_pad = [Y - mean(Y); zeros(M,1)];

X_hat = fft(X_pad);
Y_hat = fft(Y_pad);

acv = ifft( conj( X_hat ) .* Y_hat );%the imaginary part is due to float precision errors.
acv = real( acv(1:M-1) ) ./ (M - (1:1:M-1))';