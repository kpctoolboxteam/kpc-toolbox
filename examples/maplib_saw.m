function MAP=maplib_saw()
% MAP=maplib_saw() -  MAP(2) process with saw-like autocorrelation function
% (-1 eigenvalue in embedded process matrix)
%
%  Output:
%  MAP: pre-fitted MAP process
%

D0=[-1     0
     0    -2];
D1=[ 0     1
     2     0];
MAP=map_scale({D0,D1},1);
end