function MAP=maplib_hmap4()
% MAP=maplib_hmap4() -  Hyperexponential MMPP(2), mean=1, scv=4, p=0.99,
% weak correlations (lag-1 autocorrelation coefficient = 0.1, quick decay
% rate)
%
%  Output:
%  MAP: pre-fitted MAP process
%

MAP=map_mmpp2(1,4,-1,0.1);
end