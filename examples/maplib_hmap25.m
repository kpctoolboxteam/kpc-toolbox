function MAP=maplib_hmap25()
% MAP=maplib_hmap25() -  Hyperexponential MMPP(2), mean=1, scv=25, p=0.99,
% strong correlations (lag-1 autocorrelation coefficient = 0.48, very slow
% decay rate)
%
%  Output:
%  MAP: pre-fitted MAP process
%
% MAP Queueing Networks Toolbox
% Version 1.0 	 15-Apr-2008
MAP=map_mmpp2(1,25,-1,-1);
end