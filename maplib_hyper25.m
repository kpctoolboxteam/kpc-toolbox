function MAP=maplib_hyper25()
% MAP=maplib_hyper25() -  Hyperexponential process, mean=1, scv=25, p=0.99,
% no correlations
%
%  Output:
%  MAP: pre-fitted MAP process
%
% MAP Queueing Networks Toolbox
% Version 1.0 	 15-Apr-2008
MAP=map_hyperexp(1,25,0.99);
end