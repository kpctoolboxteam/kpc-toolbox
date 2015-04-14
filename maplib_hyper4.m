function MAP=maplib_hyper4()
% MAP=maplib_hyper4() -  Hyperexponential process, mean=1, scv=4, p=0.99, no
% correlations
%
%  Output:
%  MAP: pre-fitted MAP process
%
% MAP Queueing Networks Toolbox
% Version 1.0 	 15-Apr-2008
MAP=map_hyperexp(1,4,0.99);
end