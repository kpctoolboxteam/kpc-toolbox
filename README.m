% MAP Queueing Networks Toolbox
% Version 1.0, 15-Apr-2008, last release available at [1]
%
% MAP PROCESSES: Analysis
% map_acf         -  Compute autocorrelation coefficients
% map_cdf         -  Compute cumulative distribution function
% map_embedded    -  Embedded discrete-time process
% map_idc         -  Compute asymptotic index of dispersion
% map_infgen      -  Infinitesimal generator of the underlying continuous-time process
% map_isfeasible  -  Evaluate feasibility of a MAP process
% map_kurt        -  Compute kurtosis
% map_lambda      -  Compute mean arrival rate
% map_mean        -  Compute mean interarrival time
% map_moment      -  Compute moments of interarrival times
% map_pie         -  Equilibrium distribution of the embedded discrete-time process
% map_prob        -  Equilibrium distribution of the underlying continuous-time process
% map_sample      -  Generate a random sample of inter-arrival times
% map_scv         -  Compute squared coefficient of variation
% map_skew        -  Compute skewness
% map_var         -  Compute variance
%
% MAP PROCESSES: Fitting
% map_erlang      -  Fit an Erlang-K process as a MAP
% map_exponential -  Fit a Poisson process as a MAP
% map_hyperexp    -  Fit a two-phase Hyperexponential process as a MAP
% map_mmpp2	      -  Fit a MMPP(2) process as a MAP
%
% MAP PROCESSES: Manipulation 
% map_normalize   -  Try to make a MAP feasible
% map_renewal     -  Remove all correlations
% map_scale       -  Rescale mean inter-arrival time
%
% UTILITY:
% dtmc_solve      -  Equilibrium distribution of a discrete-time Markov chain
% ctmc_solve      -  Equilibrium distribution of a continuous-time Markov chain
%
% MAP PROCESS LIBRARY
% maplib_bcaug89  - MAP(16) process obtained by KPC fitting of Bellcore Aug89 Trace, see [2]
% maplib_circul4  - MAP(4) with circulant embedded process (complex eigenvalues, obscillating autocorrelation)
% maplib_erlang2  - Erlang-2 process, mean=1
% maplib_erlang3  - Erlang-3 process, mean=1
% maplib_exp      - Poisson process, mean=1
% maplib_hyper4   - Hyperexponential process, mean=1, scv=4, p=0.99, no correlations
% maplib_hyper25  - Hyperexponential process, mean=1, scv=25, p=0.99, no correlations
% maplib_hmap4    - Hyperexponential MMPP(2), mean=1, scv=4, p=0.99, weak correlations (lag-1 autocorrelation coefficient = 0.1, quick decay rate)
% maplib_hmap25   - Hyperexponential MMPP(2), mean=1, scv=25, p=0.99, strong correlations (lag-1 autocorrelation coefficient = 0.48, very slow decay rate)
% maplib_saw      - MAP(2) process with saw-like autocorrelation function (-1 eigenvalue in embedded process matrix)
%
% MAP QUEUEING NETWORKS: Analysis
% mapqn_ezsolve   - Define and solve a MAP queueing network
% mapqn_hashstate - Determine position of a state in the equilibrium probability vector of the MAP queueing network
%
% References:
% [1] http://www.cs.wm.edu/MAPQN
% [2] G.Casale, E.Zhang, and E.Smirni. Interarrival times characterization  and fitting for markovian traffic analysis, 2008.
%     Department of Computer Science, College of William and Mary, Tech. Rep. WM-CS-2008-02. Available at:
%     http://www.wm.edu/computerscience/techreport/2008/WM-CS-2008-02.pdf

% *************************************************
% *************************************************
% *************************************************
% type "help README" to display this help in MATLAB 
% *************************************************
% *************************************************
% *************************************************
