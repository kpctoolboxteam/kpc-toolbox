% MAP Queueing Networks Toolbox
% Version 1.0 	 15-Apr-2008
%
% MAP PROCESSES: Analysis
% map_acf         -  Compute autocorrelation coefficients
% map_cdf         -  Compute cumulative distribution function
% map_embedded    -  Embedded discrete-time process
% map_idc         -  Compute asymptotic index of dispersion
% map_infgen      -  Infinitesimal generator of the underlying continuous-time process
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
% map_isfeasible  -  Evaluate feasibility of a MAP process
% map_normalize   -  Try to make a MAP feasible
% map_renewal     -  Remove all correlations
% map_scale       -  Rescale mean inter-arrival time
%
% UTILITY:
% dtmc_solve      -  Equilibrium distribution of a discrete-time Markov chain
% ctmc_solve      -  Equilibrium distribution of a continuous-time Markov chain
%
% MAP QUEUEING NETWORKS
% mapqn_ezsolve   -  Define and solve a MAP queueing network
%
% WEB PAGE: http://www.cs.wm.edu/MAPQN