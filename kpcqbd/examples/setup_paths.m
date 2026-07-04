function setup_paths()
%SETUP_PATHS Add the KPC Toolbox to the path for the KPC-QBD examples.
%   The KPC-QBD code is self-contained and does NOT require the LINE solver.

% kpc-toolbox root is two levels up from this examples/ folder
here    = fileparts(mfilename('fullpath'));   % .../kpcqbd/examples
kpcroot = fileparts(fileparts(here));         % .../kpc-toolbox
addpath(genpath(kpcroot));

% Suppress warnings during fitting
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

fprintf('Paths configured for KPC-QBD examples (LINE not required).\n');
end
