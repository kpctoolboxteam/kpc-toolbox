function setup_paths()
%SETUP_PATHS Add the KPC Toolbox to the path for the KPC-QBD examples.

here    = fileparts(mfilename('fullpath'));
kpcroot = fileparts(fileparts(here));
addpath(genpath(kpcroot));

warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

fprintf('Paths configured for KPC-QBD examples.\n');
end
