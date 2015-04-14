function Q=map_infgen(MAP)
% Q=map_infgen(MAP) -  Infinitesimal generator of the underlying
%  continuous-time process of the MAP describing the state dynamics before
%  absorption
%
%  Input: 
%  MAP: a MAP in the form of {D0,D1}
%       
%  Output: 
%  Q: infinitesimal generator (=D0+D1)
%

Q=MAP{1}+MAP{2};