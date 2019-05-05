function I=map_idc(MAP)
% I=map_idc(MAP) - Compute the asymptotic index of dispersion. I equals
%  
%  I=SCV(1+2\sum_{k=1}^{\infty} \rho_k)
%
%  where SCV is the squared coefficient of variation and \rho_k is the
%  lag-k autocorrelation coefficient of inter-arrival times. I is also the 
%  limiting value of the index of dispersion for counts and of the index of 
%  dispersion for intervals.
%
%  Input: 
%  MAP: a MAP in the form of {D0,D1}
%       
%  Output: 
%  I: asymptotic index of dispersion
%
%  Examples:
%  - I=map_idc(map_renewal(MAP)) equals the SCV of the MAP
%


e=ones(size(MAP{1},1),1);
I=1+2*(map_lambda(MAP)-map_pie(MAP)*inv(map_infgen(MAP)+e*map_prob(MAP))*MAP{2}*e);
end