function [SCVj,acfCoeff]=kpcfit_sub_eval_iat2(SCV,G2,acfLag)
%  function [SCVj,acfCoeff]=kpcfit_sub_eval_iat2(SCV,G2,acfLag)
%
%   Input:
%       SCV = a vector of squared coefficient of variation for each
%             composing map
%       G2  = a vector of gamma2 for each composing map
%       acfLag
%           = a vector of autocorrelations lags
%
%   Output:
%       SCVj= the final squared coefficient of variation of J composing
%       maps
%       acfCoeff
%           = a vector of autocorrelation coefficients of J composing maps
%
J=length(G2);
SCVj=SCV(1);
acfLag=acfLag(:);
acfCoeff=0.5.*(1-1/SCV(1)).*G2(1).^acfLag;
for j=2:J
    SCVj_1=SCVj;
    SCVj=(1+SCVj)*(1+SCV(j))/2-1;
    r0j=0.5*(1-1/SCV(j));
    X=SCV(j)*r0j*G2(j).^acfLag;
    acfCoeff=(X+SCVj_1.*acfCoeff.*(1+X))./SCVj;
end

end
