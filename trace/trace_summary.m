function [MEAN,SCV,MAD,SKEW,KURT,QUART,TAILS1PERC,MINMAX,IQR,ACF,IDC]=trace_summary(m,fid)
if nargin==1
    fid=1;
end
MEAN=mean(m);
SCV=scv(m);
SKEW=skewness(m);
KURT=kurtosis(m);
QUART=[prctile(m,25),prctile(m,50),prctile(m,75)];
TAILS1PERC=[prctile(m,95),prctile(m,(1-1e-6)*100)];
MINMAX=[min(m),max(m)];
MAD=mad(m,1); %median based mad
ACF=trace_acf(m,1:10).';
IQR=iqr(m);
IDC=trace_idc(m);
LEN =length(m);
fprintf(fid,'length=%d mean=%f scv=%f cv=%f mad=%f, skew=%f kurt=%f p25=%f p50=%f p75=%f p95=%f min=%f max=%f iqr=%d ',LEN,MEAN,SCV,sqrt(SCV),MAD,SKEW,KURT,QUART(1),QUART(2),QUART(3),TAILS1PERC(1),MINMAX(1),MINMAX(2),IQR);
%fprintf(fid,'iqr/med=%f min=%f max=%f ',length(m),MEAN,SCV,sqrt(SCV),MAD,SKEW,IQR/QUART(2),MINMAX(1),MINMAX(2));
fprintf(fid,'acf=%f %f %f %f idc(1000)/scv=%f',ACF(1),ACF(2),ACF(3),ACF(4),IDC/SCV);
fprintf(fid,'\n');
end
