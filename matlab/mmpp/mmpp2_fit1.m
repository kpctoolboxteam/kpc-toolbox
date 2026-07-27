function MAP=mmpp2_fit1(mean,scv,skew,idc)
% MAP=mmpp2_fit1(mean,scv,skew,idc)
E1=mean;
E2=(1+scv)*E1^2;
g2=-(scv-idc)/(-1+idc);
if skew==-1
    E3=-1;
else
    E3=-(2*E1^3-3*E1*E2-skew*(E2-E1^2)^(3/2));
end
MAP=map2_fit(E1,E2,E3,g2);
%MAP=mmpp2_fit3(E1,E2,map_e3(MAP),g2);
end
