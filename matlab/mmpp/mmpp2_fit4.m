function MAP=mmpp2_fit4(mean,scv,skew,acf1)
E1=mean;
E2=(1+scv)*E1^2;
if skew==-1
    E3=-1;
else
    E3=-(2*E1^3-3*E1*E2-skew*(E2-E1^2)^(3/2));
end
rho0=(1-1/scv)/2;
g2=acf1/rho0;
MAP=mmpp2_fit3(E1,E2,E3,g2);
if ~issym(E1)
    FEAS=map_isfeasible(MAP);
    if FEAS<12
        warning('model infeasible');
    end
end
end
