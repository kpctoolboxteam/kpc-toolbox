function MAP=mmpp2_fit2(mean,scv,skew,g2)
if scv==1
    MAP=mmpp2_exponential(mean);
    return
end

E1=mean;
E2=(1+scv)*E1^2;
E3=-(2*E1^3-3*E1*E2-skew*(E2-E1^2)^(3/2));
MAP=mmpp2_fit3(E1,E2,E3,g2);
if ~issym(E1)
    FEAS=map_isfeasible(MAP);
    if FEAS<12
        warning('model infeasible');
    end
end
end
