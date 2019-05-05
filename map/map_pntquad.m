function [Pnt,S]=map_pntquad(MAP,na,t)
D0=MAP{1};
D1=MAP{2};
Ki=length(MAP{1});
P0=zeros(Ki^2*(na+1),1);
P0(1:Ki^2)=reshape(eye(Ki),1,Ki^2);
for n=1:na
    P0(n*Ki^2+1:(n+1)*Ki^2)=0*eye(Ki);
end
options = odeset('RelTol',1e-10,'MaxStep',1e-3,'AbsTol',1e-10,'NonNegative',1:length(P0),'Refine',10);
[t,P]=ode45(@pnt_ode,[0 t],P0,options);
P=P(end,:);
Pnt=[];
S=zeros(Ki);
for n=0:na
    Pnt{n+1}=reshape(P(n*Ki^2+1:(n+1)*Ki^2),Ki,Ki);
    S=S+Pnt{n+1};
end
Pnt=Pnt{end};
    function dP = pnt_ode(t,P)
        na=round(length(P)/Ki^2)-1;
        dP(1:Ki^2,1)=reshape(reshape(P(1:Ki^2),Ki,Ki)*D0,Ki^2,1);
        for n=1:na
            Pn_1t=reshape(P(((n-1)*Ki^2+1):n*Ki^2),Ki,Ki);
            Pnt=reshape(P((n*Ki^2+1):(n+1)*Ki^2),Ki,Ki);
            dP((n*Ki^2+1):(n+1)*Ki^2,1)=reshape(Pnt*D0+Pn_1t*D1,Ki^2,1);
        end
    end
end
