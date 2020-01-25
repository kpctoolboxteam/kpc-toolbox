function [MAP,ERR]=map2_fit(e1,e2,e3,g2)
% [MAP,ERR]=map2_fit(e1,e2,e3,g2)
% A. Heindl, G.Horvath, K. Gross "Explicit inverse characterization of
% acyclic MAPs of second order"
if nargin==3
    g2=e3;
    e3=-1;
end
ERR=0;
r1=e1; r2=e2/2;
h2=(r2-r1^2)/r1^2;
if e3==-1
    % select e3 that maximizes the range of correlations
    scv=(e2-e1^2)/e1^2;
    if 1<=scv && scv<3
        if g2<0
            h3=h2-h2^2; %b=0
            e3=12*e1^3*h2+6*e1^3*h3+6*e1^3*(1+h2^2);            
        else
            e3=(3/2+1e-3)*e2^2/e1;
        end
    elseif 3<=scv
        e3=(3/2+1e-3)*e2^2/e1;
    elseif 0<scv && scv<1
        e3=(1+1e-10)*(12*e1^3*h2+6*e1^3*(h2*(1-h2-2*sqrt(-h2)))+6*e1^3*(1+h2^2));        
    end
end
if e3==-2
    % select the minimum e3
    scv=(e2-e1^2)/e1^2;
    if 1<=scv 
        e3=(3/2+1e-6)*e2^2/e1;
    elseif 0<scv && scv<1
        h3=h2*(1-h2-2*sqrt(-h2));
        e3=6*e1^3*(h2^2 + h3);        
    end
end
if e3==-3
    % select the maximum e3
    scv=(e2-e1^2)/e1^2;
    if 1<=scv 
        e3=10^6;
    elseif 0<scv && scv<1
        h3=(-h2)^2;
        e3=6*e1^3*(h2^2 + h3);        
    end
end
if e3==-4
    % select a random e3
    scv=(e2-e1^2)/e1^2;
    r = rand;
    if 1<=scv 
        e3=r*(3/2+1e-6)*e2^2/e1+(1-r)*10^6;        
    elseif 0<scv && scv<1
        h3=r*(-h2)^2+(1-r)*h2*(1-h2-2*sqrt(-h2));
        e3=6*e1^3*(h2^2 + h3);        
    end
end
if e3>-1 && e3<0    
    % use a custom random e3
    scv=(e2-e1^2)/e1^2;
    r = abs(e3);    
    if 1<=scv 
        e3=r*(3/2+1e-6)*e2^2/e1+(1-r)*10^6;        
    elseif 0<scv && scv<1
        h3=r*h2*(1-h2-2*sqrt(-h2))+(1-r)*(-h2)^2;
        e3=6*e1^3*(h2^2 + h3);        
    end
end
r3=e3/6;
h3=(r3*r1-r2^2)/r1^4;
b=h3+h2^2-h2;
c=sqrt(b^2+4*h2^3);

if r1<=0
    MAP={};
    ERR=10; % mean out of bounds
    return
end

if h2==0
    if h3==0 && g2==0
        MAP=map_exponential(e1);
    else
        MAP={};
        ERR=20; % correlated exponential
        return
    end
end

if h2>0 && h3>0
    MAP=map_fit_hyper();
    if length(MAP)>0 && map_isfeasible(MAP)==0 
        ERR=-1;
    end
    return
elseif -1/4<=h2 && h2<0 && h2*(1-h2-2*sqrt(-h2))<=h3 && h3<=-h2^2
    MAP=map_fit_hypo();
    if length(MAP)>0 && map_isfeasible(MAP)==0
        ERR=-1;
    end
    return
else
    if ~(-1/4<=h2 && h2<0 && h2*(1-h2-2*sqrt(-h2))<=h3 && h3<=-h2^2)
        ERR=30; % h2 out of bounds
        MAP={};
        return
    elseif (h2>0 && h3<0) || h2*(1-h2-2*sqrt(-h2))>h3 || h3<=-h2^2
        ERR=40; % h3 out of bounds
        MAP={};
        return
    else
        error('I lost an error')
    end
end

    function MAP=map_fit_hyper()
        if b>=0
            if (b-c)/(b+c)<=g2 && g2<1
                MAP{1}=(1/(2*r1*h3))*[-(2*h2+b-c),0;0,-(2*h2+b+c)];
                MAP{2}=(1/(4*r1*h3))*[(2*h2+b-c)*(1-b/c+g2*(1+b/c)),(2*h2+b-c)*(1+b/c)*(1-g2); (2*h2+b+c)*(1-b/c)*(1-g2),(2*h2+b+c)*(1+b/c+g2*(1-b/c))];
                return
            else
                ERR=51; % g2 out of bounds
                MAP={};
                return
            end
        elseif b<0
            if 0<=g2 && g2<1
                MAP{1}=(1/(2*r1*h3))*[-(2*h2+b-c),0;0,-(2*h2+b+c)];
                MAP{2}=(1/(4*r1*h3))*[(2*h2+b-c)*(1-b/c+g2*(1+b/c)),(2*h2+b-c)*(1+b/c)*(1-g2); (2*h2+b+c)*(1-b/c)*(1-g2),(2*h2+b+c)*(1+b/c+g2*(1-b/c))];
                return
            elseif -(h3+h2^2)/h2 <= g2 && g2<0
                a=(h3+h2^2)/h2;
                d1=((1-a)*(2*h2*g2+b-c)+g2*(b+c)-(b-c))/((1-a)*(2*h2+b-c)+2*c);
                d2=((g2-1)*(b-c))/((1-a)*(2*h2+b-c)+2*c);
                MAP{1}=(1/(2*r1*h3))*[-(2*h2+b-c),(2*h2+b-c)*(1-a);0,-(2*h2+b+c)];
                MAP{2}=(1/(2*r1*h3))*[(2*h2+b-c)*d1,(2*h2+b-c)*(a-d1); (2*h2+b+c)*d2, (2*h2+b+c)*(1-d2)];
                return
            else
                ERR=52; % g2 out of bounds
                MAP={};
                return
            end
        end
    end
    function MAP=map_fit_hypo()
        if g2>=0
            if g2<=-(h2+sqrt(-h3))^2/h2
                a=(2*h2+b-c)*(h2+sqrt(-h3))/(2*h2*sqrt(-h3));
                c=-c;
                d1=((1-a)*(2*h2*g2+b-c)+g2*(b+c)-(b-c))/((1-a)*(2*h2+b-c)+2*c);
                d2=((g2-1)*(b-c))/((1-a)*(2*h2+b-c)+2*c);
                MAP{1}=(1/(2*r1*h3))*[-(2*h2+b-c),(2*h2+b-c)*(1-a);0,-(2*h2+b+c)];
                MAP{2}=(1/(2*r1*h3))*[(2*h2+b-c)*d1,(2*h2+b-c)*(a-d1); (2*h2+b+c)*d2, (2*h2+b+c)*(1-d2)];
            else
                ERR=53; % g2 out of bounds
                MAP={};
                return
            end
        elseif g2<0
            if g2 >= -(h3+h2^2)/h2
                a=(h3+h2^2)/h2;
                c=-c;
                d1=((1-a)*(2*h2*g2+b-c)+g2*(b+c)-(b-c))/((1-a)*(2*h2+b-c)+2*c);
                d2=((g2-1)*(b-c))/((1-a)*(2*h2+b-c)+2*c);
                MAP{1}=(1/(2*r1*h3))*[-(2*h2+b-c),(2*h2+b-c)*(1-a);0,-(2*h2+b+c)];
                MAP{2}=(1/(2*r1*h3))*[(2*h2+b-c)*d1,(2*h2+b-c)*(a-d1); (2*h2+b+c)*d2, (2*h2+b+c)*(1-d2)];
            else
                ERR=54; % g2 out of bounds
                MAP={};
                return
            end
        end
    end
end
