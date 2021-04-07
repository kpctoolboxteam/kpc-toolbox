function [alpha,T] = aph_simplify(a1,T1,a2,T2,p1,p2,pattern)
    % convolve two ME distributions with parameters a1,T1,a2,T2
    % The input variable pattern switch this function to different pattern convolution
    % p1 and p2 are the specified branch probabilities for the third pattern (branch)
    if pattern == 1 % sequence structure
        Size1 = size(a1); 
        order1 = Size1(2); 
        Size2 = size(a2); 
        order2 = Size2(2);
        e1 = ones(order1,1); 
        alpha = [a1,(1-a1*e1)*a2];
        T = [T1,(-T1*e1)*a2;zeros(order2,order1),T2];

    elseif pattern == 2 % parallel strucutre
        Size1 = size(a1); 
        order1 = Size1(2); 
        Size2 = size(a2); 
        order2 = Size2(2); 
        e1 = ones(order1,1); 
        e2 = ones(order2,1); 
        alpha = [kron(a1,a2),(1-a2*e2)*a1,(1-a1*e1)*a2];
        Tr1 = [kron(T1,eye(order2))+kron(eye(order1),T2),kron(eye(order1),-T2*e2),kron(-T1*e1,eye(order2))];
        Tr2 = [zeros(order1,order1*order2),T1,zeros(order1,order2)];
        Tr3 = [zeros(order2,order1*order2),zeros(order2,order1),T2];
        T = [Tr1;Tr2;Tr3];

    elseif pattern == 3 % branch structure
        Size1 = size(a1); 
        order1 = Size1(2); 
        Size2 = size(a2); 
        order2 = Size2(2);
        alpha = [p1*a1,p2*a2];
        T = [T1,zeros(order1,order2);zeros(order2,order1),T2];

    else % loop structure
        %alpha = a1;
        %T = [T1(1,1),T1(1,2);-a1(1)*p1*T1(2,2),T1(2,2)-T1(2,2)*a1(2)*p1];
    end
end