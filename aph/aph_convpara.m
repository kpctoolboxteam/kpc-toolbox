function [alpha,T] = aph_convpara(p)
    % perform convolution on parallel structure with any number of elements
    if length(p) == 2
        alpha = p{1};
        T = p{2};
    else
        [alpha,T] = aph_simplify(p{1},p{2},p{3},p{4},1,1,2);
        for i = 1:1:(length(p)/2)-2
            [alpha,T] = aph_simplify(alpha,T,p{3+2*i},p{4+2*i},1,1,1); 
        end
    end
end