function [alpha,T] = aph_convseq(s)
    % perform convolution on sequential with any number of elements
    if length(s) == 2
        alpha = s{1};
        T = s{2};
    else
        [alpha,T] = aph_simplify(s{1},s{2},s{3},s{4},1,1,1); 
        for i = 1:1:(length(s)/2)-2
            [alpha,T] = aph_simplify(alpha,T,s{3+2*i},s{4+2*i},1,1,1); 
        end
    end
end
