function MAP=map_mixture(alpha,MAPs)
% probabilistic mixture of MAPs - MAP=map_mixture(alpha,MAPs)
D0 = [];
D1 = [];
for i = 1:length(MAPs)
    if i==1        
        D0 = MAPs{i}{1};
    else
        D0 = blkdiag(D0,MAPs{i}{1});
    end
    D1i = alpha(1)*MAPs{i}{2}*ones(length(MAPs{i}{2}),1)*map_pie(MAPs{1});
    for j=2:length(MAPs)
        D1i = horzcat(D1i,MAPs{i}{2}*alpha(j)*ones(length(MAPs{i}{2}),1)*map_pie(MAPs{j}));
    end
    D1 = vertcat(D1,D1i);
end

MAP = map_normalize({D0,D1});

end