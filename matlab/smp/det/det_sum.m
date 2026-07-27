function DET=det_sum(DETs)
% for i=2:length(DETs)
%     DETs{1}{1} = inv(inv(DETs{1}{1}) + inv(DETs{i}{1}));
%     DETs{1}{2} = DETs{1}{2} + DETs{i}{2};
% end
% DET=DETs{1};
DET=map_sum(DETs);
end