function [cmap, re_cmap] = cl2parent(pre_cmap, remain_cmap, p_ids)
%generate a new color map for current frame with consistent color coding as
%the previous frame
%INPUT:
% pre_cmap: color map from previous frame
% remain_cmap: colors that not used in pre_camp
% p_ids: parent ids
%OUTPUT:
% cmap: new color map should be used in current frame
% re_cmap: color that not used

% contact: ccwang@vt.edu 02/25/2020

p_c = find(p_ids ~= 0);
op_c = find(p_ids == 0); % orphans

cmap = zeros(length(p_ids), 3);
cmap(p_c,:) = pre_cmap(p_ids(p_c), :);

pre_cmap(p_ids(p_c), :) = [];
pre_cmap = pre_cmap(randperm(length(pre_cmap)), :); % re-shuffle
remain_cmap = cat(1, remain_cmap, pre_cmap);
cmap(op_c,:) = remain_cmap(1:length(op_c),:);

re_cmap = remain_cmap(length(op_c)+1:end,:);
end