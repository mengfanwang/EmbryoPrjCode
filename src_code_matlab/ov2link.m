function [parents, ious] = ov2link(pre_idmap, cur_idmap, cell_ids)
%find the parent of current cell in previous frame
% INPUT:
% pre_idmap: id map of the previous frame
% cur_idmap: id map of the current frame
% cell_ids: the cells to be considered
%OUTPUT?
% parents: the parent cells' id
% ious: their overlapping ratio between current cell and its parent

% contact: ccwang@vt.edu, 02/25/2020

if abs(numel(pre_idmap)-numel(cur_idmap))>0
    parents = cell_ids * 0;
    ious = parents;
    return;
end
parents = zeros(length(cell_ids),1);
ious = zeros(length(cell_ids),1);
cur_voxlist = regionprops3(cur_idmap, 'VoxelIdxList');
pre_voxlist = regionprops3(pre_idmap,'VoxelIdxList');
for i=1:length(cell_ids)
    voxs = cur_voxlist.VoxelIdxList{cell_ids(i)};
    ids = pre_idmap(voxs);
    if ~isempty(ids(ids~=0))
        parent = mode(ids(ids~=0));
        parents(i) = parent;
        voxs_p = pre_voxlist.VoxelIdxList{parent};
        u = length(voxs) + length(voxs_p) - sum(ids==parent);
        ious(i) = sum(ids==parent) / u;
    end
end
end