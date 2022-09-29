function [newIdMap, cnt, map] = rearrange_id_oldversion(idMap)
% re-arrange the id in idMap: this one is slow
% INPUT:
% idMap: the original id map, the ids may not be consecutive
% OUTPUT:
% newIdMap: new id map, the ids are consecutive now

% ccwang@vt.edu, 01/22/2020

max_id = max(idMap(:));
if max_id == 0
    newIdMap = idMap;
    cnt = 0;
    map = [];
    return;
end
if max_id < 1e4 
    newIdMap = zeros(size(idMap));
    s = regionprops3(idMap, {'VoxelIdxList'});
    cnt = 0;
    
    if ~iscell(s.VoxelIdxList)
        cnt = cnt + 1;
        newIdMap(s.VoxelIdxList) = cnt;
        map = [idMap(s.VoxelIdxList(1)) 1];
    else
        map = zeros(numel(s.VoxelIdxList), 2);
        fprintf('Rearrange old version:');
        for i=1:numel(s.VoxelIdxList)
            fprintf('%d ', i);
            if ~isempty(s.VoxelIdxList{i})
                cnt = cnt + 1;
                newIdMap(s.VoxelIdxList{i}) = cnt;
                map(cnt,:) = [i cnt];
            end
        end
        map = map(1:cnt, :);
    end
else % need to change to avoid bug
    error('Too much seeds!')
    l_map = bwlabeln(idMap>0, 26);
    newIdMap = zeros(size(l_map));
    s = regionprops3(l_map, {'VoxelIdxList'});
    cnt = 0;
    map = zeros(numel(s.VoxelIdxList)*10, 2);
    for i=1:numel(s.VoxelIdxList)
        cur_ls = idMap(s.VoxelIdxList{i});
        ids = unique(cur_ls);
        if length(ids) >= 1
            for j = 1:length(ids)
                cnt = cnt + 1;
                newIdMap(s.VoxelIdxList{i}(cur_ls==ids(j))) ...
                    = cnt;
                map(cnt,:) = [ids(j) cnt];
            end
        end
    end
    map = map(1:cnt, :);
end

end