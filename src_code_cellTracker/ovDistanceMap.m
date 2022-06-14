function [neighbors, distances, distances_dirwise] = ovDistanceMap(curMap,...
    nextMap, curVoxIdxCells, nextVoxIdxCells, drift)
% calculate the distance between regions across frames based on overlapping ratio
%INPUT:
% curMap: current label map consisting all the current regions
% nextMap: next label map consisting all neighbor regions
% curVoxIdxCells: the voxel idxs of all regions in curMap --cells
% nextVoxIdxCells: the voxel idxs of all regions in nextMap --cells
%OUTPUT:
% neighbors: cells containing the neighbor id
% distances: cells overlapping distances between current region and its
% corresponding neighbors

% contact: ccwang@vt.edu, 02/26/2020
if nargin == 4
    drift = nan;
end
neighbors = cell(numel(curVoxIdxCells),1);
distances = cell(numel(curVoxIdxCells),1);
distances_dirwise = cell(numel(curVoxIdxCells),1);
[h,w,z] = size(curMap);
r = zeros(1e7,1);
r_cnt = 0;
for cur_i = 1:numel(curVoxIdxCells)
%     disp(cur_i);
    curVoxIdx = curVoxIdxCells{cur_i};
    [curVox_y, curVox_x, curVox_z] = ind2sub([h,w,z], curVoxIdx);
    
    cur_Vox_yxz = [curVox_y, curVox_x, curVox_z];
    if isnan(drift)
        ids = nextMap(curVoxIdx);
    else
        shifted_vox = cur_Vox_yxz + round(drift([2, 1, 3]));
        shifted_vox = verifyInData(shifted_vox, size(nextMap));
        shifted_voxIdx = sub2ind(size(nextMap), shifted_vox(:, 1),...
            shifted_vox(:, 2), shifted_vox(:, 3));
        ids = nextMap(shifted_voxIdx);
    end
    neighbors{cur_i} = double(unique(ids(ids~=0)));
    distances{cur_i} = nan(length(neighbors{cur_i}), 1);
    distances_dirwise{cur_i} = nan(length(neighbors{cur_i}), 2); % i2j and j2i
    
%     if length(curVoxIdx) > 500
%         rd_ids = randperm(length(curVoxIdx));
%         curVoxIdx = curVoxIdx(rd_ids(1:500));
%     end
    for j=1:length(neighbors{cur_i})
        %disp(j);
        neiId = neighbors{cur_i}(j);
        nextVoxIdx = nextVoxIdxCells{neiId};
%         if length(nextVoxIdx) > 500
%             rd_ids = randperm(length(nextVoxIdx));
%             nextVoxIdx = nextVoxIdx(rd_ids(1:500));
%         end
        r_cnt = r_cnt + 1;
        [nextVox_y, nextVox_x, nextVox_z] = ind2sub([h,w,z], nextVoxIdx);
        if isnan(drift)
            [distances{cur_i}(j), ~, distances_dirwise{cur_i}(j,:), r(r_cnt)] = ...
                ovDistanceRegion(cur_Vox_yxz, [nextVox_y, nextVox_x, nextVox_z], [0 0 0]);
        else
            [distances{cur_i}(j), ~, distances_dirwise{cur_i}(j,:), r(r_cnt)] = ...
                ovDistanceRegion(cur_Vox_yxz, [nextVox_y, nextVox_x, nextVox_z], drift([2, 1, 3]));
        end
        
    end
end
end