function [gt_tracks, gt_tracks_ov] = loc2id(gt_mat, detections)
% transfer the cell centers to corresponding detections
%INPUT:
% gt_mat: cells, each of which corresponds to a track, consisting of frame,
% y-loc, x-loc and z-loc
% detections: cells, each of which corresponds to a label map of a frame
%OUTPUT:
% gt_tracks: cells each of which corresponds to a track, consisting of the
% detection ids ordered by their original label and time information

% ccwang@vt.edu, 02/03/2020

gt_tracks = cell(numel(gt_mat), 1);
gt_tracks_ov = cell(numel(gt_mat), 1);
cells_fr = cellfun(@(x) max(x(:)), detections);
[hh,ww,zz] = size(detections{1});
used = [];
for i=1:numel(gt_mat)
    cur_tr = gt_mat{i};
    gt_tracks{i} = nan(size(cur_tr,1), 1);
    gt_tracks_ov{i} = nan(size(cur_tr,1), 1);
    for j=1:size(cur_tr,1)
        fr = cur_tr(j,1);
        coo = cur_tr(j,2:4);
        y = [floor(coo(1)) ceil(coo(1))];
        
        y(y<1 | y>hh) = [];
        x = [floor(coo(2)) ceil(coo(2))];
        x(x<1 | x>ww) = [];
        z = [floor(coo(3)) ceil(coo(3))];
        z(z<1 | z>zz) = [];
        coos = combvec(y, x, z);
        idx = sub2ind([hh,ww,zz],...
            coos(1,:)',coos(2,:)',coos(3,:)');
        ids = detections{fr}(idx);
        if isempty(ids(ids~=0))
            id = nan;
        else
            id = mode(ids(ids~=0));
        end
        id_ = id + sum(cells_fr(1:fr-1));
        if isempty(find(used==id_,1))
            gt_tracks{i}(j) = id_;
            gt_tracks_ov{i}(j) = id_;
            used(end+1) = id_;
        else
            gt_tracks{i}(j) = nan; % shared node, we will view as in-correct
            gt_tracks_ov{i}(j) = id_;
        end
            
    end
end
end
