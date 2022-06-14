function [nei_segments, fg_locs] = gap_voxel_search(slabel, connect)
% given a 3d int label map, find the gap voxels between two different regions
%INPUT:
% slabel: the int label map
% fg_locs: the voxel indexes of foreground voxels
% connect: the way to determine neighbors
%OUTPUT:
% nei_segments: the two different regions that the voxels lay between

% contact: ccwang@vt.edu, 02/10/2020

if nargin<3
    connect = 124;
end
[h,w,z] = size(slabel);
label_di = zeros(h+6, w+6, z+6);
label_di(4:end-3, 4:end-3, 4:end-3) = slabel;

fg_locs = find(label_di>0);
[h,w,z] = size(label_di);
if connect == 26
    neighbors = [-h-1, -h, -h+1, -1, 1, h-1, h, h+1];
    neighbors = [neighbors-h*w, -h*w, neighbors, h*w, neighbors+h*w];
elseif connect == 124
    neighbors = [-2*h-2:-2*h+2, -h-2:-h+2, -2,-1,1,2, h-2:h+2, 2*h-2:2*h+2];
    neighbors = [neighbors-2*h*w, neighbors-h*w, -2*h*w, -h*w, neighbors...
        , 2*h*w, h*w, neighbors+h*w, neighbors+2*h*w];
else
    %connect == 49
    neighbors = [-3*h-3:-3*h+3, -2*h-3:-2*h+3, -h-3:-h+3, -3,-2,-1,1,2,3,...
        h-3:h+3, 2*h-3:2*h+3, 3*h-3:3*h+3];
    neighbors = [neighbors-3*h*w, neighbors-2*h*w, neighbors-h*w,...
        -3*h*w, -2*h*w, -h*w, neighbors, 3*h*w, 2*h*w, h*w, ...
        neighbors+h*w, neighbors+2*h*w, neighbors+3*h*w];
end

nei_mat = repmat(fg_locs, 1, length(neighbors)) +...
    repmat(neighbors, length(fg_locs), 1);

% tmp_mat = repmat(fg_locs, 1, length(neighbors));
% nei_mat(nei_mat<=0) = tmp_mat(nei_mat<=0);

fg_neighbors = label_di(nei_mat);
%     fg_neighbors(fg_neighbors==0) = nan;
%     find(nanmin(fg_neighbors, 2) < nanmean(fg_neighbors, 2));
%     neiCnt = sum(diff(sort(A,2),1,2)~=0,2)+1;
%     fg_neighbors = fg_neighbors();
nei_segments = nan(size(fg_neighbors,1),2);
for j=1:size(fg_neighbors,1)
    cur_neis = fg_neighbors(j,:);
    cur_neis(cur_neis == 0) = [];
    uni_neis = sort(unique(cur_neis));
    if length(uni_neis)==2
        nei_segments(j,:) = uni_neis;
    end
end

[cy,cx,cz] = ind2sub(size(label_di), fg_locs);

fg_locs = sub2ind(size(slabel), cy-3, cx-3, cz-3);
end