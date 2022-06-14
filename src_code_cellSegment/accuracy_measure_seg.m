function [iou, recall, best_match] = accuracy_measure_seg(gt, id_map, display_flag)

if nargin <3
    display_flag = false;
end
num_cell = max(gt(:));
iou = zeros(num_cell, 1);
best_match = nan(num_cell, 1);
gt_voxIdx = regionprops3(gt, 'VoxelIdxList');
gt_voxIdx = gt_voxIdx.VoxelIdxList;
test_voxIdx = regionprops3(id_map, 'VoxelIdxList');
test_voxIdx = test_voxIdx.VoxelIdxList;
for i=1:num_cell
    cur_voxidx = gt_voxIdx{i};
    if isempty(cur_voxidx)
        iou(i) = nan;
        continue;
    end
    ids = unique(id_map(cur_voxidx));
    ids = ids(ids>0);
    for j=1:length(ids)
        tmp_iou = length(intersect(test_voxIdx{ids(j)}, cur_voxidx)) / ...
            length(union(test_voxIdx{ids(j)}, cur_voxidx));
        if tmp_iou > iou(i)
            iou(i) = tmp_iou;
            best_match(i) = ids(j);
        end
    end
end
iou_thres = 0.1:0.05:0.8;
missing_rate = iou_thres*0;
for i=1:length(iou_thres)
    missing_rate(i) = sum(iou<iou_thres(i)) / length(iou);
end
recall = 1-missing_rate;
if display_flag
   figure; plot(iou_thres, 1-missing_rate);
   xlabel('iou');ylabel('recall');
end