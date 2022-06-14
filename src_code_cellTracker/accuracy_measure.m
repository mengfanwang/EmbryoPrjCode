function track_id_head = accuracy_measure(movieInfo, refine_res, gt_tracks, metric)
% pick a way to measure the accuracy of current model
if nargin == 3
    metric = 2;
end
if metric == 1 % way one: iou
    [~, gt_tracks_ov] = loc2id(gt_tracks, refine_res);
    [~, ratio_mat, track_id_head, ~] = validate_res(movieInfo, gt_tracks_ov);
    bestOvRatio = nanmax(ratio_mat,[],2);
    disp(nanmean(bestOvRatio));
elseif metric == 2 % recall
    accuracy = validate_link_acc(movieInfo, [], gt_tracks);
    accuracy(cellfun(@isempty, gt_tracks),:) = [];
    disp([sum(accuracy(:,1))./ sum(accuracy(:,2)) mean(accuracy(:,3))])
end
end