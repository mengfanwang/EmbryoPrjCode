function [forward_mat, ratio_mat, track_head, ratio_mat2] = validate_res(movieInfo, gt)
% validate the tracking results
%INPUT:
% movieInfo: a structure saving all the essential information for tracking
% and from tracking. assume we have m traces
% gt: the ground truth locations of (part of )cells. assume we have n cells
%OUTPUT:
% forward_mat: a matrix with size n-by-m. each row represents a cell in gt.
% each element represents which trace have this detection
% ratio_mat: forward_mat divided by the length of tracks
% track_head: the track id which occupies the first element of the groundtruth
% trajectory

% ccwang@vt.edu, 03/02/2020

n_gt = numel(gt);
m_tr = numel(movieInfo.tracks);

forward_mat = nan(n_gt, m_tr);
track_head = nan(n_gt,1);
for i=1:n_gt
    if isempty(find(~isnan(gt{i}),1))
        continue;
    end
    trs = nan(length(gt{i}),1);
    trs(~isnan(gt{i})) = movieInfo.particle2track(gt{i}(~isnan(gt{i})),1);
    if ~isempty(find(~isnan(trs),1))
        track_head(i) = trs(find(~isnan(trs),1));
        uni_trs = unique(trs(~isnan(trs)));
        for j=1:length(uni_trs)
            forward_mat(i, uni_trs(j)) = sum(trs==uni_trs(j))/...
                (length(gt{i}) + length(movieInfo.tracks{uni_trs(j)})...
                - sum(trs==uni_trs(j)));
        end
    end
end

ratio_mat = forward_mat;%./cellfun(@length, gt);
ratio_mat2 = forward_mat;%./[cellfun(@length, movieInfo.tracks)]';
end