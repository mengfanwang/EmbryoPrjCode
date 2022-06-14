function res = validateMovieInfo(gt_mat, out_refine_res, movieInfo)

res = cell(5,1);
cnt = 0;
[gt_tracks, gt_tracks_ov] = loc2id(gt_mat, out_refine_res);
[~, ratio_mat, ~, ratio_mat2] = validate_res(movieInfo, gt_tracks);
zz = round(nanmax(ratio_mat,[],2),2);
zz(isnan(zz)) = 0;
zz1 = round(nanmax(ratio_mat2,[],2),2);
zz1(isnan(zz1)) = 0;
cnt = cnt + 1;
res{cnt,1} = [zz, zz1];
% if consider overlapped nodes: e.g. two gt_tracks overlapped with one
% detected track, they are both viewed as correct
[~, ratio_mat, ~, ratio_mat2] = validate_res(movieInfo, gt_tracks_ov);
zz = round(nanmax(ratio_mat,[],2),2);
zz(isnan(zz)) = 0;
zz1 = round(nanmax(ratio_mat2,[],2),2);
zz1(isnan(zz1)) = 0;
res{cnt,2} = [zz, zz1];

fprintf('If consider no overlap, normorlized by gt and tracking: %.2f, %.2f; %d, %d\n',...
    mean(res{cnt,1}), sum(res{cnt,1}==1));
fprintf('If consider overlap, normorlized by gt and tracking: %.2f, %.2f; %d, %d\n',...
    mean(res{cnt,2}), sum(res{cnt,2}==1));

end