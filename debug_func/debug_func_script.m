%% debug:
q = initial_q(2, true);

%% debug add missing cell
it = numel(movieInfoAll)-1;
movieInfo = movieInfoAll{it};
refine_res = refine_resAll{it};
threshold_res = threshold_resAll{it};

profile on;
[movieInfo, refine_res, thresholdMaps, added_cell_id] = ...
    add_missing_cell(movieInfo, refine_res, embryo_vid, g, ...
    threshold_res, eigMaps, varMaps, q);
profile viewer;
profile off
[movieInfo1, refine_res1, g1] = movieInfo_update(movieInfo, ...
    refine_res, added_cell_id, g);

%% debug regionRefresh
it = numel(movieInfoAll);

movieInfo = movieInfoAll{it};%
refine_res = refine_resAll{it};
threshold_res = threshold_resAll{it};
[~, gt_tracks_ov] = loc2id(g.gt_mat, refine_res);
[~, ratio_mat, track_id_head, ~] = validate_res(movieInfo, gt_tracks_ov);
max_ov_vec = nanmax(ratio_mat,[],2);
disp(nanmean(max_ov_vec));

accuracy = validate_link_acc(movieInfo, [], g.gt_mat);
accuracy(cellfun(@isempty, g.gt_mat),:) = [];
disp(sum(accuracy(:,1))./ sum(accuracy(:,2)));
disp([mean(accuracy(:,3)), std(accuracy(:,3))])

regionRefresh(movieInfo, refine_res, embryo_vid, g, eigMaps, q);
rule_based_division_detect(movieInfo, refine_res, embryo_vid, g);
if exist('vidMap','var')
    embryo_vid = vidMap;
else
    vidMap = embryo_vid;
end
%test_voxIdx = movieInfo.voxIdx{3237};
test_seg_id = sum(movieInfo.n_perframe(1:28))+...
    unique(refine_res{29}(test_voxIdx));
test_track_id = movieInfo.particle2track(2831,1);
test_track_id = track_id_head(4);
zz = display_track_link(movieInfo, test_track_id,refine_res, embryo_vid);
costs = arc_cost_in_track(movieInfo, test_track_id);

idsAtFrame = movieInfo.tracks{test_track_id}...
    (movieInfo.frames(movieInfo.tracks{test_track_id})==10);

test_seg_id = 1490;
zz = display_seg_res(movieInfo, refine_res, embryo_vid,  test_seg_id, threshold_res);
zz = display_seg_res(movieInfo, refine_res, embryo_vid,  test_seg_id);

tifwrite(zz, fullfile(ppt_dir, 'b'));
figure;imshow(max(zz,[],4))
display_multi_regions(movieInfo, refine_res, test_seg_id,[]);
display_track_from_root(movieInfo, test_seg_id,refine_res, embryo_vid);

principla_3d = true;
display_gap_res(movieInfo, refine_res, embryo_vid, eigMaps, test_seg_id,principla_3d);
principla_3d = false;
display_gap_res(movieInfo, refine_res, embryo_vid, eigMaps, test_seg_id,principla_3d);

for i=1:18
    display_track_link(movieInfo, track_id_head(i),refine_res,embryo_vid);
    title(num2str(i));
end

for i=[7 28 38 43]%[5 12 14 18 19 22 26 28 30 32 36 43]
    display_track_link(movieInfo, track_id_head(i),refine_res,embryo_vid);
    title(num2str(i));
end

c1 = [1843 1844];%[1402 1403];%2653;%
c2 = 1964;
[mc1, mc2] = voxIdx2cost(cat(1, movieInfo.voxIdx{c1}), ...
    cat(1, movieInfo.voxIdx{c2}), ...
    [movieInfo.frames(c1(1)) movieInfo.frames(c2(1))], ...
    movieInfo, size(refine_res{1}), movieInfo.jumpRatio)

%% debug handleMergeSplitRegions.m
movieInfo = handleMergeSplitRegions(movieInfo, g);
%% vox cost
c1 = [5784 5785];%2653;%
c2 = 5591;
[mc1, mc2] = voxIdx2cost(cat(1, movieInfo.voxIdx{c1}), ...
    cat(1, movieInfo.voxIdx{c2}), ...
    [movieInfo.frames(c1(1)) movieInfo.frames(c2(1))], ...
    movieInfo, size(refine_res{1}), movieInfo.jumpRatio)
movieInfo.preNei{2434}
%% debug drift correction
movieInfo = driftFromTracks(movieInfo,g);
%% debug separateRegion.m
%sepReg = {[960 1189 1203]};
sepReg = {[1225 1095 1203]};
separateRegion(sepReg, embryo_vid, refine_res, movieInfo, eigMaps, g, q);


%% others
total_detection = 0;
for i=1:21
    total_detection = total_detection + ...
        range(movieInfo.frames(movieInfo.tracks{track_id_head(i)}))+1;
end

%08/02 ==> [18->cell split(noisy cell show) 12&31 30->data_crop(cell on boundary) 28->cell split]
% 
for i=1:100%:20%[18 30]%[12 28 30 31 32] %
    test_track_id = track_id_head(i); % -8 is ok, -7 broken, -6,-5...are disaster
    display_track_link(movieInfo, test_track_id,refine_res, embryo_vid);
end
% check more tracks: from brighter to dimmer
leftover_ids = setdiff(1:numel(movieInfo.tracks), track_id_head);
brightness = leftover_ids * 0;
for i=1:length(leftover_ids)
    if ~isempty(movieInfo.tracks{leftover_ids(i)})
        head_id = movieInfo.tracks{leftover_ids(i)}(1);
        brightness(i) = mean(embryo_vid{movieInfo.frames(head_id)}(movieInfo.voxIdx{head_id}));
    end
end
[~, od] = sort(brightness, 'descend');
leftover_ids = leftover_ids(od);
for i=31:50%:20%[18 30]%[12 28 30 31 32] %
    test_track_id = leftover_ids(i);
    display_track_link(movieInfo, test_track_id,refine_res, embryo_vid);
end

%% save the 106 trajectories for further check
intersted_track = track_id_head;
intersted_track([54 58 84 86 89 90 91 92 97 98]) = [];% remove noisy tracks
more_valid_track_ids = [1:21, 23:25, 27:30, 32, 34 ,40:42, 47];
intersted_track = unique(cat(1, intersted_track, leftover_ids(more_valid_track_ids)'));

g_weak = [105,61,247,40,237,10,20,194, 102, 162];
b_undecidable = [266 307 236];
b_multipleCell = [5 7 222 22 181 315];
g_ones = setdiff(intersted_track, cat(2, g_weak, b_undecidable, b_multipleCell));
g_weak(end) = [];% too weak indeed

save_ids = cat(2, b_multipleCell, b_undecidable, g_ones', g_weak);
reference_track_mat = cell(length(save_ids), 1);
for i=1:length(save_ids)
    cur_track = movieInfo.tracks{save_ids(i)};
    track_info = zeros(length(cur_track),4);
    for j=1:length(cur_track)
        [yy, xx, zz] = ind2sub(size(refine_res{1}), movieInfo.voxIdx{cur_track(j)});
        track_info(j,1) = movieInfo.frames(cur_track(j));
        track_info(j,2:3) = mean([yy, xx]);% 
        track_info(j,4) = mean(zz);
    end
    reference_track_mat{i} = track_info;
end
%save('../data/crop_embryo_data_500x500x30x40/reference_track.mat', 'reference_track_mat');


%% test reference tracks
load('../data/crop_embryo_data_500x500x30x40/reference_track.mat', 'reference_track_mat');
for i=1:numel(reference_track_mat)
    track_info = reference_track_mat{i};
    track_info(:,2:3) = track_info(:,2:3)./2;
    reference_track_mat{i} = track_info;
end

it = numel(movieInfoAll)-1;
movieInfo = movieInfoAll{it};%
refine_res = refine_resAll{it};
threshold_res = threshold_resAll{it};
[~, gt_tracks_ov] = loc2id(g.gt_mat, refine_res);%reference_track_mat, refine_res);
[~, ratio_mat, track_id_head, ~] = validate_res(movieInfo, gt_tracks_ov);
bestOvRatio = nanmax(ratio_mat,[],2);
disp(mean(bestOvRatio));

for i=[5 12 14 18 19 22 26 28 30 32 36 43]
    display_track_link(movieInfo, track_id_head(i),refine_res,embryo_vid);
    title(num2str(i));
end
display_track_link(movieInfo, 257,refine_res, embryo_vid);
display_track_link(movieInfo, [2 29],refine_res, embryo_vid);
display_rgb_time_course(refine_res)
display_multi_regions(movieInfo, refine_res, 2373,[]);
display_seg_res(movieInfo, refine_res, embryo_vid,  2527);
for i=1:6
    it = i;%numel(movieInfoAll)-1;
    movieInfo = movieInfoAll{it};%
    refine_res = refine_resAll{it};
    threshold_res = threshold_resAll{it};
    [~, gt_tracks_ov] = loc2id(reference_track_mat, refine_res);
    [~, ratio_mat, track_id_head, ~] = validate_res(movieInfo, gt_tracks_ov);
    bestOvRatio = nanmax(ratio_mat,[],2);
    %disp(mean(bestOvRatio));

    test_track_id = track_id_head(1);
    zzMax = display_track_link4PE(movieInfo, test_track_id,refine_res, ...
        embryo_vid);
    %tifwrite(zzMax, 'illustration_split_merge_3');
end

%% save tracks for display
for i= [3 8 22, 23, 25, 32, 34, 42 58]
    %further test needed--splitting traces: 3, 8, 22, 23, 25, 32, 34, 42, 58...
    % 57(splitted (s)), 56(s), 66(s), 76(s), 75(s), 74(s), 73, 86(s),
    % 85(s), 82(s), 92(s), 91(s);
    display_track_link4GT(embryo_vid, gt_tracks{i});title(num2str(i));
end


for i= [57 56 66 76 75 74 73 86 85 82 92 91]
    %further test needed--splitting traces: 3, 8, 22, 23, 25, 32, 34, 42, 58...
    % 57(splitted (s)), 56(s), 66(s), 76(s), 75(s), 74(s), 73, 86(s),
    % 85(s), 82(s), 92(s), 91(s);
    display_track_link4GT(embryo_vid, gt_tracks{i});title(num2str(i));
end


rr = [57 56 66 76 75 74 73 86 85 82 92 91, 61 , 83, 104];
gt_tracks(rr) = {[]};

gt_track_split_parent = gt_tracks([3 8 22, 23, 25, 32, 34, 42 58]);
gt_tracks([3 8 22, 23, 25, 32, 34, 42 58]) = {[]};
gt_track_split_kids = cell(numel(gt_track_split_parent),2);
g_kid_single_column = gt_track_split_kids';
g_kid_single_column = g_kid_single_column(:);

gt_tracks_no_split = gt_tracks(~cellfun(@isempty,gt_tracks));
gt_tracks_all = cat(1, gt_tracks_no_split, ...
    gt_track_split_parent, g_kid_single_column);

% save('groundtruth_tracks_tyxz_with_split.mat',...
%     'gt_tracks_all', 'gt_tracks_no_split', 'gt_track_split_parent',...
%     'gt_track_split_kids');




