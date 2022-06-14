function validate_TGMM_res(embryo_vid, sc_f)
addpath('./TGMM_wrapper/');
% use cropped data to validate the results of TGMM
load('../data/crop_embryo_data_500x500x30x40/embryo_image_detection.mat','dets');

%sc_f = 2; % scaling factor, we resize the data to [h/sc_f, w/sc_f, z, t]

for i=1:numel(dets)
    dets{i}(:,1:2) = dets{i}(:,1:2)/sc_f;
end

%% build movieInfo_TGMM
det_maps = cell(numel(dets),1);
for i=1:numel(det_maps)
    det_maps{i} = zeros(size(embryo_vid{i}));
    for j=1:size(dets{i},1)
        f_pt = max(floor(dets{i}(j,1:3)), 1);
        c_pt = min(ceil(dets{i}(j,1:3)), size(embryo_vid{i}));
%         f_pt = max(floor(dets{i}(j,1:3)), 1);
%         c_pt = min(ceil(dets{i}(j,1:3)), size(embryo_vid{i}));
        [yy, xx, zz] = meshgrid(f_pt(1):c_pt(1), f_pt(2):c_pt(2), f_pt(3):c_pt(3));
        tmp_idx = sub2ind(size(det_maps{i}), yy, xx, zz);
        det_maps{i}(tmp_idx) = j;
    end
end
movieInfo_TGMM = tree2tracks_cell(det_maps, [], dets);

[~, gt_tracks_ov] = loc2id(reference_track_mat, det_maps);
[~, ratio_mat, track_id_head, ~] = validate_res(movieInfo_TGMM, gt_tracks_ov);
bestOvRatio = nanmax(ratio_mat,[],2);

%%
p = track_disp_para();
%if ~exist('dataWithParticles', 'var')
movieInfo_TGMM.orgCoord = [movieInfo_TGMM.xCoord, movieInfo_TGMM.yCoord, movieInfo_TGMM.zCoord];
dataWithParticles = generate3DImWithParticle(...
    cat(4, embryo_vid{:}), movieInfo_TGMM, p);
%end
write_id = unique(track_id_head);
%[5 31 47 102 258 254 162];%find(res{cnt,1}<0.8); % consider no overlap
drawTracksStart1st(cat(4, embryo_vid{:}), fullfile(res_folder,'TGMM_res'),...
    movieInfo_TGMM, dataWithParticles, write_id);

%% 
for i=11:20%:20%[18 30]%[12 28 30 31 32] %
    test_track_id = track_id_head(i); % -8 is ok, -7 broken, -6,-5...are disaster
    display_track_link(movieInfo_TGMM, test_track_id,det_maps, embryo_vid);
end
% check more tracks: from brighter to dimmer
leftover_ids = setdiff(1:numel(movieInfo_TGMM.tracks), track_id_head);
brightness = leftover_ids * 0;
for i=1:length(leftover_ids)
    if ~isempty(movieInfo_TGMM.tracks{leftover_ids(i)})
        head_id = movieInfo_TGMM.tracks{leftover_ids(i)}(1);
        brightness(i) = mean(embryo_vid{movieInfo_TGMM.frames(head_id)}(movieInfo_TGMM.voxIdx{head_id}));
    end
end
[~, od] = sort(brightness, 'descend');
leftover_ids = leftover_ids(od);
for i=31:50%:20%[18 30]%[12 28 30 31 32] %
    test_track_id = leftover_ids(i);
    display_track_link(movieInfo_TGMM, test_track_id,det_maps, embryo_vid);
end
tested = cat(2, track_id_head', voxes(1:50))';
end