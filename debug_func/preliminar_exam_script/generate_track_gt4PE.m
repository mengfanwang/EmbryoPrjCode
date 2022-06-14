% build the gt of linkings
load('../data/crop_embryo_data_500x500x30x40/reference_track.mat', 'reference_track_mat');
for i=1:numel(reference_track_mat)
    track_info = reference_track_mat{i};
    track_info(:,2:3) = track_info(:,2:3)./2;
    reference_track_mat{i} = track_info;
end


it = numel(movieInfoAll);
movieInfo = movieInfoAll{it};%
refine_res = refine_resAll{it};
threshold_res = threshold_resAll{it};
[~, gt_tracks_ov] = loc2id(reference_track_mat, refine_res);
[~, ratio_mat, track_id_head, ~] = validate_res(movieInfo, gt_tracks_ov);
bestOvRatio = nanmax(ratio_mat,[],2);
disp(mean(bestOvRatio));

test_centers = cell(9,1);
for i=1:9
    test_track_id = track_id_head(i);
    test_centers{i} = display_track_link4GT(movieInfo, test_track_id,refine_res, ...
        embryo_vid, reference_track_mat{i});%embryo_vid);
    title(num2str(i));
end

datapath = 'C:\Users\Congchao\Desktop\cell_detection_samples\crop_embryo_data_500x500x30x40\images_downsample';
if ~isfolder(datapath)
    mkdir(datapath);
end
for i=1:40
    tifwrite(embryo_vid{i}/255, fullfile(datapath,num2str(i)));
end
[hh,ww,zz] = size(embryo_vid{1});
data =zeros(hh,zz, 40);
for i=1:40
    loc = find(yx_locs(:,1)==i);
    if isempty(loc)
        data(:,:,i) = squeeze(max(embryo_vid{i},[],2));
    else
        xx = floor(yx_locs(loc,3)-1): ceil(yx_locs(loc,3)+1);
        data(:,:,i) = squeeze(max(embryo_vid{i}(:,xx,:),[],2));
    end
end
zzshow(scale_image(data,0,1)*5);

data =zeros(ww,zz, 40);
for i=1:40
    loc = find(yx_locs(:,1)==i);
    if isempty(loc)
        data(:,:,i) = squeeze(max(embryo_vid{i},[],1));
    else
        yy = floor(yx_locs(loc,2)-1): ceil(yx_locs(loc,2)+1);
        data(:,:,i) = squeeze(max(embryo_vid{i}(yy,:,:),[],1));
    end
end
zzshow(scale_image(data,0,1)*5);



%% build the ground truth
