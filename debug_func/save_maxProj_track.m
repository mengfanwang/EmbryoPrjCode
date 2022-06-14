function saved = save_maxProj_track(movieInfo, track_id, ...
    refine_res, embryo_vid, save_path, save_track_id)
% debug use: display the linking condition of a track
%
if nargin == 5
    save_track_id = track_id;
end
saved = false;
cur_track = movieInfo.tracks{track_id};
if length(cur_track) == 1
    fprintf('track length is 1 with node %d\n', cur_track);
    return;
end
blue_label_scale = 3;
[h,w,zslice] = size(refine_res{1});
frs = movieInfo.frames(cur_track);
unifrs = unique(frs);
segMaxProj = uint8(zeros(h,w,3,numel(refine_res)));
for i=1:numel(refine_res)
    tmp_g = max(embryo_vid{i},[],3);
    segMaxProj(:,:,2,i) = tmp_g;
end
cur_t_y = zeros(length(unifrs),2);
cur_t_y(:,1) = unifrs;
temporal_change = zeros(length(unifrs), 2);
for i=1:length(unifrs)
    voxIdx = cat(1, movieInfo.voxIdx{cur_track(frs==unifrs(i))});
    temporal_change(i,1) = mean(embryo_vid{unifrs(i)}(voxIdx));
    temporal_change(i,2) = length(voxIdx);
    cur_t_y(i,2) = length(voxIdx);
    [yy,xx,~] = ind2sub([h,w,zslice], voxIdx);
    tmp_fr = false(h,w);
    tmp_fr(sub2ind([h,w], yy, xx)) = true;
    tmp_g = max(embryo_vid{unifrs(i)},[],3);
    %segMaxProj(:,:,2,unifrs(i)) = tmp_g;
    tmp_b = uint8(zeros(h,w));
    tmp_b(tmp_fr) = tmp_g(tmp_fr)*blue_label_scale;
    segMaxProj(:,:,3,unifrs(i)) = tmp_b;
end

tifwrite(segMaxProj, fullfile(save_path, ['max_proj_track_',...
    num2str(save_track_id)]));
saved = true;