function centers = display_track_link4GT_MIBased(movieInfo, track_id, ...
    refine_res, embryo_vid, ref_locs)
% debug use: display the linking condition of a track
%

cur_track = movieInfo.tracks{track_id};
if length(cur_track) == 1
    fprintf('track length is 1 with node %d\n', cur_track);
    return;
end
blue_label_scale = 2;
if nargin > 3 && ~isempty(refine_res) && ~isempty(embryo_vid)
    [h,w,zslice] = size(refine_res{1});
    frs = movieInfo.frames(cur_track);
    unifrs = unique(frs);
    segMaxProj = uint8(zeros(h,w,3,numel(refine_res)));
    for i=1:numel(refine_res)
        tmp_g = max(embryo_vid{i},[],3);
        segMaxProj(:,:,2,i) = tmp_g*3;
    end
    centers = zeros(length(unifrs), 4);
    centers(:,1) = unifrs;
    for i=1:length(unifrs)
        voxIdx = cat(1, movieInfo.voxIdx{cur_track(frs==unifrs(i))});
        [yy,xx,zz] = ind2sub([h,w,zslice], voxIdx);
        centers(i,2:4) = mean([yy,xx,zz]);
        tmp_fr = false(h,w);
        tmp_fr(sub2ind([h,w], yy, xx)) = true;
        tmp_g = max(embryo_vid{unifrs(i)},[],3);
        %segMaxProj(:,:,2,unifrs(i)) = tmp_g;
        tmp_b = uint8(zeros(h,w));
        tmp_b(tmp_fr) = tmp_g(tmp_fr)*blue_label_scale;
        segMaxProj(:,:,3,unifrs(i)) = tmp_b;
    end
    
    % add ref_locs
    for i=1:size(ref_locs,1)
        tmp_locs = ref_locs(i,:);
        tmp_locs(:,4) = 1;% max project no need z loc
        idx = get_index_from_center(tmp_locs(2:4), embryo_vid{1});
        cur_f = segMaxProj(:,:,:,tmp_locs(1));
        cur_f(idx) = 255;
        segMaxProj(:,:,:,tmp_locs(1)) = cur_f;
    end
    zzshow(segMaxProj);
end