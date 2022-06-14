function display_track_link4GT(embryo_vid, cur_track_locs)
% cur_track: frame, y, x, z
%

%blue_label_scale = 4;
[h,w,zslice] = size(embryo_vid{1});
frs = cur_track_locs(:,1);
unifrs = unique(frs);
segMaxProj = uint8(zeros(h,w,3,numel(embryo_vid)));
for i=1:numel(embryo_vid)
    tmp_g = max(embryo_vid{i},[],3);
    segMaxProj(:,:,2,i) = tmp_g*5;
end

for i=1:length(unifrs)
    cur_centers = cur_track_locs(frs==unifrs(i), 2:4);
    cur_centers(:,end) = 1; %max project, so z is set to one
    for j=1:size(cur_centers,1)
        idx = get_index_from_center(cur_centers(j,:), embryo_vid{1}, [1 1 0]);

        tmp_b = uint8(zeros(h,w));
        tmp_b(idx) = 255;
        segMaxProj(:,:,3,unifrs(i)) = tmp_b;
    end
end

zzshow(segMaxProj);
