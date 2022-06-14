function display_multi_regions(movieInfo, refine_res, cur_id, bi_ids, display_flag)
% debug use: a temporary function to display infor of one region and its kid regions 

idx = movieInfo.voxIdx{cur_id};
if isempty(bi_ids)
    bi_ids = movieInfo.nei{cur_id};
end
idx_bi = cell(length(bi_ids), 1);
for i=1:length(bi_ids)
    idx_bi{i} = movieInfo.voxIdx{bi_ids(i)};
    fprintf('size of region %d: %d, region %d: %d, overlap: %d; ', cur_id, ...
        length(idx), bi_ids(i), length(idx_bi{i}), length(intersect(idx, idx_bi{i})));
    if bi_ids(i) > cur_id
        locs = find(movieInfo.nei{cur_id}==bi_ids(i));
        fprintf('CDist: %.3f, distance: %.3f\n', movieInfo.CDist{cur_id}(locs),...
            movieInfo.Cij{cur_id}(locs));
    else
        locs = find(movieInfo.nei{bi_ids(i)}==cur_id);
        fprintf('CDist: %.3f, distance: %.3f\n', movieInfo.CDist{bi_ids(i)}(locs),...
            movieInfo.Cij{bi_ids(i)}(locs));
    end
end
if nargin<5
    display_flag = false;
end
if display_flag
    zzmap = zeros(size(refine_res{1}));
    zzmap(idx) = 1;
    zzmap(cat(1, idx_bi{:})) = 2 + zzmap(cat(1, idx_bi{:}));
    zzshow(label2rgb3d(zzmap));
end
end