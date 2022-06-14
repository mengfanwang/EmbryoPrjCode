function ov_reg_id = search_tail_from_center(pt_xyz, f, movieInfo, ...
    refine_res, g)
% given a point, find the overlapped regions in previous several frames


idx = get_index_from_center(pt_xyz([2 1 3]), size(refine_res{1}));

ov_reg_id = [];
for i=1:g.k
    f = f-1;
    if f<1
        break;
    end
    ids = refine_res{f}(idx);
    
    if f > 1
        tmp = unique(ids(ids>0)) + sum(movieInfo.n_perframe(1:f-1));
    else
        tmp = unique(ids(ids>0));
    end
    if ~isempty(tmp)
        valid_ones = ~isnan(movieInfo.particle2track(tmp)) & ...
            cellfun(@isempty, movieInfo.kids(tmp));
        ov_reg_id = cat(1, ov_reg_id, tmp(valid_ones));
    end
end

end