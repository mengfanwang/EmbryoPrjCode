function out_idsAndLocs = extractCellLocsInMovieInfo(idMap, movieInfo, frame)
% extract the ids in moveiInfo of the real ids in idMap

real_ids = unique(idMap(idMap>0));
if isempty(real_ids)
    out_idsAndLocs = [];
    return;
end
out_locs = real_ids;
if frame > 1
    old_ids = out_locs<=movieInfo.n_perframe(frame);
    out_locs(old_ids) = sum(movieInfo.n_perframe(1:frame-1)) + ...
        out_locs(old_ids);
end
out_idsAndLocs = [real_ids out_locs];
end