function division_candidate_ids = extract_division_candidate(movieInfo, refine_res, ...
    embryo_vid, g)
% For a cell, if it 
% 1) happens suddenly in the field of view
% 2) does not touch the boundary of the field of view
% then it is highly likely a cell from cell division, we call them candidate.
% 3) more criterion???

trace_heads = cellfun(@(x) x(1), movieInfo.tracks);

trace_heads = trace_heads(movieInfo.frames(trace_heads) > 1); % rule 1)
sz_fov = size(refine_res{1});
for i=1:length(trace_heads)
    xyz = movieInfo.vox{trace_heads(i)};
    leftup = min(xyz,[], 1);
    leftup = leftup([2 1 3]);
    rigthbottom = max(xyz,[], 1);
    rigthbottom = rigthbottom([2 1 3]);
    if ~isempty(find(leftup == 1, 1)) || ...
            ~isempty(find(rigthbottom == sz_fov, 1)) % rule 2): boundary
        trace_heads(i) = nan;
    end
end


division_candidate_ids = trace_heads(~isnan(trace_heads));