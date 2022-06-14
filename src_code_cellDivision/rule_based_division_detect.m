function movieInfo = rule_based_division_detect(movieInfo, refine_res, embryo_vid, g)
% For a cell, if it 
% 1) happens suddenly in the field of view
% 2) does not touch the boundary of the field of view
% then it is highly likely a cell from cell division, we call them candidate.

% We search the sourrouding area of the cell to see if a peer cell exist.
% Such peer cell should be either a candidate too or has a meaningful
% parent.
% We move these two cells together to their center, they jointly should be
% good region from their parent.

% ccwang@vt.edu


division_candidate_ids = extract_division_candidate(movieInfo, refine_res, ...
    embryo_vid, g);

for i=1:length(division_candidate_ids)
    if isnan(division_candidate_ids(i))
        continue;
    end
    
    [parent, peer_node] = peer_cell_search(division_candidate_ids(i), ...
        division_candidate_ids, movieInfo, refine_res, embryo_vid, g);
    
    if ~isnan(parent)
        movieInfo.kids{parent} = [division_candidate_ids(i); peer_node];
        
        movieInfo.parents{division_candidate_ids(i)} = parent;
        movieInfo.parents{peer_node} = parent;
        parent_track_id = movieInfo.particle2track(parent);
        
        kid_track_id = movieInfo.particle2track([division_candidate_ids(i) ...
            peer_node]);
        
        kid_track_id(kid_track_id == parent_track_id) = [];
        
        movieInfo.tracks{parent_track_id} = unique(cat(1, ...
            movieInfo.tracks{parent_track_id}, movieInfo.tracks{kid_track_id}));
        
        movieInfo.tracks(kid_track_id) = {[]};
        
        movieInfo.particle2track(movieInfo.tracks{parent_track_id},1) = parent_track_id;
        movieInfo.particle2track(movieInfo.tracks{parent_track_id},2) = ...
            1:length(movieInfo.tracks{parent_track_id});
        display_track_link(movieInfo, parent_track_id,refine_res, embryo_vid);

        division_candidate_ids(division_candidate_ids==peer_node) = nan;
    end
end

end