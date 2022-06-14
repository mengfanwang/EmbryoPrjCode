function newLabel = region_refineV2(inLabel, comMaps, q, bg2sinkLink)
% second version: refine the region label with shrink and graph-cut
% INPUT:
% inLabel: the label of regions segmented by principal curvature
% L:
% contact: ccwang@vt.edu, 05/11/2020
if nargin==3
    bg2sinkLink = false; % not link background to sink
end

fg = comMaps.fmapComp;
seedRegion = inLabel;
%% there will be no other seed regions                                      % mask other id as append id
% all other seed regions are assigned to the same id
seedRegion(~fg) = 0; % it is possible the seed region is not fully covered
other_id_map = comMaps.idComp>0 & ~comMaps.regComp;
other_id = other_id_map & fg;
append_id = nan;
if ~isempty(find(other_id, 1))
    append_id = max(seedRegion(:)) + 1;
    seedRegion(other_id) = append_id;
end
%% test if there is an early stop                                           % implement region growing if multilabels
val_id_loc = find(seedRegion,1);
if isempty(val_id_loc) % no seed region or fg is blank
    newLabel = double(fg);
elseif isempty(find(seedRegion~=seedRegion(val_id_loc) & seedRegion>0,1))
    newLabel = double(fg); % only one valid seed, so all fg is assigned to it
else % there is multiple seeds
    newLabel = regionGrow(seedRegion, ...
        comMaps.score2dMap+comMaps.score3dMap,...
        fg, q.growConnectInRefine, ...
        q.cost_design, bg2sinkLink);
    
    newLabel(~fg) = 0;
    % here for extra voxelsc
    if ~isempty(find(fg & (newLabel==0), 1))
        %disp('we find extra voxels to re-assign');
        newLabel = extraVoxReassign(newLabel, fg);
    end
    if ~isnan(append_id)
        newLabel(newLabel == append_id) = 0;
    end
end
%% shrink the region to see if more can be found
if q.shrinkFlag
    shrinkScale = [2,1];
    [outLabel, reg_split] = shrinkTest3d(newLabel, shrinkScale, comMaps, q);
    if reg_split
        newLabel = outLabel;
    end

end

end