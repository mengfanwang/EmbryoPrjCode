function [newLabel, fg] = region_refine(inLabel, L, ...
     vidComp, idComp, reg_id, fmapComp,  scoreComp, ...
    OrSt, connect, cost_design, bg2sinkLink)
% refine the region label with shrink and graph-cut
% INPUT:
% newLabel: the label of regions segmented by principal curvature
% idComp: the labels of regions
% regTest: gaps detected
% scoreComp: score that used to grow the region, can be principal curvature
% or the gradient. Here we use gradient
% fMap: the foreground map to indicate the region can be grown on
% connect: 6, 26 (TODO: 10(8+2)) connection for connecting edges in the graph
% cost_design: how to design the cost in the Graph-Cut
%OUTPUT:
% newLabel: new segmentation after region growing
%
% contact: ccwang@vt.edu, 02/13/2020
seedRegion = double(inLabel);
%% use synQuant to refine the region
if nargin == 10
    bg2sinkLink = false; % this means we separate fg based on seed and pc
end
if mean(vidComp(seedRegion>0)) > mean(vidComp(L>0))
    L = seedRegion;
end
other_id_map = idComp;
other_id_map(idComp==reg_id) = 0;
fg = ordstat4fg(vidComp, L, other_id_map, fmapComp, OrSt);
%% other seed regions
% all other seed regions are assigned to the same id
seedRegion(~fg) = 0; % it is possible the seed region is not fully covered
other_id = other_id_map>0 & fg;
append_id = nan;
if ~isempty(find(other_id, 1))
    append_id = max(seedRegion(:)) + 1;
    seedRegion(other_id) = append_id;
end
%% test if there is an early stop
val_id_loc = find(seedRegion,1);
if isempty(val_id_loc) % no seed region or fg is blank
    newLabel = double(fg);
    %n = 0;
elseif isempty(find(seedRegion~=seedRegion(val_id_loc) & seedRegion>0,1))
    newLabel = double(fg); % only one valid seed, so all fg is assigned to it
    %n = 1;
else % there is multiple seeds
    fgIn = fg;%imdilate(fg, strel('cube', 3));
    newLabel = regionGrow(seedRegion, scoreComp, fgIn, connect, cost_design,bg2sinkLink);
    %newLabel = regionGrow(seedRegion, scoreComp, fg, connect, cost_design);
    %return;
    newLabel(~fg) = 0;
    if ~isempty(find(fg & (newLabel==0), 1))
        disp('we find extra voxels to re-assign');
        newLabel = extraVoxReassign(newLabel, fg);
    end

    if ~isnan(append_id)
        newLabel(newLabel == append_id) = 0;
    end
end
%% shrink the region to see if more can be found
if nargout == 1 % >1 means we are for foreground detection
    shrinkScale = [3,0];
    [newLabel, ~] = shrinkTest3d(newLabel, shrinkScale);
end


end