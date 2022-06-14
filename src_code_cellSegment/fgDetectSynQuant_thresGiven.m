function comMaps = fgDetectSynQuant_thresGiven(comMaps, threshold, q)
%Detect the foreground within a cropped region
% INPUT: comMaps which contains all infor about the cropped region
% contact: ccwang@vt.edu, 05/11/2020
comMaps.pickedThreshold = threshold;
fgInit = comMaps.fmapCompInit & comMaps.vid_sm >= comMaps.pickedThreshold;
seedRegion = double(comMaps.regComp);
fg = refine_with_seedRegion(fgInit, seedRegion, q.minSize);
if isempty(find(fg & seedRegion,1))% the seedregion no longer related to fg
    comMaps.pickedThreshold = nan;
    comMaps.fmapComp = false(size(fg));
    return;
end
%% after getting fg, update score maps
min_pv = min(comMaps.score3dMap(fg));
max_pv = max(comMaps.score3dMap(fg));
comMaps.score3dMap = scale_image(comMaps.score3dMap, 1e-3,1, min_pv, max_pv);
min_pv = min(comMaps.score2dMap(fg));
max_pv = max(comMaps.score2dMap(fg));
comMaps.score2dMap = scale_image(comMaps.score2dMap, 1e-3,1, min_pv, max_pv);
%% other seed regions
% if isfield(comMaps, 'real_regComp') % over-merged region
%     % only happen in redetectCellinTrackingwithSeed.m
%     other_id_map = comMaps.idComp>0 & ~comMaps.real_regComp;
% else
%     other_id_map = comMaps.idComp>0 & ~comMaps.regComp;
% end
other_id_map = comMaps.idComp>0 & ~comMaps.regComp;

seedRegion(~fg) = 0; % it is possible the seed region is not fully covered.

% the boundary of foreground should not be included in current region, it
% is too large. We can either: 
% (1) exclude it by labeling it as other id
% (2) re-define the foreground again
[h,w,zslice] = size(fg);
fgboundary = false(h,w,zslice);
for i=1:size(fgboundary,3)% does not consider the real boundary
    B = bwboundaries(comMaps.fmapCompInit(:,:,i));
    bnd = cat(1, B{:});
    if isempty(bnd)
        continue;
    end
    idx = sub2ind([h,w], bnd(:,1), bnd(:,2)) + (i-1)*h*w;
    fgboundary(idx) = true;
end
fgboundary(seedRegion>0) = false;
% consider boundary as sink: this foreground should be large enough
other_id = (other_id_map>0 | fgboundary) & fg;
% all other seed regions are assigned to the same id
append_id = nan;
if ~isempty(find(other_id, 1))
    append_id = max(seedRegion(:)) + 1;
    seedRegion(other_id) = append_id;
end
%% test if there is an early stop
comMaps = splitFGintoCells(fg, fgboundary, seedRegion, comMaps, append_id, q);
% val_id_loc = find(seedRegion,1);
% if ~isempty(val_id_loc)
%     reg1stId = seedRegion(val_id_loc);
% else
%     reg1stId = 0;
% end
% 
% if ~isempty(find(seedRegion~=reg1stId & seedRegion>0,1))
%     % there is multiple seeds
%     fgIn = fg;
%     % 2d first?
%     newLabel = regionGrow(seedRegion, ...
%         comMaps.score2dMap,...
%         fgIn, q.growConnectInTest, ...
%         q.cost_design, true);
%     newLabel = regionGrow(newLabel, ...
%         comMaps.score2dMap+comMaps.score3dMap,...
%         fgIn, q.growConnectInRefine, ...
%         q.cost_design, false);
%     
%     newLabel(~fg) = 0;
%     % here for extra voxels
%     if ~isempty(find(fg & (newLabel==0), 1))
%         %disp('we find extra voxels to re-assign');
%         newLabel = extraVoxReassign(newLabel, fg);
%     end
%     if ~isnan(append_id)
%         newLabel(newLabel == append_id) = 0;
%     end
%     fg = newLabel>0;
% end
% comMaps.regComp = fg;
% 
% comMaps.fmapComp = fg;
% comMaps.newIdComp = fg;
% comMaps.newIdComp(comMaps.eigPosMap) = 0;


end