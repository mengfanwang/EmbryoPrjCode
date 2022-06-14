function [newLabel, comMaps, fgReDo] = refineOneRegion_with_seed(seed_id, ...
    yxz, vid, vid_stb, idMap, eig2d, eig3d, OrSt, q)
% refine the seed region indicate by index in yxz
% NOTE: it is possible the yxz is not exactly the same as seed_id label
% indicates, e.g. the seed is relabeled. We view yxz as the correct one.
%% crop needed information
[comMapsInit, OrSt] = get_local_area(vid, vid_stb, idMap, seed_id,...
    eig2d, eig3d, yxz, OrSt, q);
%     debug
%     locs = find(comMaps3tps{2}.fmapCompInit);
%     [~,~,zz] = ind2sub(size(comMaps3tps{2}.fmapCompInit), locs);
%     zz = unique(zz(zz>0));
%     if min(zz) > 1 || max(zz) < size(comMaps3tps{2}.fmapCompInit,3)
%         keyboard;
%     else
%         continue;
%     end
%% grow the region firstly immediately after synQuant
[comMaps, repeatedSeed] = fgDetectSynQuant(comMapsInit, OrSt, q);
if isempty(find(comMaps.fmapComp, 1)) % no valid threshold/region can be found
    newLabel = seed_id;
    fgReDo = false;
    return;
end
if strcmp(q.fgBoundaryHandle, 'repeat') && ~isempty(repeatedSeed)
    % fg is too small, we choose way 1 ('repeat'): to re-define fg based on
    % the new boundary-touched seed. NOTE: this could be too liberal
    % because the boundary may be touched by another cell's seed, so I
    % prefer use way 3: leaveAloneFirst
    yxz_enlarged = comMaps.linerInd(repeatedSeed>0);
    [comMapsInit, OrSt] = get_local_area(vid, vid_stb, idMap, seed_id,...
        eig2d, eig3d, yxz_enlarged, OrSt, q);
    comMaps = fgDetectSynQuant(comMapsInit, OrSt, q);
end
if iscell(OrSt.stbVarCropMap)
    OrSt.stbVarCropMap = OrSt.stbVarCropMap{2};
    OrSt.NoStbVarCropMap = OrSt.NoStbVarCropMap{2};
end
%     if true
%         writeRefineCell(comMaps.vid_sm, double(comMaps.fmapComp),...
%             comMaps.regComp, i, save_folder);
%         continue;
%     end
%% over-segment regions first and try merge them
[newLabel, ~] = segmentCurrentRegion(comMaps, seed_id, q, OrSt);
%% refine segmentation results (TODO: foreground detection in 'region_refine' is duplicated)
newLabel = region_refineV2(newLabel, comMaps, q);
%% remove regions unrelated with region i and small regions
seed_ids = newLabel(comMaps.idComp == seed_id);
seed_ids = unique(seed_ids(seed_ids>0));
% there also can be pixels not belonging to foreground
newLabel(~ismember(newLabel, seed_ids) | ~comMaps.fmapCompInit) = 0;
[newLabel, ~] = region_sanity_check(newLabel, q.minSize); % previously use 20

fgReDo = false;
z_not_enough = false;

if ~isempty(find(newLabel(:,:,1) > 0, 1)) || ...                            % label touch z boundary
        ~isempty(find(newLabel(:,:,end) > 0, 1))
    z_not_enough = true;
end
xy_not_enough = false;                                                      % label touch initial fmap boundary
if isfield(comMaps, 'fmapCompInitBndIdx')
    if ~isempty(find(newLabel(comMaps.fmapCompInitBndIdx)>0,1))
        xy_not_enough = true;
    end
end
%% test if the initial foreground is too small                              % if foreground is not enough, enlarge it and redo
if  strcmp(q.fgBoundaryHandle, 'leaveAloneFirst') && ...
        (z_not_enough || xy_not_enough)
    fgReDo = true;
    % fg is too small, we choose way 3 : to re-define fg based on
    % the new boundary-touched seed. NOTE:one time is enough
    %yxz_enlarged = comMaps.linerInd(newLabel>0);
    if z_not_enough
        q.shift(3) = q.shift(3)*2;
    end
    if xy_not_enough
        q.shift(1) = q.shift(1)*2;
        q.shift(2) = q.shift(2)*2;
    end
    if iscell(vid) % mutliple frames together
        [comMapsNew, OrSt] = get_local_area(vid{2}, vid_stb{2}, idMap{2}, ...
            seed_id, eig2d, eig3d, yxz, OrSt, q);
        OrSt.stbVarCropMap = OrSt.stbVarCropMap{2};
        OrSt.NoStbVarCropMap = OrSt.NoStbVarCropMap{2};
    else
        [comMapsNew, OrSt] = get_local_area(vid, vid_stb, idMap, seed_id,...
            eig2d, eig3d, yxz, OrSt, q);
    end
    % we still based on the old threshold to get the foreground region;
    % other steps are the same
    comMaps = fgDetectSynQuant_thresGiven(comMapsNew, ...
        comMaps.pickedThreshold, q);
    % over-segment regions first and try merge them
    [newLabel, ~] = segmentCurrentRegion(comMaps, seed_id, q, OrSt);
    % refine segmentation results (TODO: foreground detection in 'region_refine' is duplicated)
    newLabel = region_refineV2(newLabel, comMaps, q);
    % remove regions unrelated with region i and small regions
    seed_ids = newLabel(comMaps.idComp == seed_id);
    seed_ids = unique(seed_ids(seed_ids>0));
    % there also can be pixels not belonging to foreground
    newLabel(~ismember(newLabel, seed_ids) | ~comMaps.fmapCompInit) = 0;
    [newLabel, ~] = region_sanity_check(newLabel, q.minSize); % previously use 20
end
