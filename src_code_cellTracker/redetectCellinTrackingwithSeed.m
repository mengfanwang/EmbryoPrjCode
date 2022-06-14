function [comMaps, pseudo_seed_label] = redetectCellinTrackingwithSeed(...
    movieInfo, refine_res, embryo_vid,eigMaps, varAllMaps, frame, yxz, q)
% detect a missing cell given a seed region, this function is quite similar
% to the function 'refineOneRegion_with_seed.m', but does not have the gap
% testing and shrinking steps
if iscell(yxz) % we have multiple seeds
    seed_yxz = cat(1, yxz{:});
    pseudo_seed_label = zeros(numel(yxz), 1);
    for i=1:numel(yxz)
        pseudo_seed_label(i) = numel(movieInfo.voxIdx) + i;
        refine_res{frame}(yxz{i}) = pseudo_seed_label(i);
    end
else
    seed_yxz = yxz;
    pseudo_seed_label = numel(movieInfo.voxIdx) + 1;
    refine_res{frame}(seed_yxz) = pseudo_seed_label;
end

if isempty(seed_yxz)
    comMaps.pickedThreshold = nan;
    comMaps.fmapComp = [];
    return;
end
eig2d = eigMaps{frame}{1};
eig3d = eigMaps{frame}{2};
if q.multi_frames_flag
    vid = cell(3,1);
    varMap = cell(3,1);
    vid_stb = cell(3,1); % no use
    idMap = cell(3,1);
    for j=frame-1:frame+1
        tmp_idx = j-frame+2;
        if j>0 && j<=numel(refine_res)
            vid{tmp_idx} = embryo_vid{j};
            varMap{tmp_idx} = varAllMaps{j};
            idMap{tmp_idx} = refine_res{j};
        end
    end
else
    vid = embryo_vid{frame};
    varMap = varAllMaps{frame};
    idMap = refine_res{frame};
    
    vid_stb = [];
end

OrSt = inital_Orst(varMap);
if strcmp(OrSt.imProcMethod, 'stb')
    if iscell(vid)
        vid_stb = cellfun(@(x) sqrt(x+3/8), vid,'UniformOutput', false);
    else
        vid_stb = sqrt(vid+3/8);
    end
end
[comMapsInit, OrSt] = get_local_area(vid, vid_stb, ...
    idMap, pseudo_seed_label(1),...
    eig2d, eig3d, seed_yxz, OrSt, q);
if iscell(comMapsInit)
    init_seed_map = comMapsInit{2}.regComp;
else
    init_seed_map = comMapsInit.regComp;
end
% this may be not needed
% if iscell(comMapsInit)
%     if isfield(comMapsInit{2}, 'seed_map4fg')
%         comMapsInit{2}.real_regComp = comMapsInit{2}.regComp;
%         comMapsInit{2}.regComp = comMapsInit{2}.seed_map4fg;
%     end
% else
%     if isfield(comMapsInit, 'seed_map4fg')
%         comMapsInit.real_regComp = comMapsInit.regComp;
%         comMapsInit.regComp = comMapsInit.seed_map4fg;
%     end
% end

% should we remove other well dectected cells?
%     comMaps3tps{2}.fmapCompInit = comMaps3tps{2}.fmapCompInit &...
%         ~(comMaps3tps{2}.idComp>0 & ~comMaps3tps{2}.regComp);
comMaps = fgDetectSynQuant(comMapsInit, OrSt, q);
if isnan(comMaps.pickedThreshold) || isempty(find(comMaps.fmapComp, 1))
    return;
end

% test if the foreground is too small (this is very simplified version, if 
% we have extremely crowded cells, the fg will have to be the whole FOV. 
% If so, this way to address it is not enough.)
z_not_enough = false;
if ~isempty(find(comMaps.regComp(:,:,1) > 0, 1)) || ...
        ~isempty(find(comMaps.regComp(:,:,end) > 0, 1))
    z_not_enough = true;
end
xy_not_enough = false;
if ~isempty(find(comMaps.regComp(comMaps.fmapCompInitBndIdx)>0,1))
    xy_not_enough = true;
end
if  strcmp(q.fgBoundaryHandle, 'leaveAloneFirst') && ...
        (z_not_enough || xy_not_enough)
    if z_not_enough
        q.shift(3) = q.shift(3)*2;
    end
    if xy_not_enough
        q.shift(1) = q.shift(1)*2;
        q.shift(2) = q.shift(2)*2;
    end
    if iscell(vid)
        comMapsNew = get_local_area(vid{2}, vid_stb{2}, idMap{2}, pseudo_seed_label(1),...
            eig2d, eig3d, seed_yxz, OrSt, q);
    else
        comMapsNew = get_local_area(vid, vid_stb, idMap, pseudo_seed_label(1),...
            eig2d, eig3d, seed_yxz, OrSt, q);
    end
    init_seed_map = comMapsNew.regComp;
    %     if iscell(yxz) % re-define seed region as all
    %         comMapsNew.real_regComp = comMapsNew.regComp;
    %         comMapsNew.regComp = comMapsNew.seed_map4fg;
    %     end
    % we still based on the old threshold to get the foreground region;
    % other steps are the same
    comMaps = fgDetectSynQuant_thresGiven(comMapsNew, ...
        comMaps.pickedThreshold, q);
    if isnan(comMaps.pickedThreshold) || isempty(find(comMaps.fmapComp, 1))
        return;
    end
end
% if ~iscell(yxz)
%     comMaps.idComp(~comMaps.fmapComp) = 0;
% end

comMaps.init_seed_map = init_seed_map;
end