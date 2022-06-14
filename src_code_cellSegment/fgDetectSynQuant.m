function [comMaps, repeatedSeed] = fgDetectSynQuant(comMapsIn, OrSt, q)
%Detect the foreground within a cropped region
% INPUT: comMaps which contains all infor about the cropped region
% contact: ccwang@vt.edu, 05/11/2020

repeatedSeed = [];
if ~iscell(comMapsIn)
    if strcmp(OrSt.imProcMethod, 'stb')
        vidIn = comMapsIn.vidStbComp;
    else
        vidIn = comMapsIn.vidComp;
    end
    % use synQuant to refine the region
    seedRegion = double(comMapsIn.regComp);
    L = seedRegion;
    other_id_map = comMapsIn.idComp>0 & ~comMapsIn.regComp;
    fmapCompInit = comMapsIn.fmapCompInit;
    [fg, pickedThreshold] = ordstat4fg(vidIn, ...
        comMapsIn.vid_sm, L, other_id_map, fmapCompInit, OrSt);
    
    comMaps = comMapsIn;
else
    if strcmp(OrSt.imProcMethod, 'stb')
        vidIn = cat(1, comMapsIn{1}.vidStbComp, ...
            comMapsIn{2}.vidStbComp, comMapsIn{3}.vidStbComp);
    else
        vidIn = cat(1, comMapsIn{1}.vidComp, ...
            comMapsIn{2}.vidComp, comMapsIn{3}.vidComp);
    end
    % use synQuant to refine the region
    seedRegion = double(cat(1, zeros(size(comMapsIn{1}.vidComp)), ...
        comMapsIn{2}.regComp, zeros(size(comMapsIn{3}.vidComp))));
    L = seedRegion;
    
    other_id_map = comMapsIn{2}.idComp>0 & ~comMapsIn{2}.regComp;
    other_id_map = cat(1, comMapsIn{1}.idComp>0, other_id_map, ...
        comMapsIn{3}.idComp>0);
    
    vid_sm = cat(1, comMapsIn{1}.vid_sm, comMapsIn{2}.vid_sm,...
        comMapsIn{3}.vid_sm);
    
    fmapCompInit = comMapsIn{2}.fmapCompInit;
    if ~isempty(comMapsIn{1}.vidComp)
        fmapCompInit = cat(1, comMapsIn{2}.fmapCompInit, fmapCompInit);
    end
    if ~isempty(comMapsIn{3}.vidComp)
        fmapCompInit = cat(1, fmapCompInit, comMapsIn{2}.fmapCompInit);
    end
    tmpOrSt = OrSt;
    tmpOrSt.stbVarCropMap = cat(1, OrSt.stbVarCropMap{:});
    tmpOrSt.NoStbVarCropMap = cat(1, OrSt.NoStbVarCropMap{:});
    if ~isfield(q, 'minSize')
        [fg, pickedThreshold] = ordstat4fg(vidIn, vid_sm,...
            L, other_id_map, fmapCompInit, tmpOrSt);
    else
        [fg, pickedThreshold] = ordstat4fg(vidIn, vid_sm,...
            L, other_id_map, fmapCompInit, tmpOrSt, q.minSize);
    end
    [h1,~,~] = size(comMapsIn{1}.vidComp);
    [h2,~,~] = size(comMapsIn{2}.vidComp);
    fg = fg(h1+1:h1+h2,:,:);
    %other_id_map = other_id_map(h1+1:h1+h2,:,:);
    seedRegion = seedRegion(h1+1:h1+h2,:,:);
    fmapCompInit = fmapCompInit(h1+1:h1+h2,:,:);
    comMaps = comMapsIn{2};
end
% if isfield(comMaps, 'real_regComp') % over-merged region
%     % only happen in redetectCellinTrackingwithSeed.m
%     other_id_map = comMaps.idComp>0 & ~comMaps.real_regComp;
% else
%     other_id_map = comMaps.idComp>0 & ~comMaps.regComp;
% end
other_id_map = comMaps.idComp>0 & ~comMaps.regComp;
%seedRegion(seedRegion>0) = comMaps.idComp(seedRegion>0); % re-label seed regions
comMaps.pickedThreshold = pickedThreshold;
if isnan(pickedThreshold) || isempty(find(fg, 1)) % no valid threshold/region can be found
    comMaps.fmapComp = fg;
    return;
end
%% after getting fg, update score maps
min_pv = min(comMaps.score3dMap(fg));
max_pv = max(comMaps.score3dMap(fg));
comMaps.score3dMap = scale_image(comMaps.score3dMap, 1e-3,1, min_pv, max_pv);
min_pv = min(comMaps.score2dMap(fg));
max_pv = max(comMaps.score2dMap(fg));
comMaps.score2dMap = scale_image(comMaps.score2dMap, 1e-3,1, min_pv, max_pv);

% fg = imdilate(seedRegion,strel("sphere",2));                   % debug
%% other seed regions
% all other seed regions are assigned to the same id
seedRegion(~fg) = 0; % it is possible the seed region is not fully covered
% the boundary of foreground should not be included in current region, it
% is too large. We can either: 
% (1) exclude it by labeling it as other id
% (2) re-define the foreground again
%fgboundary = (fmapCompInit - imerode(fmapCompInit, strel('disk',1)))>0;
[h,w,zslice] = size(fmapCompInit);
fgboundary = false(h,w,zslice);
for i=1:size(fgboundary,3)% does not consider the real boundary
    B = bwboundaries(fmapCompInit(:,:,i));
    bnd = cat(1, B{:});
    if isempty(bnd)
        continue;
    end
    idx = sub2ind([h,w], bnd(:,1), bnd(:,2)) + (i-1)*h*w;
    fgboundary(idx) = true;
    %figure;imshow(L>0);
    %         fgboundary(1,:,i) = false;
    %         fgboundary(:,end,i) = false;
    %         fgboundary(end,:,i) = false;
    %         fgboundary(:,1,i) = false;
end
fgboundary(seedRegion>0) = false;
comMaps.fmapCompInitBndIdx = find(fgboundary);

if isempty(find(fg(comMaps.fmapCompInitBndIdx),1))                          % if the foreground is in fmap
    other_id = other_id_map>0 & fg;                                         
else                                                                        % fg too large and touch the boundary
    if strcmp(q.fgBoundaryHandle, 'leaveAloneFirst')% way 3
        other_id = other_id_map>0 & fg;
    elseif strcmp(q.fgBoundaryHandle, 'compete') || nargout == 1 % way 1
        other_id = (other_id_map>0 | fgboundary) & fg;
    elseif strcmp(q.fgBoundaryHandle, 'repeat') % way 2
        repeatedSeed = fg;
        return;
    end
end
append_id = nan;
if ~isempty(find(other_id, 1))
    append_id = max(seedRegion(:)) + 1;
    seedRegion(other_id) = append_id;
end
%% test if there is an early stop
% debug
% curr_idMap = comMaps.idComp(fg);
% if length(unique(curr_idMap(curr_idMap>0))) == 1 && ...
%         length(unique(seedRegion(seedRegion>0)))>1
%     keyboard;
% end
    
comMaps = splitFGintoCells(fg, fgboundary, seedRegion, comMaps, append_id, q);


end