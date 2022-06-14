function comMaps = cropNeedRegions(vid, vid_stb, idMap, reg_id, eig2d, ...
    eig3d, yxz, improc, shift)
% crop the info we need that only contains the targeted region(s)
% region is defined by yxz and shift
comMaps = [];
if ~isempty(vid)
    [comMaps.vidComp,comMaps.linerInd, yxz_edge_Flag, loc_org_xyz] = crop3D(vid, yxz, shift);
    comMaps.yxz_edge_Flag = yxz_edge_Flag;
    if ~isempty(vid_stb)
        comMaps.vidStbComp = crop3D(vid_stb, yxz, shift);
    else
        comMaps.vidStbComp=[];
    end
    
    comMaps.idComp = crop3D(idMap, yxz, shift);
    sm_term = [1 1 1];
    if strcmp(improc,'stb') % stablized smoothed data
        comMaps.vid_sm = imgaussfilt3(comMaps.vidStbComp, sm_term);
    else% original smoothed data
        comMaps.vid_sm = imgaussfilt3(comMaps.vidComp, sm_term);
    end
else
    comMaps.vidStbComp = [];
    comMaps.vidComp = [];
    comMaps.idComp = [];
    comMaps.vid_sm = [];
end
if isnan(reg_id)
    return;
end
% L = bwlabeln(comMaps.idComp == reg_id, 26);
% s = regionprops3(L, 'Volume');
% if numel(s) == 1
%     comMaps.regComp = comMaps.idComp == reg_id;
% else
%     %fprintf('the seed region contains %d connected components!\n', numel(s));
%     [~, od] = max([s.Volume]);
%     comMaps.regComp = L==od;
%     comMaps.idComp(L>0 & L~=od) = 0;
% end
reg_map = comMaps.idComp == reg_id;
[comMaps.regComp, rm_flag] = pickLargestReg(reg_map, 26);%
% if rm_flag % reg_id contains multiple connected components
%     keyboard;
% end
% 
%% build the foreground map                                                % This part include two regions, a larger one(yxz) and a smaller
seed_idx_in_cropReg = coordinate_transfer(yxz, size(vid), loc_org_xyz, ... % one(regComp). The candidate foreground is between the two
    size(comMaps.idComp));
seed_map4fg = comMaps.regComp;
if ~isempty(find(comMaps.idComp(seed_idx_in_cropReg)~=reg_id, 1))            % exist region ~= reg_id
    % the idx yxz is larger the area covered by reg_id in label map, we add  
    % yxz as our standard to detect the foreground (and only for foreground)
    seed_map4fg(seed_idx_in_cropReg) = true;                                 % if exist other region, we select all but pick the largest
    seed_map4fg = pickLargestReg(seed_map4fg, 26);
    if ~isempty(find(seed_map4fg & comMaps.regComp,1))...
        && ~isempty(find(seed_map4fg & ~comMaps.regComp,1))
        %comMaps.seed_map4fg = seed_map4fg;
        comMaps.regComp = seed_map4fg;                                       % update regComp
    end
%     
%     tmp = comMaps.idComp(seed_map4fg);
%     comMaps.idComp(seed_idx_in_cropReg) = 0;
%     comMaps.idComp(seed_map4fg) = tmp;
end
radius = shift(1);                                                           % dilate fmap
sph = strel('sphere', radius+1);
cycleSe = strel(sph.Neighborhood(:,:,radius+2));
fmap1moreCycle = imdilate(seed_map4fg, cycleSe);%crop3D(fmap, yxz, shift);

sph = strel('sphere', radius);
cycleSe = strel(sph.Neighborhood(:,:,radius+1));
fmap = imdilate(seed_map4fg, cycleSe);%crop3D(fmap, yxz, shift);
locs = find(seed_map4fg);
[~,~,zz] = ind2sub(size(seed_map4fg), locs);
zStack = unique(zz);
minZ = zStack(1); maxZ = zStack(end);
for i=zStack(1)-1:-1:max(1, zStack(1)-shift(3))   %z direction dilate
    fmap1moreCycle(:,:,i) = fmap1moreCycle(:,:,zStack(1));
    fmap(:,:,i) = fmap(:,:,zStack(1));
    minZ = i;
end
for i=zStack(end)+1:min(size(seed_map4fg,3), zStack(end)+shift(3))
    fmap1moreCycle(:,:,i) = fmap1moreCycle(:,:,zStack(end));
    fmap(:,:,i) = fmap(:,:,zStack(end));
    maxZ = i;
end
% remove >5% low-intensity pixels
comMaps.fmapCompInit = false(size(fmap));
for i=minZ:maxZ
    curfMap = fmap1moreCycle(:,:,i);
    rmLen = length(find(curfMap & ~fmap(:,:,i)));
    curVidSm = comMaps.vid_sm(:,:,i);
    vals = sort(curVidSm(curfMap));%ascending
    if isempty(vals) % use previous frame's result
        %warning('seed region is not z-continuous!');
        if i>minZ
            comMaps.fmapCompInit(:,:,i) = fmap(:,:,i-1);
        elseif i<maxZ % never reach
            comMaps.fmapCompInit(:,:,i) = fmap(:,:,i+1);
        else % never reach
            comMaps.fmapCompInit(:,:,i) = fmap(:,:,i);
        end
    else
        curThre = vals(max(rmLen, round(0.05*length(vals))));
        comMaps.fmapCompInit(:,:,i) = ...
            imopen(curVidSm>curThre & curfMap, strel('disk', 1));
    end
end
%comMaps.fmapCompInit = comMaps.vid_sm>0;
%% build 3d/2d principal curveture score map
eig3dComp = crop3D(eig3d, yxz, shift);
comMaps.eigPosMap = eig3dComp>0; % gaps detected by principal curvature
eig3dComp(comMaps.eigPosMap<0) = 0;

comMaps.score3dMap = eig3dComp;%scale_image(eig3dComp, 1e-3,1);
comMaps.score3dMap(isnan(comMaps.score3dMap)) = 0;

eig2dComp = crop3D(eig2d, yxz, shift);
comMaps.eig2dPosMap = eig2dComp>0; % gaps detected by principal curvature
eig2dComp(eig2dComp<0) = 0;
comMaps.score2dMap = eig2dComp;%scale_image(eig2dComp, 1e-3,1);
comMaps.score2dMap(isnan(comMaps.score2dMap)) = 0;


comMaps.newIdComp = comMaps.regComp;
comMaps.newIdComp(comMaps.eigPosMap) = 0; % if we segment current region further

end