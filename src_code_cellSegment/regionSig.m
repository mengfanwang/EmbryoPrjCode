function [zscore,z_debias, sz_nei, sigma] = regionSig(regMap, vidComp, fmap, validNeiMap, OrSt)
%GETORDERSTAT Compute order statistics
% d: transformed data. smask: synapse mask.
% bmask: background mask, MUST be connected. imgVar: image variance
bnd_nei_radius = 3;
sph = strel('sphere', bnd_nei_radius);
se = strel(sph.Neighborhood(:,:,bnd_nei_radius+1));
fg_locs = regMap & ~imerode(regMap,se); % h*w*slice                         % out 3 slices of foreground

% find balanced neighbors
% absdif = 0;
% v=0;
if true
    bnd_nei_gap = 0;
    nei_locs = findValidNeiReg(regMap, validNeiMap, bnd_nei_gap,...
        bnd_nei_radius, true, se);% h*w*slice
    fg_locs = fg_locs & imdilate(nei_locs, se);                             % too far to background are removed
%     vv = vidComp(:,:,6);
%     vv = vidComp(:,:,6);
%     rr = regMap(:,:,);
%     disp([mean(fg_neighbors)-mean(fg), length(fg_neighbors), length(fg)]);
%     curVar = var(fg);%(var(fg)*length(fg) + var(fg_neighbors)*length(fg_neighbors))/(length(fg)+length(fg_neighbors));
%     
%     rr = length(fg) / length(find(imdilate(regMap, ones(3,3))));
%     lower = 1 - rr;%length(fg)/(length(fg)+length(fg_neighbors));
%     [~, ~, varCorrectionTerm] = truncatedGauss(0, 1, lower, inf);
%     curVar = curVar/varCorrectionTerm;
% 
%     %curVar = min(var(fg), var(fg_neighbors));
% %     disp([length(fg), length(fg_neighbors), length(fg)/length(fg_neighbors)]);
% %     disp([sqrt(curVar),mean(fg)-mean(fg_neighbors), (mean(fg)-mean(fg_neighbors))/sqrt(curVar)]);
%     M = max(length(fg), 10);
%     N = max(length(fg_neighbors), 10);
%     absdif = mean(fg)-mean(fg_neighbors);
%     if M>10*N
%         zscore = nan;
%     else
%         zscore = (mean(fg)-mean(fg_neighbors))/sqrt(curVar);
%     end
%     v(1) = var(fg);
%     v(2) = var(fg_neighbors);
%     v(3) = curVar;
%     v(4) = length(fg) / length(find(imdilate(regMap, ones(3,3))));
%     return;
%     zz = zeros(size(regMap));
%     zz(regMap>0) = 1;
%     zz(fg_locs) = 2;
%     zz(nei_locs) = 3;
%     zzshow(label2rgb3d(zz));
else
    [h,w,z] = size(regMap);
    % dilate the data (x,y direction) to remove boundary testing
    dilate_sc = 2;
    vid_di = nan(h+dilate_sc*2, w+dilate_sc*2, z+dilate_sc*2);
    vid_di(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc,...
        dilate_sc+1:end-dilate_sc) = vidComp;
    regMap_di = zeros(size(vid_di));
    regMap_di(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc,...
        dilate_sc+1:end-dilate_sc) = regMap;
    NeiMap_di = false(size(vid_di));
    NeiMap_di(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc,...
        dilate_sc+1:end-dilate_sc) = validNeiMap;

    [h,w,~] = size(vid_di);
    surfaceVoxels = find(regMap_di - imerode(regMap_di, strel('cube',3)));
    % zz = zeros(size(vid_di));
    % zz(regMap_di>0) = 1;
    % zz(surfaceVoxels) = zz(surfaceVoxels) + 2;
    % zzshow(label2rgb3d(zz));
    % fg_neighbors8 = [];
    % nei_mat8 = [];
    nei8 = [-h-1, -h, -h+1, -1, 1, h-1, h, h+1];                                                          
    neighbors = [nei8-h*w, -h*w, nei8, h*w, nei8+h*w];
    % nei25 = [-2*h-2:-2*h+2, -h-2, -h+2, -2,2, h-2,h+2, 2*h-2:2*h+2];
    % neighbors = [nei8-2*h*w,-2*h*w, nei25-h*w, nei25...
    %     , nei25+h*w, 2*h*w, nei8+2*h*w];
    % neighbors = [nei8-h*w, nei25, nei8+h*w];
    % neighbors = nei25;

    nei_mat = repmat(surfaceVoxels, 1, length(neighbors)) +...
        repmat(neighbors, length(surfaceVoxels), 1);
    nei_mat = nei_mat(:);
    nei_mat = nei_mat(NeiMap_di(nei_mat));
    nei_mat8 = unique(nei_mat);
    fg_neighbors8 = vid_di(nei_mat8);
    % 
    nei25 = [-2*h-2:-2*h+2, -h-2, -h+2, -2,2, h-2,h+2, 2*h-2:2*h+2];
    neighbors = [nei8-2*h*w,-2*h*w, nei25-h*w, -h*w, nei25, nei25+h*w, h*w, 2*h*w, nei8+2*h*w];

    nei_mat = repmat(surfaceVoxels, 1, length(neighbors)) +...
        repmat(neighbors, length(surfaceVoxels), 1);
    nei_mat = nei_mat(:);
    nei_mat = nei_mat(NeiMap_di(nei_mat));
    nei_mat25 = unique(nei_mat);
    %fg_neighbors25 = vid_di(nei_mat25);

    nei_locs = false(size(vid_di));
    nei_locs(nei_mat25) = true;
    if ~isempty(fg_neighbors8)
        nei_locs(nei_mat8) = false;
    %     zz = zeros(size(vid_di));
    %     zz(regMap_di>0) = 1;
    %     zz(tmp_l) = zz(tmp_l) + 2;
    %     zzshow(label2rgb3d(zz));
    end
    %fg_neighbors = vid_di(nei_locs);
    nei_locs = nei_locs(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc,...
                    dilate_sc+1:end-dilate_sc);
    % zz = zeros(size(vid_di));
    % zz(regMap_di>0) = 1;
    % zz(nei_mat8) = zz(nei_mat8) + 2;
    % zz(nei_mat25) = zz(nei_mat25) + 3;
    % zzshow(label2rgb3d(zz));
end

% [fgBgCells, fgBgLocCells] = localizeFgBgPairs(vidComp, regMap, ...
%     validNeiMap, 3, 20, 4);

fg = vidComp(fg_locs);
fg_neighbors = vidComp(nei_locs);
sz_nei = length(fg_neighbors);
% approximation of mu and sigma
M = max(length(fg), 10);
N = max(length(fg_neighbors), 10);

if M>10*N
    zscore = nan;
    z_debias = nan;
    sz_nei = nan;
    sigma = nan;
else
    if OrSt.fgVarSamplewiseEst
        % we need to estimate variance sample wise
        if strcmp(OrSt.imProcMethod,'stb')
            varMap = OrSt.stbVarCropMap;
        else
            varMap = OrSt.NoStbVarCropMap;
        end
        %OrSt.curStbVar = nanmean(varMap(regMap>0 | nei_locs)); %
        %OrSt.curStbVar = nanmean(varMap(fg_locs | nei_locs));
        OrSt.curStbVar = nanmean(varMap(fmap));
        %disp(sqrt(OrSt.curStbVar));
    end
    %fmap(regMap>0 | nei_locs) = false;
    fmap(fg_locs | nei_locs) = false;
    nanVec = vidComp(fmap);
    %nanVec = cat(1, nanVec, nanVec, nanVec);
    %nanVec = vidComp(regMap - fg_locs > 0);
    %nanVec = [];
    [zscore,z_debias, sigma] = ordStats(fg, fg_neighbors, nanVec, OrSt);
end
% 
% disp([mu sigma]);
% zscore2 = (sum_st - mu*sqrt(OrSt.NoStbVar))...
%     / (sigma*sqrt(OrSt.NoStbVar));
