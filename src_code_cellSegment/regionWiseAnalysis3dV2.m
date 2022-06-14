function newIdMap = regionWiseAnalysis3dV2(idMap, eigMap,vid, ...
    fmap, varMap, save_folder)
% task of this function:
% 1. refine the region from overall synQuant
% 2. test the gap inside one given region, if it is significnat, cut the
% region into two cells
%INPUT:
% idMap: the label map of h*w*z
% eigMap: the principal curvature values of foreground in idMap include 2d
% and 3d (cells with two elements)
% eigThres: the threshold of principal curvature for eigMap
% vid: the orgriginal data
% save_folder: the folder to save images for sanity check
% over_seg_flag: true, oversegment the region and then merge, false
% otherwise
%OUTPUT:
% newIdMap: detected regions with each region one unique id
%
% contact: ccwang@vt.edu, 05/05/2020

if nargin == 4
    save_folder = [];
end
%% parameters for region grow
q.cost_design = [1, 2];
q.growConnectInTest = 4;
q.growConnectInRefine = 6;
q.edgeConnect = 48; % three circle of neighbors
% deal with one component each time
q.neiMap = 26;
q.shift = [10 10 2];
q.shrinkFlag = true;
%% parameters for fg detection and gap significance test
% data stablization
vid_stb = sqrt(vid+3/8);
% % temporal variance
%[sigma2, xx1] = temporalVar(vid, vid2, 1);
% use noise variance
OrSt.gapTestWay = 'localOrdStats';%'orderStats';%
%OrSt.stbVarMap = calVarianceStablizationBY(vid_stb, 0.8);
OrSt.stbVarMap = varMap{1,2};

fmap = fmap & ~isnan(OrSt.stbVarMap);
OrSt.corrTerm = max(vid_stb(:))/255;
% figure;plot(xx1), hold on; plot(xx2);legend('t', 's');
%OrSt.noiseVar = calVarianceStablizationBY(vid_stb, 0.8);
OrSt.mu = [];
OrSt.sigma = [];
OrSt.p_thres = 0.01; % pvalue threshold
OrSt.imProcMethod = 'noStb';% 'stb' or 'noStb'
OrSt.fgVarSamplewiseEst = true;% or 'noStb'
OrSt.fgTestWay = 'KSecApprox';
%[OrSt.NoStbVarMap, OrSt.NoStbVar] = calVarianceStablizationBY(vid, 0.8);
OrSt.NoStbVarMap = varMap{1,1};
OrSt.NoStbVar = varMap{2,1};
fmap = fmap & ~isnan(OrSt.NoStbVarMap);
if strcmp(OrSt.fgTestWay, 'lookupTable')
    muSigma = paraP3D(0.05, 0,0,255,4, 500, 0.1, 10);
    OrSt.NoStbMu = muSigma.mu;
    OrSt.NoStbSig = muSigma.sigma;
end
%% eigenvalue map 2d and 3d
eig2d = eigMap{1};
eig3d = eigMap{2};
%% start the main functions
[idMap2, redundant_flag] = region_sanity_check(idMap, 10);
if redundant_flag
    idMap = idMap2;
end
s = regionprops3(idMap, {'VoxelIdxList'});
loc_cells = cell(numel(s.VoxelIdxList),1);
gapTestSamples = false(numel(s.VoxelIdxList),1);% for debug
regCnt = 0;

for i=40:numel(s.VoxelIdxList) %[1-29 2-40 3-38]
    yxz = s.VoxelIdxList{i};
    fprintf('process %d out of %d\n', i, numel(s.VoxelIdxList));
    %% crop needed information
    comMaps = cropNeedRegions(vid, vid_stb, idMap, i, eig2d, eig3d, yxz, OrSt.imProcMethod, q.shift);
    OrSt.stbVarCropMap = crop3D(OrSt.stbVarMap, yxz, q.shift);
    OrSt.NoStbVarCropMap = crop3D(OrSt.NoStbVarMap, yxz, q.shift);
    %% grow the region firstly immediately after synQuant
    comMaps = fgDetectSynQuant(comMaps, OrSt, q);
    %% over-segment regions first and try merge them
    [newLabel, edgeLabel, gapTestSamples(i)] = segmentCurrentRegion(comMaps, i, q, OrSt);
    %% refine segmentation results (TODO: foreground detection in 'region_refine' is duplicated)
    newLabel = region_refineV2(newLabel, comMaps, q);
    %     newLabel = region_refine(newLabel,...
    %         comMaps.vidComp, comMaps.idComp, i, comMaps.fmapComp, comMaps.score2dMap+comMaps.score3dMap,...
    %         OrSt, q.growConnectInRefine, q.cost_design);
    %% remove regions unrelated with region i and small regions
    seed_ids = newLabel(comMaps.idComp == i);
    seed_ids = unique(seed_ids(seed_ids>0));
    % there also can be pixels not belonging to foreground
    newLabel(~ismember(newLabel, seed_ids) | ~comMaps.fmapCompInit) = 0;
    [newLabel, ~] = region_sanity_check(newLabel, 20);
    %% save cell locations
    valid_newLabel = newLabel(:)>0;
    loc_cells{i} = [comMaps.linerInd(valid_newLabel), newLabel(valid_newLabel)];
    regCnt = regCnt + max(loc_cells{i}(:,2));
%     if ~isempty(find(loc_cells{i}(:,1)==4070913,1))
%         keyboard;
%     end
    %% write data into images for error check
    if ~isempty(save_folder)
        writeRefineCell(comMaps.vidComp, newLabel, comMaps.regComp, i, save_folder);
    end
end
gapTestSamples = find(gapTestSamples); 
% label the new ID Map
newIdMap = zeros(size(idMap));
regCnt = 0;
for i=1:numel(loc_cells)
    if ~isempty(loc_cells{i})
        cur_locs = loc_cells{i}(:,1);
        cur_labels = loc_cells{i}(:,2);
        newIdMap(cur_locs) = cur_labels + regCnt;
        regCnt = regCnt + max(cur_labels);
    end
        
end
