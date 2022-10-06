function [newIdMap, thresholdMap] = m_regionWiseAnalysis4d(idMap, ...
    eigMap, vid, varMap, test_ids, tif_id)
% region grow: shrink the 3d gaps detected by 3d principal curvature
% considering the 4-d information (pre- and post frames)
% the idea is based on boyu's graph-cut (max-flow)
% INPUT:
% idMap: cells of 3 components, each of which is the label map of h*w*z
% eig2dMaps: cells of 3 components, each of which is the 2d principal curvature
% values of foreground in idMap
% eigThres: the threshold of principal curvature for eigMap
% vid: cells of 3 components, each of which is the orgriginal data
% save_folder: the folder to save images for sanity check
% over_seg_flag: true, oversegment the region and then merge, false
% otherwise
% OUTPUT:
% newIdMap: detected regions with each region one unique id
%
% contact: ccwang@vt.edu, 09/15/2020
if nargin == 4
    test_ids = [];
    tif_id = [];
    %save_folder = [];
end
if nargin == 5
    %save_folder = [];
    tif_id = [];
end
%% parameters for region grow
q = initial_q();
%% parameters for fg detection and gap significance test
OrSt = inital_Orst(varMap);
% data stablization
if iscell(vid)
    vid_stb = cellfun(@(x) sqrt(x+3/8), vid,'UniformOutput', false);
    if strcmp(OrSt.imProcMethod, 'stb')
        OrSt.corrTerm = max(vid_stb{2}(:))/255;
    end
    % eigenvalue map 2d and 3d: we only use current principal curvature
    eig2d = eigMap{2}{1};
    eig3d = eigMap{2}{2};
    idMap_current = idMap{2};
    vid_current = vid{2};
else
    vid_stb = sqrt(vid+3/8);
    if strcmp(OrSt.imProcMethod, 'stb')
        OrSt.corrTerm = max(vid_stb(:))/255;
    end
    % eigenvalue map 2d and 3d: we only use current principal curvature
    eig2d = eigMap{1};
    eig3d = eigMap{2};
    idMap_current = idMap;
    vid_current = vid;
end
idMap_init = idMap_current;

%% start the main functions
if isempty(test_ids)
    [idMap2, redundant_flag] = region_sanity_check(idMap_current, q.minSeedSize); % remove small objects
    if redundant_flag
        idMap_current = idMap2;
    end
    s = regionprops3(idMap_current, {'VoxelIdxList'});
else
    idMap2 = idMap_current;
    idMap2(~ismember(idMap2, test_ids)) = 0;
    s = regionprops3(idMap2, {'VoxelIdxList'});
    
    idMap2 = idMap_current;
    idMap2(ismember(idMap2, test_ids)) = 0;
    s_1st = regionprops3(idMap2, {'VoxelIdxList','Volume'});
    % build prior-knowledge information from detected cells
    cell_sz = [s_1st.Volume];
    cell_levels = cellfun(@(x) mean(vid_current(x)), s_1st.VoxelIdxList);
    cell_sz(test_ids) = []; % remove cells from 2nd iteration
    cell_levels(test_ids) = []; % remove cells from 2nd iteration
    od = round(0.05*length(cell_sz));
    st_levels = sort(cell_levels, 'ascend');
    min_level = st_levels(od);
    st_szs = sort(cell_sz, 'ascend');
    minSz = st_szs(od);
end
loc_cells = cell(numel(s.VoxelIdxList),1);
% process from brightest seed to dimmest ones
seed_levels = cellfun(@(x) mean(vid_current(x)), s.VoxelIdxList);
[~, seed_proc_order] = sort(seed_levels,'descend');
fgReDoCnt = 0;
for i=1:numel(s.VoxelIdxList)%[7 9 13 15 17 18 19 20 21 22]%
    fprintf('%d/%d ', i, numel(s.VoxelIdxList));
%     if i==9
%         fprintf('process %d out of %d\n', i, numel(s.VoxelIdxList));
%     end
    seed_id = seed_proc_order(i);
%     seed_id = 1347;                           % debug mode
    yxz = s.VoxelIdxList{seed_id};
    ids = idMap_current(yxz);
    yxz = yxz(ids==seed_id);
    if iscell(idMap)
        idMap{2} = idMap_current;
    else
        idMap = idMap_current;
    end
    [newLabel, comMaps, fgReDo] = m_refineOneRegion_with_seed(seed_id, yxz, vid,...
        vid_stb, idMap, eig2d, eig3d, OrSt, q);
    if fgReDo
        fgReDoCnt = fgReDoCnt + 1;
        %fprintf('Initial fg are too small \n');
    end
    %% save cell locations
    % modify: change threshold
    valid_newLabel = newLabel(:)>0;
    loc_cells{i} = cell(2,1);
    loc_cells{i}{1} = [comMaps.linerInd(valid_newLabel), newLabel(valid_newLabel)];
    loc_cells{i}{2} = round(quantile(comMaps.vidComp(valid_newLabel),0.5));
    
%     loc_cells{i}{2} = comMaps.pickedThreshold;
    
    %% do a simple update to the idMap
    idMap_current(yxz) = 0;
    idMap_current(loc_cells{i}{1}(:,1)) = seed_id;

    %% a further check for cells in 2nd round, remove those dim and small ones
    if ~isempty(test_ids) % this is appended cells
        vals = vid_current(loc_cells{i}{1}(:,1));
        ids = unique(loc_cells{i}{1}(:,2));
        for j=1:length(ids)
            tmp_cell = find(loc_cells{i}{1}(:,2)==ids(j));
            if length(tmp_cell) < minSz ||  mean(vals(tmp_cell)) < min_level % smaller than min size or 5% minimal
                loc_cells{i}{1}(tmp_cell,1) = 0;
            end
        end
        loc_cells{i}{1}(loc_cells{i}{1}(:,1)==0,:) = [];
    end
    %% write data into images for error check
%     if ~isempty(save_folder)
%         writeRefineCell(comMaps.vidComp, newLabel, comMaps.regComp, i, save_folder);
%     end
end
% label the new ID Map
newIdMap = zeros(size(idMap_current));
thresholdMap = zeros(size(idMap_current));
regCnt = 0;
for i=1:numel(loc_cells)
    if ~isempty(loc_cells{i}{1})
        cur_locs = loc_cells{i}{1}(:,1);
        cur_labels = loc_cells{i}{1}(:,2);
        newIdMap(cur_locs) = cur_labels + regCnt;
        thresholdMap(cur_locs) = loc_cells{i}{2};
        regCnt = regCnt + max(cur_labels);
    end
end

if ~isempty(tif_id)
    [H1,W1,D1] = size(idMap_init);
    overlay_im = zeros([H1,W1,3,D1], 'uint8');
    %vid_current = scale_image(vid_current, 0, 255);
    if ~isempty(test_ids)
        cell_1st = ~ismember(idMap_init, test_ids) & idMap_init>0;
        for i=1:D1
            overlay_im(:,:,1,i) = uint8(cell_1st(:,:,i) * 100);
            overlay_im(:,:,2,i) = uint8(vid_current(:,:,i)*2);
            overlay_im(:,:,3,i) = uint8((newIdMap(:,:,i)>0) * 100);
        end
    else
        for i=1:D1
            overlay_im(:,:,1,i) = uint8((newIdMap(:,:,i)>0) * 100);
            overlay_im(:,:,2,i) = uint8(vid_current(:,:,i)*2);
        end
    end
    
    svf = 'C:\Users\Congchao\Desktop\cell_detection_samples\crop_embryo_data_500x500x30x40\segmentation_from_synQuant\';
    tifwrite(overlay_im, fullfile(svf, ['overlay_synQuant_', num2str(tif_id)]));
end
%% repeat the refinement if there is remaining seed
% practically not applicable due to
% 1. large number of false positive seeds (e.g. the boundary of cells)
% 2. fg detection bias to the nearby regions (that's why these seeds were left)
% idMap{2}(newIdMap>0) = 0;
% idMap2 = region_sanity_check(idMap{2}, q.minSeedSize, true);
% s = regionprops3(idMap2, {'VoxelIdxList'});
% idMap2(idMap2>0) = idMap2(idMap2>0) + regCnt;
% idMap{2} = idMap2 + newIdMap;
% loc_cells = cell(numel(s.VoxelIdxList),1);
% for i=regCnt+1 : regCnt+numel(s.VoxelIdxList)%[7 9 13 15 17 18 19 20 21 22]%
%     yxz = s.VoxelIdxList{i-regCnt};
% %     if mod(i, 10) == 0
% %         fprintf('process %d out of %d\n', i, numel(s.VoxelIdxList));
% %     end
%     [newLabel, comMaps, fgReDo] = refineOneRegion_with_seed(i, yxz, vid,...
%         vid_stb, idMap, eig2d, eig3d, OrSt, q);
%     if fgReDo
%         fgReDoCnt = fgReDoCnt + 1;
%         %fprintf('Initial fg are too small \n');
%     end
%     %% save cell locations
%     valid_newLabel = newLabel(:)>0;
%     loc_cells{i-regCnt} = cell(2,1);
%     loc_cells{i-regCnt}{1} = [comMaps.linerInd(valid_newLabel), newLabel(valid_newLabel)];
%     loc_cells{i-regCnt}{2} = comMaps.pickedThreshold;
%     %% write data into images for error check
%     if ~isempty(save_folder)
%         writeRefineCell(comMaps.vidComp, newLabel, comMaps.regComp, i, save_folder);
%     end
% end
% fprintf('Totally %d cells has small initial foreground\n', fgReDoCnt);
% % label the new ID Map
% for i=1:numel(loc_cells)
%     if ~isempty(loc_cells{i}{1})
%         cur_locs = loc_cells{i}{1}(:,1);
%         cur_labels = loc_cells{i}{1}(:,2);
%         newIdMap(cur_locs) = cur_labels + regCnt;
%         thresholdMap(cur_locs) = loc_cells{i}{2};
%         regCnt = regCnt + max(cur_labels);
%     end
% end

end