function cell_2ndIter = m_Synquant4Embryo_2iter(imgIn, cellIdMap_1stIter, posEigMap)
% re-detect cells using synQuant given results from 1st iteration

% posEigMap: map of positive principal curvature
% contact: ccwant@vt.edu 02/04/2020
%% get the cell information
s_1st = regionprops3(cellIdMap_1stIter, {'VoxelIdxList','Volume','SurfaceArea'});
cell_sz = [s_1st.Volume];
cell_levels = cellfun(@(x) mean(imgIn(x)), s_1st.VoxelIdxList);
min_level = quantile(cell_levels,0.05);
minSz = 50; %min(cell_sz); %quantile(cell_sz,0.01); not reasonable. should use seed size
maxSz = quantile(cell_sz,0.95);
imgIn_org = imgIn;
%% mask out all the detected cells and re-do everything
imgIn = imgIn_org;
cell_map_dilate = imdilate(cellIdMap_1stIter>0, strel('disk',1));
cell_map_dilate = imdilate(cell_map_dilate, strel('sphere',1));
mask_out = posEigMap | cell_map_dilate;%

seedMap = ~mask_out; 
seedMap_dilate = imdilate(seedMap, strel('disk',1));% gap < 2pixels will be adjacent

synId_append = bwlabeln(seedMap_dilate, 26);
synId_append(mask_out) = 0;

s_2nd = regionprops3(synId_append, {'VoxelIdxList','Centroid'});
cnt = 0;
cell_2ndIter = zeros(size(imgIn),'uint32');
cell_level = zeros(numel(s_2nd.VoxelIdxList), 2);
for i=1:numel(s_2nd.VoxelIdxList)
    tmp_level = mean(imgIn(s_2nd.VoxelIdxList{i}));
    %center_xyz = s_2nd.Centroid(i,:);
    %center_ind = get_index_from_center(center_xyz([2 1 3]), imgIn);
    %valid_center = ismember(center_ind,s_2nd.VoxelIdxList{i}); % center should inside region

    %if sum(valid_center)==length(valid_center) && ...
    if length(s_2nd.VoxelIdxList{i}) >= minSz && ...
            length(s_2nd.VoxelIdxList{i}) <= maxSz && tmp_level > min_level
        cnt = cnt + 1;
        cell_2ndIter(s_2nd.VoxelIdxList{i}) = cnt;
        cell_level(cnt,1) = tmp_level;
        cell_level(cnt,2) = length(s_2nd.VoxelIdxList{i});
    end
end

end





