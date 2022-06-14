function [newLabel, edgeLabel, testFlag] = segmentCurrentRegion(comMaps, reg_id, q, OrSt)
% test the gaps in the region to see if it can be segmented

[L, ~] = bwlabeln(comMaps.newIdComp, q.neiMap);                             % remove small objects
S = regionprops3(L, 'Volume');
inval_ids = find([S.Volume] <= q.minSeedSize);
if ~isempty(inval_ids)
    invalidMap = ismember(L, inval_ids);
    L(invalidMap) = 0;
    L = rearrange_id(L);
end
% for small region, do not grow them, merge them to principal
% curvature maps
v_l = regionprops3(L,'VoxelIdxList');
if numel(v_l.VoxelIdxList) > 1
    [newLabel, edgeLabel] = gapTest3dV2(L, comMaps, reg_id, q, OrSt);
    testFlag = true;
else
    testFlag = false;
    edgeLabel = comMaps.regComp*0;
    newLabel = double(comMaps.regComp);
end