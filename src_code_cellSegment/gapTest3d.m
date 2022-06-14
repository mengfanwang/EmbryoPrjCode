function [newLabel, regTest] = gapTest3d(L, reg_id, vidComp, idComp...
    , scoreMap,...
    fmap, connect, cost_design, edgeConnect, p_thres, OrSt)
% test if the gap is significant enough in a 3D manner
%INPUT:
% L: the label of regions segmented by principal curvature
% n: number of regions segmented by principal curvature
% vidComp: the original image data (Comp means alreay cropped as all other inputs)
% idComp: the label of regions before segmentation by principal curvature
% eigComp: eigen values of all foreground voxels
% fmap: map of valid eigen value, (eigen value > threshold)
% connect: 6, 26 (TODO: 10(8+2)) connection for connecting edges in the graph
% cost_design: how to design the cost in the Graph-Cut
% edgeConnect: in edgeTest, the way to determine neighbors, (default 124 as 5*5*5)
% p_thres: pvalue threshold to define significant gaps
% OrSt: structure from order statistics, containing data noise variance, mu
% and sigma from order statistics (maybe empty)
%OUTPUT:
% newLabel: new segmentation which merges some segments by removing all
% fake gaps
% regTest: 1=>pixels of edge between two cells, 2=>pixels of cell for testing,
% 3=>pixels of the other cell for testing

fmapIn = fmap;
other_id = idComp~=reg_id & idComp>0 & fmapIn;
append_id = nan;
if ~isempty(find(other_id, 1))
    append_id = max(L(:)) + 1;
    L(other_id) = append_id;
end
[h, w, z] = size(fmapIn);
tmpMap = false(h,w);
for i=1:2 % the two cycles of boundary are background
    tmpMap(i,:) = true;
    tmpMap(end-i+1,:) = true;
    tmpMap(:,i) = true;
    tmpMap(:, end-i+1) = true;
end
bndIdx = find(tmpMap);
for i=1:z
    cur_bidx = bndIdx + (i-1)*w*h;
    cur_bidx = cur_bidx(idComp(cur_bidx)==0);
    fmapIn(cur_bidx) = 0;
end
newLabel = regionGrow(L, scoreMap, fmapIn, connect, cost_design);
if ~isnan(append_id)
    newLabel(newLabel == append_id) = 0;
end

% merge over-segmented regions
inLabel = newLabel;

[newLabel, regTest] = edgeTest3d(vidComp, inLabel, fmap, edgeConnect, p_thres, OrSt);


end