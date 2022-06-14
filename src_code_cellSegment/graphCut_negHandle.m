function [dat_in, src_node, sink_node] = graphCut_negHandle(vid,fMap,sMap,...
    tMap, connect, cost_design, bg2sink)
% build a graph for graph cut, note that we only test the gaps inside the
% one region, if there is some obvious gaps, we will segment it; for those
% segmented regions, we will not grow the gaps between themshift
% INPUT:
% vid: 3D (TODO: 2D) matrix indicating the similarity of voxels, it can be
% principal curvature or the gradient. if it is principal curvature,
% there may be negative values.
% fMap: the valid foreground map for segmenting
% sMap and tMap: the voxel that definitely belongs to src or sink
% connect: 6, 26 (TODO: 10(8+2)) connection for edges in the graph
% cost_design: 1 means average; 2 means sqrt;
% OUTPUT:
% dat_in: the graph with n+1 indicating src and n+2 indicating sink; for a voxel at
% location ind, its node index should be ind+1

% contact: ccwang@vt.edu

if nargin == 6
    % if we did not link background to sink, no need to consider such
    % nodes, otherwise bg2sink is true
    bg2sink = true;
end
vox_num = numel(vid);
% 
valid_vox = find(fMap);
valid_vox_num = length(valid_vox);
adjMap = zeros(size(vid)); % map of neighbors

dat_in = zeros(valid_vox_num*connect*2 + vox_num, 3);
k_dat = 0;
% connections inside principal map or gradient map
for i=1:valid_vox_num
    ind = valid_vox(i);
    p1 = vid(ind);
    [~, ~, nei_ids] = neighbours(ind, vid, connect);
    if ~bg2sink
        nei_ids = nei_ids(fMap(nei_ids));
    end
    
    p2 = vid(nei_ids);
    if cost_design(1)==1
        costs = (2./(p1+p2)).^cost_design(2);
    elseif cost_design(1)==2
        p2(p2==0) = p1; % or inf?
        costs = (1./sqrt(p1.*p2)).^cost_design(2);
%         costs = 2./(p1.*p2);
    end
    % negative principal curvature linked to fg for sure
    % but if this score is based on gradient, no negative value exists
%     if p1 < 0 
%         costs(sMap(nei_ids)) = inf;
%         costs(tMap(nei_ids)) = 0;
%         
%         costs(p2<0) = inf;
%     end
%     %costs(p1+p2<0) = 0;
%     if sMap(ind)
%         costs(p2<0 & tMap(nei_ids)) = inf;
%     end
%     if tMap(ind)
%         costs(p2<0) = 0;
%     end
%     costs(costs < 0) = 0;
    
%     if sum(isinf(costs) & ((tMap(nei_ids) & sMap(ind)) | (sMap(nei_ids) & tMap(ind))))>0
%         error('inf flow');
%     end

    k_dat = k_dat + length(nei_ids);
    if k_dat > size(dat_in, 1)
        fprintf('memory allocate error\n');
    end
    dat_in(k_dat-length(nei_ids)+1:k_dat,:) = [ind+nei_ids*0, nei_ids, costs];

    k_dat = k_dat + length(nei_ids);
    if k_dat > size(dat_in, 1)
        fprintf('memory allocate error\n');
    end
    dat_in(k_dat-length(nei_ids)+1:k_dat,:) = [nei_ids, ind+nei_ids*0, costs];

    adjMap(nei_ids) = 1;
end
% connections to src or sink
src_node = vox_num+1;
srcIds = find(sMap);

k_dat = k_dat + length(srcIds);
dat_in(k_dat-length(srcIds)+1:k_dat,:) = [src_node + srcIds*0, srcIds, inf + srcIds*0];

% for i=1:srcIds
%     k_dat = k_dat + 1;
%     cost = inf; 
%     dat_in(k_dat,:) = [src_node, srcIds(i), cost];
% end
    
sink_node = vox_num+2;
sinkIds = find(adjMap & tMap);
if isempty(sinkIds)
    sink_node = nan;
end
k_dat = k_dat + length(sinkIds);
dat_in(k_dat-length(sinkIds)+1:k_dat,:) = [sinkIds, sink_node + sinkIds*0, inf + sinkIds*0];

% for i=1:sinkIds
%     k_dat = k_dat + 1;
%     cost = inf; 
%     dat_in(1:k_dat,:) = [sinkIds(i), sink_node, cost];
% end

dat_in = dat_in(1:k_dat,:);
end