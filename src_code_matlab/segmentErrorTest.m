function [p_k_vec, minmax_ids, voxIdx_min_id, voxIdx_max_id] = segmentErrorTest(parent_id, kid_id, g, movieInfo, refine_res, embryo_vid)
% two kinds of error we are testing
% 1: large size difference
% 2: large location difference

%ovSz = movieInfo.ovSize{parent_id}(movieInfo.nei{parent_id}==kid_id);
ovSz = length(intersect(movieInfo.voxIdx{parent_id}, movieInfo.voxIdx{kid_id}));
p_Sz = size(movieInfo.vox{parent_id}, 1);
k_Sz = size(movieInfo.vox{kid_id}, 1);

if p_Sz > k_Sz
    min_id = kid_id;
    max_id = parent_id;
    size_ratio = p_Sz / k_Sz;
else
    max_id = kid_id;
    min_id = parent_id;
    size_ratio = k_Sz / p_Sz;
end
minmax_ids = [min_id, max_id];
p_k_vec = [];
if k_Sz==0 || p_Sz==0
    voxIdx_min_id = [];
    voxIdx_max_id = [];
    return;
end
% if the breaking reason is "large size difference"
if size_ratio > 2
    % if this is a highly overlapped small region
    if ovSz/min(p_Sz, k_Sz)>0.5
        [~, minDistance] = ovDistanceRegion(movieInfo.vox{max_id}, movieInfo.vox{min_id});
        minDistance = overlap2cost(minDistance, movieInfo.ovGamma);
        if minDistance < min(g.c_ex, g.c_en)
            p_k_vec = [max_id, min_id];
            voxIdx_min_id = [];
            voxIdx_max_id = [];
            return;
        end
    end
    
end
% if the previous criteria is not met, breaking reason is location
% difference. We will region grow to see
% shift = [3 3 1];
% punish = 0.1;
% cost_design = [1, 2];
% connect = 6;
% %seed_idxs = intersect(movieInfo.voxIdx{max_id}, movieInfo.voxIdx{min_id});
% %[sy, sx, sz] = ind2sub(size(refine_res{1}), seed_idxs);
% 
% yxz = union(movieInfo.voxIdx{max_id}, movieInfo.voxIdx{min_id});
% [uy, ux, uz] = ind2sub(size(refine_res{1}), yxz);
% 
% [dumbReg, correspIdx, ~, loc_org_xyz] = crop3D(embryo_vid{1}, yxz, shift);
% %seed_vox = [sx, sy, sz] - (loc_org_xyz-1);
% %seed_idxs = sub2ind(size(dumbReg), seed_vox(:,2), seed_vox(:,1), seed_vox(:,3));
% union_vox = [ux, uy, uz] - (loc_org_xyz-1);
% union_idxs = sub2ind(size(dumbReg), union_vox(:,2), union_vox(:,1), union_vox(:,3));
% ids = [max_id, min_id];
% frs = [movieInfo.frames(max_id), movieInfo.frames(min_id)];
% scoreMap = cell(3, 1);
% fMap = cell(3, 1);
% sMap = cell(3, 1);
% tMap = cell(3, 1);
% for i=1:length(frs)
%     vidComp = crop3D(embryo_vid{frs(i)}, yxz, shift);
%     idComp = crop3D(refine_res{frs(i)}, yxz, shift);
%     sMap{i} = idComp == refine_res{frs(i)}(movieInfo.voxIdx{ids(i)}(1));
%     
% %     fMap{i} = idComp > 0;
% %     sMap{i} = zeros(size(idComp));
% %     sMap{i}(seed_idxs) = 1;
% %     tMap{i} = ~sMap{i} & idComp>0;
%     fMap{i} = zeros(size(idComp));
%     fMap{i}(union_idxs) = 1;
%     tMap{i} = ~fMap{i} | (~sMap{i} & idComp>0);
%     
%     eig2dComp = principalCv2d(vidComp, idComp, 3);%
%     eig2dComp(eig2dComp<0) = 0;
%     scoreMap{i} = scale_image(eig2dComp, 1e-3,1);% remove 0s
% end
% 
% [dat_in, src, sink] = graphCut_4d(scoreMap, fMap, sMap, tMap, punish, connect, cost_design);
% 
% G = digraph(dat_in(:,1),dat_in(:,2),dat_in(:,3));
% [~,~,cs,~] = maxflow(G, src, sink); % cs: fg, ct: bg
% cs = cs(cs<src);
% % [y, x, z] = ind2sub(size(sMap{1}), cs(cs <= numel(sMap{1})));
% % max_vox = [x, y, z] + (loc_org_xyz-1);
% [y, x, z] = ind2sub(size(sMap{1}), cs(cs > numel(sMap{1})) - numel(sMap{1}));
% min_vox = [x, y, z] + ((loc_org_xyz-1) - ...
%     movieInfo.drift(movieInfo.frames(min_id),:));

[voxIdx_min_id, ~, voxIdx_max_id] = growHighOvReg(max_id, min_id, movieInfo, refine_res, embryo_vid);
if isempty(voxIdx_min_id) || isempty(voxIdx_min_id)
    maxDistance = inf;
else
    [y, x, z] = ind2sub(size(refine_res{1}), voxIdx_min_id);
    min_vox = [x, y, z] - movieInfo.drift(movieInfo.frames(min_id),:);

    [y, x, z] = ind2sub(size(refine_res{1}), voxIdx_max_id);
    max_vox = [x, y, z] - movieInfo.drift(movieInfo.frames(max_id),:);


    [maxDistance, ~] = ovDistanceRegion(max_vox, min_vox);
    maxDistance = overlap2cost(maxDistance, movieInfo.ovGamma);
end
if maxDistance >= abs(g.observationCost)
%     movieInfo.vox{min_id} = min_vox;
%     min_vox = [x, y, z] + (loc_org_xyz-1);
%     voxIdx_min_id = sub2ind(size(refine_res{1}), ...
%         min_vox(:,2), min_vox(:,1), min_vox(:,3));
% else
    voxIdx_min_id = [];
    voxIdx_max_id = [];
end
