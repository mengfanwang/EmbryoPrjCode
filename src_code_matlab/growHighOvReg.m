function [min_id_voxIdx, loc_org_xyz, max_id_voxIdx] = growHighOvReg(...
    max_id, min_ids, movieInfo, refine_res, embryo_vid, ...
    eigMaps, q)
% if the previous criteria is not met, breaking reason is location
% difference. We will region grow to see

% shift = [3 3 1];
punish = 1;
% cost_design = [1, 2];
% connect = 6;
% smooth_sc = [5 3];
seed_idxs = cat(1, movieInfo.voxIdx{min_ids});
[sy, sx, sz] = ind2sub(size(refine_res{1}), seed_idxs);

yxz = union(movieInfo.voxIdx{max_id}, cat(1, movieInfo.voxIdx{min_ids}));
[uy, ux, uz] = ind2sub(size(refine_res{1}), yxz);

[dumbReg, correspIdx, ~, loc_org_xyz] = crop3D(refine_res{1}, yxz, q.shift);
seed_vox = [sx, sy, sz] - (loc_org_xyz-1);
seed_idxs = sub2ind(size(dumbReg), seed_vox(:,2), ...
    seed_vox(:,1), seed_vox(:,3));

union_vox = [ux, uy, uz] - (loc_org_xyz-1);
union_idxs = sub2ind(size(dumbReg), union_vox(:,2), union_vox(:,1), union_vox(:,3));
ids = [max_id, min_ids(1)];
frs = movieInfo.frames(ids);%[movieInfo.frames(max_id), movieInfo.frames(min_ids)];
scoreMap = cell(3, 1);
fMap = cell(3, 1);
sMap = cell(3, 1);
tMap = cell(3, 1);
for i=1:length(frs)
    vidComp = crop3D(embryo_vid{frs(i)}, yxz, q.shift);
    idComp = crop3D(refine_res{frs(i)}, yxz, q.shift);
    if i==1
        sMap{i} = idComp == refine_res{frs(i)}(movieInfo.voxIdx{ids(i)}(1));
    else
        sMap{i} = zeros(size(idComp));
        sMap{i}(seed_idxs) = 1;
    end

    score2dMap = crop3D(eigMaps{frs(i)}{1}, yxz, q.shift);
    score2dMap(score2dMap<0) = 0;
    score3dMap = crop3D(eigMaps{frs(i)}{2}, yxz, q.shift);
    score3dMap(score3dMap<0) = 0;
    
    fMap{i} = zeros(size(idComp));
    fMap{i}(union_idxs) = 1;
    
    tMap{i} = ~fMap{i} | (~sMap{i} & idComp>0);
    
    scoreMap{i} = scale_image(score2dMap, 1e-3,1) + ...
        scale_image(score3dMap, 1e-3,1);% remove 0s
end

[dat_in, src, sink] = graphCut_4d(scoreMap, fMap, sMap, tMap, ...
    punish, connect, cost_design);

G = digraph(dat_in(:,1),dat_in(:,2),dat_in(:,3));
[~,~,cs,~] = maxflow(G, src, sink); % cs: fg, ct: bg
cs = cs(cs<src);
min_id_voxIdx = correspIdx(cs(cs > numel(sMap{1})) - numel(sMap{1}));
max_id_voxIdx = correspIdx(cs(cs <= numel(sMap{1})));
if false
    zz = zeros(size(fMap{1}));
    zz(cs(cs > numel(sMap{1})) - numel(sMap{1})) = 1;
    zz(cs(cs <= numel(sMap{1})))  =zz(cs(cs <= numel(sMap{1})))  + 2;
    zzshow(label2rgb3d(zz));
end
% [y, x, z] = ind2sub(size(sMap{1}), cs(cs > numel(sMap{1})) - numel(sMap{1}));
% min_vox = [x, y, z] + ((loc_org_xyz-1) - ...
%     movieInfo.drift(movieInfo.frames(min_ids),:));
% [y, x, z] = ind2sub(size(sMap{1}), cs(cs <= numel(sMap{1})) );
% max_vox = [x, y, z] + ((loc_org_xyz-1) - ...
%     movieInfo.drift(movieInfo.frames(max_id),:));

end