function regVoxIdx = bisectRegion_gapguided(root_voxIdx, root_frame, ...
    base_voxIdx, vidMap, refine_res, eigMaps, validGapMaps, q)
% test if 
% 1. there is a gap defined by principal curvature
% 2. the gap separate the root region into two regions consistent with
% the TWO regions defined by base_voxIdx

% If the two criterions are both satisfied, we separate the root region
% into two regions, otherwise return empty elements

if numel(base_voxIdx) ~= 2 || ~isempty(find(cellfun(@isempty, base_voxIdx),1))% split into two regions
    regVoxIdx = cell(2,1);
    return;
end

yxz = union(root_voxIdx, cat(1, base_voxIdx{:}));
% crop regions
[idComp, linerInd, ~, loc_org_xyz] = crop3D(refine_res{root_frame}, ...
    yxz, q.shift);

% 2d principal curvature
eig2dComp = crop3D(eigMaps{root_frame}{1}, yxz, q.shift);
possibleGaps2d = eig2dComp>0;
eig2dComp(eig2dComp<0) = 0;
min_pv = 0;%min(eigMaps{root_frame}{1}(root_voxIdx));
max_pv = max(eigMaps{root_frame}{1}(root_voxIdx));
eig2dComp = scale_image(eig2dComp, 1e-3,1, min_pv, max_pv);
eig2dComp(isnan(eig2dComp)) = 0;


% 3d principal curvature
eig3dComp = crop3D(eigMaps{root_frame}{2}, yxz, q.shift);
possibleGaps3d = eig3dComp>0;
eig3dComp(eig3dComp<0) = 0;
min_pv = 0;
max_pv = max(eigMaps{root_frame}{2}(root_voxIdx));
eig3dComp(isnan(eig3dComp)) = 0;
eig3dComp = scale_image(eig3dComp, 1e-3,1, min_pv, max_pv);
if ~isempty(validGapMaps)
    gapComp = crop3D(validGapMaps{root_frame}, yxz, q.shift);
    possibleGaps2d = possibleGaps2d & gapComp;
    possibleGaps3d = possibleGaps3d & gapComp;
end
% [vidComp, linerInd, ~, loc_org_xyz] = crop3D(vidMap{root_frame}, yxz, ...
%     q.shift);

% build valid label map
real_id = refine_res{root_frame}(root_voxIdx(1));
base_idx_local = coordinate_transfer(base_voxIdx, ...
    size(refine_res{root_frame}), loc_org_xyz, size(idComp));
[label_map_binary, vflag] = binary_seeds_create(idComp==real_id, ...
    possibleGaps3d, [], base_idx_local);%, q.minSeedSize
if ~vflag % if 3d cannot detect gaps, we change to 2d gaps
    [label_map_binary, vflag] = binary_seeds_create(idComp==real_id, ...
        [], possibleGaps2d, base_idx_local);
    if ~vflag % indeed we can still increase treshold, but do we need?
        regVoxIdx = cell(2,1);
        return;
    end
end
% [label_map,n] = bwlabeln(idComp == real_id & ~possibleGaps, 6);
% % get the index of splitted regions
% base_idx_local = cell(numel(base_voxIdx), 1); % label current region
% base_labels = nan(n, numel(base_idx_local));
% for j=1:numel(base_idx_local)
%     [y, x, z] = ind2sub(size(vidMap{root_frame}), base_voxIdx{j});
%     base_idx_local{j} = [x, y, z] - (loc_org_xyz-1);
%     base_idx_local{j} = sub2ind(size(vidComp), base_idx_local{j}(:,2), ...
%         base_idx_local{j}(:,1), base_idx_local{j}(:,3));
%     out_freq = frequency_cnt(label_map(base_idx_local{j}));
%     for k=1:size(out_freq,1)
%         if out_freq(k,1) > 0
%             base_labels(out_freq(k,1), j) = out_freq(k,2);
%         end
%     end
% end
% label_map_binary = zeros(size(label_map));
% all_seed_used = false(numel(base_idx_local), 1);
% for i = 1:n
%     [v, od] = nanmax(base_labels(i,:));
%     if ~isnan(v)
%         label_map_binary(label_map == i) = od;
%         all_seed_used(od) = true;
%     end
% end
% % if there is a based region that has no correponding region in root region
% if ~isempty(find(~all_seed_used,1))
%     regVoxIdx = cell(2,1); % is this reasonable???
%     return;
% end


%% region grow the regions from label_map_binary
% bg2sink = false;
% newLabel = regionGrow(label_map_binary, eig2dComp + eig3dComp, ...
%     idComp == real_id, q.growConnectInRefine, q.cost_design, bg2sink);
% regVoxIdx{1} = linerInd(newLabel==1);
% regVoxIdx{2} = linerInd(newLabel==2);

fMap = idComp == real_id;
new_l = zeros(size(fMap));
scoreMap = eig2dComp + eig3dComp;
regVoxIdx = cell(numel(base_voxIdx), 1);
for i=1:numel(base_voxIdx)
    sMap = label_map_binary == i;
    tMap = label_map_binary>0;
    tMap(sMap) = false;
    [dat_in, src, sink] = graphCut_negHandle_mat(scoreMap, fMap, sMap, tMap, ...
        q.growConnectInRefine, q.cost_design, false);
    % if src/sink cannot be linked in the foreground, return
    if isempty(find(dat_in(:,1) == src,1)) || ...
            isempty(find(dat_in(:,2) == sink,1))
        regVoxIdx = cell(2,1);
        return;
    end
    if ~isempty(find(isnan(dat_in),1))
        keyboard;
    end
    G = digraph(dat_in(:,1),dat_in(:,2),dat_in(:,3));
    [~,~,cs,~] = maxflow(G, src, sink); % cs: fg, ct: bg
    reg_locs = cs(cs<numel(fMap));
    regVoxIdx{i} = linerInd(reg_locs);
    
    new_l(reg_locs) = i;
%     if all_seed_used(i) % if only exist valid seed
%         sMap = label_map_binary == i;
%         tMap = label_map_binary>0;
%         tMap(sMap) = false;
%         if ~isempty(find(tMap,1)) % if there is another seed
%             [dat_in, src, sink] = graphCut_negHandle(scoreMap, fMap, sMap, tMap, ...
%                 q.growConnectInRefine, q.cost_design, false);
%         else% if there is no other seeds, sink link to background
%             tmp_fgMap = idComp == real_id;
%             tMap(imdiate(tmp_fgMap, strel('disk',1))-tmp_fgMap>0) = true;
%             [dat_in, src, sink] = graphCut_negHandle(scoreMap, fMap, sMap, tMap, ...
%                 q.growConnectInRefine, q.cost_design, true);
%         end
%         if ~isempty(find(isnan(dat_in),1))
%             keyboard;
%         end
%         G = digraph(dat_in(:,1),dat_in(:,2),dat_in(:,3));
%         [~,~,cs,~] = maxflow(G, src, sink); % cs: fg, ct: bg
%         reg_locs = cs(cs<numel(fMap));
%         regVoxIdx{i} = linerInd(reg_locs);
%     end
end
end