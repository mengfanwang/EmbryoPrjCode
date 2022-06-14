function [newLabel, regComp] = edgeTest3dV2(vid, label, fgMap, connect, p_thres, OrSt)
% test if the edge is indeed weaker than surrounding area or simply becuase
% of noise
% INPUT:
% vid: original 3d image data
% label: current segmentation results
% fgMap: current foreground map, which contains all boundary area whose
% principla curvature is also large
% connect: the way to determine neighbors
% p_thres: the pvalue threshold to define significant edge
% OrSt: mu and sigma matrix for order statistics analysis
% OUTPUT:
% newLabel: new segmentation which merges some segments by removing all fake edges
% regTest: 1=>pixels of edge between two cells, 2=>pixels of cell for testing,
% 3=>pixels of the other cell for testing

% contact: ccwang@vt.edu, 12/03/2019

if nargin < 6
    OrSt = [];
end
% way 1: process the data slice by slice
[h,w,z] = size(vid);
% dilate the data (x,y direction) to remove boundary testing
dilate_sc = 4; % maximum neighbor size
sph = strel('sphere', dilate_sc);
fgMap = imdilate(fgMap, strel(sph.Neighborhood(:,:,dilate_sc+1)));
vid_di = nan(h+dilate_sc*2, w+dilate_sc*2, z+dilate_sc*2);
vid_di(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc) = vid;
label_di = zeros(size(vid_di));
label_di(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc) = label;
fgMap_di = zeros(size(vid_di));
fgMap_di(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc) = fgMap;
varMap = zeros(size(vid_di));
varMap(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc) = OrSt.stbVarCropMap;    % use stb varMap no options!
[h,w,z] = size(vid_di);

regComp = zeros(h,w,z);
fg_locs = find(fgMap_di>0);                         

% test each pair to see if we can merge
merged_regions = cell(100,1);
bk_regions = cell(100,1);
bk_cnt = 0;
mr_cnt = 0;
min_r = 2;
bk_perT = zeros(dilate_sc, 1);
for i=min_r:dilate_sc                                                       % neighbor radius
    gap = 2*i+2; % five circles of neighbors

    neighbors = zeros(2*i+1);
    nCnt = 0;
    for j=-i:i                                                              % 2D [-i:i,-i:i] neighbor
        nCnt = nCnt + 1;
        neighbors(nCnt,:) = j*h-i:j*h+i;
    end
    neighbors = neighbors(:)';
    neighbors(neighbors==0) = [];
    slabel = label_di;
    sComp = zeros(h,w,z);
    
    nei_mat = repmat(fg_locs, 1, length(neighbors)) +...                    % build neighbor mat for each pixel in foreground
        repmat(neighbors, length(fg_locs), 1);
    fg_neighbors = slabel(nei_mat);
    
    nei_segments = nan(size(fg_neighbors,1),2);
    for j=1:size(fg_neighbors,1)
        cur_neis = fg_neighbors(j,:);
        cur_neis(cur_neis == 0) = [];
        uni_neis = sort(unique(cur_neis));
        if length(uni_neis)==2                                              % detect segment (different label of a pair)
            nei_segments(j,:) = uni_neis;                                   
        end
    end
    if i==min_r % for regions with no overlap, not considering them         % remove duplicated region
        pair_candidates = unique(nei_segments,'rows');
        pair_candidates(isnan(pair_candidates(:,1)),:) = [];                % pairs of labels with gap
    end
    if isempty(pair_candidates)
        break;
    end
    for j=1:size(pair_candidates,1)
        duplicate_flag = false;                                                         
        for bb = 1:bk_cnt                                                   % remove duplicated pairs
            if sum(abs(pair_candidates(j,:)-bk_regions{bb}))==0
                duplicate_flag = true;
                break;
            end
        end
        if (pair_candidates(j,1) ~= pair_candidates(j,2)) && ~duplicate_flag
            edge_all = nei_segments(:,1) == pair_candidates(j,1) & ...      % logical array of edge list
                nei_segments(:,2) == pair_candidates(j,2);
            edge_locs = fg_locs(edge_all);% & label_di(fg_locs) == 0
            edge_locs = gapRefine(slabel, pair_candidates(j,:), edge_locs);
            if length(edge_locs) < 10 % means two regions are too far away:real gap
                p1 = nan;
                p2 = nan;
            else
                edge_pixels = vid_di(edge_locs);
                edge_vals = slabel(edge_locs);
                slabel(edge_locs) = 0; % remove edges from labels
                % use the pixels adjacent to the boundaries
                neighbor_cells = neighbor_pixels(edge_locs, size(slabel), gap);
                neighbors_locs = cat(1, neighbor_cells{3:gap});             % jump 2 circles
                invalid_nei_locs = cat(1, neighbor_cells{1:2});
                invalid_nei = vid_di(invalid_nei_locs(...
                    slabel(invalid_nei_locs)==pair_candidates(j,1) | ...
                    slabel(invalid_nei_locs)==pair_candidates(j,2)));
                
                r1_locs = neighbors_locs(slabel(neighbors_locs)==pair_candidates(j,1));
                r1 = vid_di(r1_locs);                                       % pixels in region 1
                r2_locs = neighbors_locs(slabel(neighbors_locs)==pair_candidates(j,2));
                r2 = vid_di(r2_locs);
                
                slabel(edge_locs) = edge_vals;
                sComp(edge_locs) = 1;
                sComp(neighbors_locs(slabel(neighbors_locs)==pair_candidates(j,1))) = 2;
                sComp(neighbors_locs(slabel(neighbors_locs)==pair_candidates(j,2))) = 3;
                % now we are testing edge_pixels and pixels in area r1, r2
                %[~, p1] = ttest2(edge_pixels, [r1;r2]);
                if isfield(OrSt,'stbVarCropMap')
                    % we need to estimate variance every time
                    curVar = nanmean(varMap(cat(1,r1_locs,r2_locs,edge_locs)));
                else
                    curVar = OrSt.noiseVar;
                end
                if strcmp(OrSt.gapTestWay, 'ttest') %~isfield(OrSt,'mu') || ~isfield(OrSt,'sigma')
                    [~, p1] = ttest2(edge_pixels, r1);
                    [~, p2] = ttest2(edge_pixels, r2);
                elseif strcmp(OrSt.gapTestWay, 'loopUpTbl') %~isempty(OrSt.mu) && ~isempty(OrSt.sigma)
                    M = max(min(length(r1), size(OrSt.mu,1)), 10);
                    N = max(min(length(edge_pixels), size(OrSt.mu,1)), 10);
                    sum_st = mean(r1)-mean(edge_pixels);
                    zscore = (sum_st - OrSt.mu(M,N)*sqrt(curVar))...
                        / (OrSt.sigma(M,N)*sqrt(curVar));
                    p1 = 1-normcdf(zscore);
                    M = max(min(length(r2), size(OrSt.mu,1)), 10);
                    sum_st = mean(r2)-mean(edge_pixels);
                    zscore = (sum_st - OrSt.mu(M,N)*sqrt(curVar))...
                        / (OrSt.sigma(M,N)*sqrt(curVar));
                    p2 = 1-normcdf(zscore);
                elseif strcmp(OrSt.gapTestWay, 'orderStats')
                    [mu, sigma] = ordStatApproxKsec(r1, edge_pixels);
                    sum_st = mean(r1)-mean(edge_pixels);
                    zscore = (sum_st - mu*sqrt(curVar))...
                        / (sigma*sqrt(curVar));
                    p1 = normcdf(zscore,'upper');
                    [mu, sigma] = ordStatApproxKsec(r2, edge_pixels);
                    sum_st = mean(r2)-mean(edge_pixels);
                    zscore = (sum_st - mu*sqrt(curVar))...
                        / (sigma*sqrt(curVar));
                    p2 = normcdf(zscore,'upper');
                elseif strcmp(OrSt.gapTestWay, 'localOrdStats')% problematic!
                    % only consider the x+2 pixels for correction
                    [mu, sigma] = ordStatApproxKsecWith0s(edge_pixels, [], invalid_nei);
                    curMu = mu*sqrt(curVar);
                    curSigma = sqrt(curVar*sigma^2 + curVar/length(r1));
                    sum_st = mean(r1)-mean(edge_pixels); 
                    zscore = (sum_st - curMu) / curSigma;
                    p1 = normcdf(zscore,'upper');
                    sum_st = mean(r2)-mean(edge_pixels);
                    % estimate the mean of neighbor and edge, mean of edge
                    % will always be negative
                    zscore = (sum_st - curMu) / curSigma;
                    p2 = normcdf(zscore,'upper');
                end
                %p2(isnan(p2)) = 0;
                %p1(isnan(p1)) = 0;
            end
            %figure;histogram(r2); hold on; histogram(r1);hold on; histogram(edge_pixels);legend;
            % zz = slabel;
            % zz(sComp>0) = sComp(sComp>0)+2; zzshow(label2rgb3d(zz))
            %fprintf('t-%d, reg: %d-vs-%d, p:%.3f, %.3f\n', i, pair_candidates(j,1),pair_candidates(j,2),p1,p2);
            if (~isnan(p1) && ~isnan(p2)) || (i==dilate_sc)
                p2(isnan(p2)) = 0;
                p1(isnan(p1)) = 0;
                if ~(p1 < p_thres/4 && p2 < p_thres/4) %p1 > p_thres% merge
                    % disp('This region should be merged');
                    mr_cnt = mr_cnt + 1;
                    merged_regions{mr_cnt} = pair_candidates(j,:);
                    %slabel(edge_locs) = pair_candidates(j,1);
                else
                    regComp(fg_locs(edge_all)) = 1; % all pixels along with two regions are edges
                    bk_cnt = bk_cnt + 1;
                    bk_regions{bk_cnt} = pair_candidates(j,:);
                    bk_perT(i) = bk_perT(i)+1;
%                     if i>min_r
%                         fprintf('!!!find a new region for broken!\n');
%                     end
                end
            end
        end
    end
end
merged_regions = merged_regions(1:mr_cnt);                                  % incase merge and break duplicated
bk_regions = bk_regions(1:bk_cnt);
for j=1:numel(merged_regions)
    for bb = 1:numel(bk_regions)
        if sum(abs(merged_regions{j}-bk_regions{bb}))==0
            merged_regions{j} = [];
            break;
        end
    end
end
for j=1:numel(merged_regions)                                               % merge merged region
    if isempty(merged_regions{j})
        continue;
    end
    for k = j+1:numel(merged_regions)
        if isempty(merged_regions{k})
            continue;
        end
        if ~isempty(intersect(merged_regions{j}, merged_regions{k}))
            merged_regions{j} = union(merged_regions{j}, merged_regions{k});
            merged_regions{k} = [];
        end
    end
end

newLabel = label_di;
for j=1:numel(merged_regions)
    if isempty(merged_regions{j})
        continue;
    end
    for k=2:length(merged_regions{j})
        newLabel(label_di==merged_regions{j}(k)) = merged_regions{j}(1);
    end
end

newLabel = newLabel(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc);
newLabel = rearrange_id(newLabel);
regComp = regComp(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc);
% remove boundary parts for further region-grow
newLabel(regComp>0) = 0;