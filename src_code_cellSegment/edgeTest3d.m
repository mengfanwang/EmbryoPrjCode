function [newLabel, regComp] = edgeTest3d(vid, label, fgMap, connect, p_thres, OrSt)
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
dilate_sc = 5;
vid_di = nan(h+dilate_sc*2, w+dilate_sc*2, z+dilate_sc*2);
vid_di(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc) = vid;
label_di = zeros(size(vid_di));
label_di(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc) = label;
fgMap_di = zeros(size(vid_di));
fgMap_di(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc) = fgMap;

[h,w,z] = size(vid_di);
newLabel = zeros(h,w,z);
regComp = zeros(h,w,z);
if connect == 26
    neighbors = [-h-1, -h, -h+1, -1, 1, h-1, h, h+1];
    neighbors = [neighbors-h*w, -h*w, neighbors, h*w, neighbors+h*w];
elseif connect == 124
    neighbors = [-2*h-2:-2*h+2, -h-2:-h+2, -2,-1,1,2, h-2:h+2, 2*h-2:2*h+2];
    neighbors = [neighbors-2*h*w, neighbors-h*w, -2*h*w, -h*w, neighbors...
        , 2*h*w, h*w, neighbors+h*w, neighbors+2*h*w];
elseif connect == 342
    neighbors = [-3*h-3:-3*h+3, -2*h-3:-2*h+3, -h-3:-h+3, -3,-2,-1,1,2,3,...
        h-3:h+3, 2*h-3:2*h+3, 3*h-3:3*h+3];
    neighbors = [neighbors-3*h*w, neighbors-2*h*w, neighbors-h*w,...
        -3*h*w, -2*h*w, -h*w, neighbors, 3*h*w, 2*h*w, h*w, ...
        neighbors+h*w, neighbors+2*h*w, neighbors+3*h*w];
elseif connect == 146
    neighbors = [-3*h-3:-3*h+3, -2*h-3:-2*h+3, -h-3:-h+3, -3,-2,-1,1,2,3,...
        h-3:h+3, 2*h-3:2*h+3, 3*h-3:3*h+3];
    neighbors = [neighbors-h*w, -h*w, neighbors, h*w, neighbors+h*w];
elseif connect == 48 % only 2d neighbors are considered
    neighbors = [-3*h-3:-3*h+3, -2*h-3:-2*h+3, -h-3:-h+3, -3,-2,-1,1,2,3,...
        h-3:h+3, 2*h-3:2*h+3, 3*h-3:3*h+3];
elseif connect == 66 % only 2d neighbors are considered
    nei8 = [-h-1, -h, -h+1, -1, 1, h-1, h, h+1];
    neighbors = [-3*h-3:-3*h+3, -2*h-3:-2*h+3, -h-3:-h+3, -3,-2,-1,1,2,3,...
        h-3:h+3, 2*h-3:2*h+3, 3*h-3:3*h+3, nei8-h*w, -h*w, h*w, nei8+h*w];
elseif connect == 80 % only 2d neighbors are considered
    neighbors = [-4*h-4:-4*h+4,-3*h-4:-3*h+4, -2*h-4:-2*h+4, -h-4:-h+4, -4, -3,-2,-1,1,2,3,4,...
        h-4:h+4, 2*h-4:2*h+4, 3*h-4:3*h+4, 4*h-4:4*h+4];
else
    error('gap connection not defined');
end
gap = 5; % two circles of neighbors

slabel = label_di;
sComp = zeros(h,w,z);

fg_locs = find(fgMap_di>0);

nei_mat = repmat(fg_locs, 1, length(neighbors)) +...
    repmat(neighbors, length(fg_locs), 1);
fg_neighbors = slabel(nei_mat);
%     fg_neighbors(fg_neighbors==0) = nan;
%     find(nanmin(fg_neighbors, 2) < nanmean(fg_neighbors, 2));
%     neiCnt = sum(diff(sort(A,2),1,2)~=0,2)+1;
%     fg_neighbors = fg_neighbors();
nei_segments = nan(size(fg_neighbors,1),2);
for j=1:size(fg_neighbors,1)
    cur_neis = fg_neighbors(j,:);
    cur_neis(cur_neis == 0) = [];
    uni_neis = sort(unique(cur_neis));
    if length(uni_neis)==2
        nei_segments(j,:) = uni_neis;
    end
end
pair_candidates = unique(nei_segments,'rows');
pair_candidates(isnan(pair_candidates(:,1)),:) = [];
% test each pair to see if we can merge
merged_regions = cell(size(pair_candidates,1),1);
mr_cnt = 0;

for j=1:size(pair_candidates,1)
    if pair_candidates(j,1) ~= pair_candidates(j,2)
        edge_all = nei_segments(:,1) == pair_candidates(j,1) & ...
            nei_segments(:,2) == pair_candidates(j,2);
        edge_locs = fg_locs(edge_all);% & label_di(fg_locs) == 0
        if length(edge_locs) < 10 % means two regions are too far away:real gap
            p1 = 0;
            p2 = 0;
        else
            edge_pixels = vid_di(edge_locs);
            
            edge_vals = slabel(edge_locs);
            slabel(edge_locs) = 0; % remove edges from labels
            % use the pixels adjacent to the boundaries
            neighbor_cells = neighbor_pixels(edge_locs, size(slabel), gap);
            neighbors_locs = cat(1, neighbor_cells{3:gap});
            r1 = vid_di(slabel==pair_candidates(j,1));%neighbors_locs(slabel(neighbors_locs)==pair_candidates(j,1)));
            r2 = vid_di(slabel==pair_candidates(j,2));%neighbors_locs(slabel(neighbors_locs)==pair_candidates(j,2)));
            
            slabel(edge_locs) = edge_vals;
            sComp(edge_locs) = 1;
            sComp(neighbors_locs(slabel(neighbors_locs)==pair_candidates(j,1))) = 2;
            sComp(neighbors_locs(slabel(neighbors_locs)==pair_candidates(j,2))) = 3;
            % now we are testing edge_pixels and pixels in area r1, r2
            %[~, p1] = ttest2(edge_pixels, [r1;r2]);
            var_flag = false;
            if isempty(OrSt.noiseVar)
                var_flag = true;
                OrSt.noiseVar = var(cat(1,r1,r2,edge_pixels));
            end
            p1all = []; p2all = [];
            if ~isfield(OrSt,'mu') || ~isfield(OrSt,'sigma')
                [~, p1] = ttest2(edge_pixels, r1);
                [~, p2] = ttest2(edge_pixels, r2);
                p1all = p1;
                p2all = p2;
            elseif ~isempty(OrSt.mu) && ~isempty(OrSt.sigma)
                M = max(min(length(r1), size(OrSt.mu,1)), 10);
                N = max(min(length(edge_pixels), size(OrSt.mu,1)), 10);
                sum_st = mean(r1)-mean(edge_pixels);
                zscore = (sum_st - OrSt.mu(M,N)*sqrt(OrSt.noiseVar))...
                    / (OrSt.sigma(M,N)*sqrt(OrSt.noiseVar));
                p1 = 1-normcdf(zscore);
                
                M = max(min(length(r2), size(OrSt.mu,1)), 10);
                sum_st = mean(r2)-mean(edge_pixels);
                zscore = (sum_st - OrSt.mu(M,N)*sqrt(OrSt.noiseVar))...
                    / (OrSt.sigma(M,N)*sqrt(OrSt.noiseVar));
                p2 = 1-normcdf(zscore);
                
                p1all = [p1all p1];
                p2all = [p2all p2];
            else
                [mu, sigma] = ordStatApproxKsec(r1, edge_pixels);
                sum_st = mean(r1)-mean(edge_pixels);
                zscore = (sum_st - mu*sqrt(OrSt.noiseVar))...
                    / (sigma*sqrt(OrSt.noiseVar));
                p1 = normcdf(zscore,'upper');
                [mu, sigma] = ordStatApproxKsec(r2, edge_pixels);
                sum_st = mean(r2)-mean(edge_pixels);
                zscore = (sum_st - mu*sqrt(OrSt.noiseVar))...
                    / (sigma*sqrt(OrSt.noiseVar));
                p2 = normcdf(zscore,'upper');
                p1all = [p1all p1];
                p2all = [p2all p2];
            end
            p2(isnan(p2)) = 0;
            p1(isnan(p1)) = 0;
            if var_flag
                OrSt.noiseVar = [];
            end
            %disp([p1all;p2all]);
        end
        %[~, p0] = ttest2(R1, r2);
        %figure;histogram(r2); hold on; histogram(r1);hold on; histogram(edge_pixels);legend;
        %zz = slabel;zz(edge_locs) = zz(edge_locs)+3;zzshow(label2rgb3d(zz))
        % disp([p1, p2]);
        if ~(p1 < p_thres && p2 < p_thres) %p1 > p_thres% merge
            % disp('This region should be merged');
            mr_cnt = mr_cnt + 1;
            merged_regions{mr_cnt} = pair_candidates(j,:);
            slabel(edge_locs) = pair_candidates(j,1);
        else
            regComp(fg_locs(edge_all)) = 1; % all pixels along with two regions are edges
        end
    end
end
merged_regions = merged_regions(1:mr_cnt);
for j=1:numel(merged_regions)
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
ids = unique(slabel(slabel>0));
newIds = ids;
for j = 1:numel(merged_regions)
    if isempty(merged_regions{j})
        continue;
    end
    for k=2:length(merged_regions{j})
        newIds(ids==merged_regions{j}(k)) = merged_regions{j}(1);
    end
end
newLabel = zeros(h,w,z);
for j=1:length(ids)
    newLabel(slabel==ids(j)) = newIds(j);
end

newLabel = newLabel(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc);
regComp = regComp(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc);
% remove boundary parts for further region-grow
newLabel(regComp>0) = 0; 