function [newLabel, regComp] = edgeTest(vid, label, fgMap, connect, p_thres, OrSt)
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
vid_di = nan(h+6, w+6, z);
vid_di(4:end-3, 4:end-3, :) = vid;
label_di = zeros(h+6, w+6, z);
label_di(4:end-3, 4:end-3, :) = label;
fgMap_di = zeros(h+6, w+6, z);
fgMap_di(4:end-3, 4:end-3, :) = fgMap;

[h,w,z] = size(vid_di);
newLabel = zeros(h,w,z);
regComp = zeros(h,w,z);
if connect == 4
    neighbors = [-1, 1, -h, h];
elseif connect == 8
    neighbors = [-h-1, -h, -h+1, -1, 1, h-1, h, h+1];
elseif connect == 24
    neighbors = [-2*h-2:-2*h+2, -h-2:-h+2, -2,-1,1,2, h-2:h+2, 2*h-2:2*h+2];
else
    %connect == 49
    neighbors = [-3*h-3:-3*h+3, -2*h-3:-2*h+3, -h-3:-h+3, -3,-2,-1,1,2,3,...
        h-3:h+3, 2*h-3:2*h+3, 3*h-3:3*h+3];
end
gap = 3; % two circles of neighbors
for i=1:z
    im = vid_di(:,:,i);
    slabel = label_di(:,:,i);
    sFg = fgMap_di(:,:,i);
    sComp = zeros(h,w);
    real_edge_map = zeros(h,w);
    fg_locs = find(sFg>0);
    
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
            edge_locs = fg_locs(nei_segments(:,1) == pair_candidates(j,1) & ...
                nei_segments(:,2) == pair_candidates(j,2));
            if length(edge_locs) < 10
                p1 = 0;
                p2 = 0;
            else
                edge_pixels = im(edge_locs);

                edge_vals = slabel(edge_locs);
                slabel(edge_locs) = 0; % remove edges from labels
                % way 1 use all the pixels in the area
    %             r1 = im(slabel==pair_candidates(j,1));
    %             r2 = im(slabel==pair_candidates(j,2));
                % way 2 use the pixels adjacent to the boundaries
                neighbor_cells = neighbor_pixels(edge_locs, size(slabel), gap);
                neighbors_locs = cat(1, neighbor_cells{2:3});
                r1 = im(neighbors_locs(slabel(neighbors_locs)==pair_candidates(j,1)));
                r2 = im(neighbors_locs(slabel(neighbors_locs)==pair_candidates(j,2)));
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
                    [mu, sigma] = conventionalIntegralv1(edge_pixels, r1);
                    sum_st = mean(r1)-mean(edge_pixels);
                    zscore = (sum_st - mu*sqrt(OrSt.noiseVar))...
                        / (sigma*sqrt(OrSt.noiseVar));
                    p1 = 1-normcdf(zscore);
                    [mu, sigma] = conventionalIntegralv1(edge_pixels, r2);
                    sum_st = mean(r2)-mean(edge_pixels);
                    zscore = (sum_st - mu*sqrt(OrSt.noiseVar))...
                        / (sigma*sqrt(OrSt.noiseVar));
                    p2 = 1-normcdf(zscore);
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
            if ~(p1 < p_thres && p2 < p_thres) %p1 > p_thres% merge
                mr_cnt = mr_cnt + 1;
                merged_regions{mr_cnt} = pair_candidates(j,:);
                slabel(edge_locs) = pair_candidates(j,1);
            else
                real_edge_map(edge_locs) = 1;
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
    new_slabel = zeros(h,w);
    for j=1:length(ids)
        new_slabel(slabel==ids(j)) = newIds(j);
    end
%     nei_pair_segments = nan(size(fg_neighbors,1),2);
%     cnt1 = 0;
%     nei_single_segments = nan(size(fg_neighbors,1),2);
%     cnt2 = 0;
%     nei_segments = nan(size(fg_neighbors,1),2);
%     for j=1:size(fg_neighbors,1)
%         cur_neis = fg_neighbors(j,:);
%         cur_neis(cur_neis == 0) = [];
%         uni_neis = sort(unique(cur_neis));
%         if length(uni_neis)==2
%             nei_segments(j,:) = uni_neis;
%             cnt1 = cnt1 + 1;
%             nei_pair_segments(cnt1,:) = uni_neis;
%         elseif length(uni_neis)==1
%             cnt2 = cnt2 + 1;
%             nei_single_segments(cnt2,:) = [uni_neis uni_neis];
%         end
%     end
%     pair_candidates = cat(1, unique(nei_pair_segments(1:cnt1,:),'rows'), ...
%         unique(nei_single_segments(1:cnt2,:),'rows'));
%     %pair_candidates = unique(nei_segments,'rows');
%     pair_candidates(isnan(pair_candidates(:,1)),:) = [];
%     % test each pair to see if we can merge
%     new_slabel = zeros(h,w);
%     org_label = slabel;
%     regCnt = 0;
%     merged_regions = cell(size(pair_candidates,1),1);
%     for j=1:size(pair_candidates,1)
%         if pair_candidates(j,1) ~= pair_candidates(j,2)
%             edge_locs = fg_locs(nei_segments(:,1) == pair_candidates(j,1) & ...
%                 nei_segments(:,2) == pair_candidates(j,2));
%             edge_pixels = im(edge_locs);
% 
%             slabel(edge_locs) = 0; % remove edges from labels
%             r1 = im(slabel==pair_candidates(j,1));
%             r2 = im(slabel==pair_candidates(j,2));
% 
%             % now we are testing edge_pixels and pixels in area r1, r2
%             [~, p1] = ttest2(edge_pixels, r1);
%             [~, p2] = ttest2(edge_pixels, r2);
%             %[~, p0] = ttest2(R1, r2);
%             if ~(p1 < p_thres && p2 < p_thres) % merge
%                 regCnt = regCnt + 1;
%                 new_slabel(slabel==pair_candidates(j,1)) = regCnt;
%                 new_slabel(slabel==pair_candidates(j,2)) = regCnt;
%                 new_slabel(edge_locs) = regCnt;
%             else % not merge
%                 regCnt = regCnt + 1;
%                 new_slabel(org_label==pair_candidates(j,1)) = regCnt;
%                 regCnt = regCnt + 1;
%                 new_slabel(org_label==pair_candidates(j,2)) = regCnt;
%             end
%         else
%             tested_segments = pair_candidates(1:j-1,:);
%             if isempty(find(tested_segments(:)==pair_candidates(j,1),1))
%                 regCnt = regCnt + 1;
%                 new_slabel(org_label==pair_candidates(j,1)) = regCnt;
%             end
%         end
% 
%     end
    newLabel(:,:,i) = new_slabel;
    %regComp(:,:,i) = sComp;
    regComp(:,:,i) = real_edge_map;
end

newLabel = newLabel(4:end-3, 4:end-3,:);
regComp = regComp(4:end-3, 4:end-3,:);