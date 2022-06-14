function [fg, threshold] = ordstat4fg(vidComp, vid_sm, ...
    seedRegion, other_id_map, fmapComp, OrSt, minSize)
% given a seed region, find the most significant region include this region
% and view it as the foreground for further segmentation
if nargin < 7
    minSize = 100;
end
if strcmp(OrSt.imProcMethod, 'stb')
    corrFactor = OrSt.corrTerm; % scale coefficient
else
    corrFactor = 1;
end
suplb = 30*corrFactor;
infub = 5*corrFactor;
step_i = -1*corrFactor;
    
ub = max(suplb, round(mean(vidComp(seedRegion>0)))); % upper bound intensity level
seedSz = length(find(seedRegion>0));
lb = max(min(vidComp(fmapComp>0)), infub);
if isempty(lb) || isempty(ub) || ub+step_i <= lb
    fg = fmapComp;
    return;
end
%validNeiMap = ~otherIdMap;
regCandidates = cell(round(abs((ub-lb)/step_i))+1, 1);
%OrSt.fgTestWay = 'ttest';
test_way = OrSt.fgTestWay;%'Integral';%
if isempty(test_way)
    z = nan(numel(regCandidates), 6);
%     z_debias = nan(numel(regCandidates), 6);
%     sz_nei = nan(numel(regCandidates), 6);
else
    z = nan(numel(regCandidates), 1);
%     z_debias = nan(numel(regCandidates), 1);
%     sz_nei = nan(numel(regCandidates), 1);
end
%vid_sm = imgaussfilt3(vidComp, [1 1 0.5]);
% absdif = [];

% since we dilate 3 times for boundary and neighbors selection
sph = strel('sphere', 2);
se = strel(sph.Neighborhood(:,:,3));% equal to strel('disk',2)
otherIdTerritory = imdilate(other_id_map,se);
otherIdTerritory = imdilate(otherIdTerritory,strel('sphere',1));
otherIdTerritory(seedRegion>0) = 0; % seedregion has higher priority

test_stat = [];sz_nei = []; s = [];
thrs = ub:step_i:lb;
for cnt = 1:length(thrs) % parfor is even slower, 2 times slower
    thr = thrs(cnt);
    fgIn = vid_sm >= thr & fmapComp;
    %% way 1 select the region containing seed region
    if strcmp(OrSt.refine_cc, 'single_seed')
        fgIn(otherIdTerritory) = 0; % remove territory of other seeds      
        l = bwlabeln(fgIn, 6);
        seed_ids = l(seedRegion>0);
        seed_ids = unique(seed_ids(seed_ids>0));
        if isempty(seed_ids)
            continue;
        end
        %all_ids = unique(l(l>0));
        regCandidates{cnt} = ismember(l, seed_ids);
        curValidNeiMap = ~otherIdTerritory & ~fgIn & fmapComp;
    %     extraSeedMap = regCandidates{cnt} & otherIdMap;
    %     if ~isempty(find(extraSeedMap, 1))
    %         labelMap = zeros(size(seedRegion));
    %         labelMap(seedRegion>0) = 1;
    %         labelMap(extraSeedMap) = 2;
    %         newLabel = regionGrow(labelMap, scoreComp, fgIn, connect, cost_design,bg2sinkLink);
    %         regCandidates{cnt} = newLabel == 1;
    %     end
    else
        %% way 2: directly use all regions
        regCandidates{cnt} = fgIn;
        curValidNeiMap = fmapComp;
    end
    %% test significance
    coverRatio = length(find(regCandidates{cnt} & seedRegion>0))/seedSz;
    if coverRatio<0.5 %less than 50% of seed region is covered
        continue;
    end
%     [z(cnt,:), z_debias(cnt,:), sz_nei(cnt,:)]= regionSig(regCandidates{cnt}, vidComp,...
%         fmapComp, curValidNeiMap, OrSt);
    [z(cnt), test_stat(cnt), sz_nei(cnt), s(cnt)]= regionSig(regCandidates{cnt}, vidComp,...
        fmapComp, curValidNeiMap, OrSt);
%     if isnan(z(cnt,1))
%         break;
%     end
end
[v_od, od] = nanmax(z, [], 1);
% if strcmp(wayChoice, 'ttest_varKnown')
%     figure;subplot(2,2,1);plot(z); title('zscore');
%     subplot(2,2,2);plot(test_stat); title('mean-diff');
%     subplot(2,2,3);plot(sz_nei); title('size-nei');
%     subplot(2,2,4);plot(s); title('sigma');
% else
%     figure;subplot(2,2,1);plot(z); title('zscore');
%     subplot(2,2,2);plot(test_stat); title('bias');
%     subplot(2,2,3);plot(sz_nei); title('size-nei');
%     subplot(2,2,4);plot(s); title('sigma');
% end

% [~,max_od] = max(z_debias);
% title(num2str(max_od));
% subplot(1,3,3);plot(sz_nei);[~,max_od] = max(sz_nei);title(num2str(max_od));
% sz_nei(isnan(sz_nei)) = [];
% title(num2str(sz_nei(1)/sz_nei(4)));
%od = pickUpBestOder(z, z_debias, sz_nei);
%disp(od);
% if strcmp(OrSt.imProcMethod,'stb')
%     figure;plot([ub:step_i:lb].^2-3/8, z_debias);
% else
%     figure;plot([ub:step_i:lb], z_debias);
% end

if size(z, 2) == 6 % compare 5 testing methods
    fprintf('pick %d, %d, %d, %d, %d, %d, out of %d thresholds\n', ...
        od(1), od(2), od(3), od(4), od(5), od(6), size(z,1));
    od = od(4);
end

if isnan(v_od)
    fg = false(size(vidComp));
    threshold = nan;
    return;
end
threshold = thrs(od);
fg_reg = regCandidates{od};
fg = refine_with_seedRegion(fg_reg, seedRegion, minSize);
if isempty(find(fg & seedRegion,1))% the seedregion no longer related to fg
    fg = false(size(vidComp));
    threshold = nan;
    return;
end