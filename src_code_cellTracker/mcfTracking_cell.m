function [movieInfo, movieInfoAll, refine_res, refine_resAll,...
    threshold_res, threshold_resAll] = ...
    mcfTracking_cell(det_maps, embryo_vid, threshold_res,...
    varMaps, eigMaps, g, q)
% tracking based on min-cost flow/circulation
% INPUT:
% movieInfo: struct consists of
% 1. xCoord, yCoord, zCoord indicate detections locations
% 2. frames indicates detections time stamps
% det_maps: cells contains the label maps of all the frames
% OUTPUT:
% movieInfo: tracking results
% movieInfoAll: cells contains results from each iteration

% contact: ccwang@vt.edu, 02/26/2020

if isfield(q,'saveInterMediateRes')
    q.saveInterMediateRes = false;
end

if isfield(g,'use_translation_as_drift')
    if g.use_translation_as_drift
        g.translation = InitialDriftEstimate(embryo_vid);
    end
end

movieInfoAll = cell(g.maxIter * 3 + 2,1);
refine_resAll = cell(g.maxIter * 3 + 2,1);
threshold_resAll = cell(g.maxIter * 3 + 2,1);
%% information initialization
[movieInfo, det_YXZ, voxIdxList] = tree2tracks_cell(det_maps, q, false);% false means we do not use existing tracking results

movieInfo.orgCoord = [movieInfo.xCoord, movieInfo.yCoord, ...
    movieInfo.zCoord];

movieInfo.Ci = zeros(g.particleNum,1)+g.observationCost; % use constant as observation cost
% initial transition cost
movieInfo = transitCostInitial_cell(det_maps, det_YXZ, voxIdxList, ...
    movieInfo,g);
% do a simple linking to estimate the drift of the data
[~, g, dat_in] = trackGraphBuilder_cell(movieInfo, g);
movieInfo4gamma = mccTracker(dat_in, movieInfo, g, false, false);
g.jumpCost = movieInfo4gamma.jumpRatio; % no need here
movieInfo = transitCostUpt_cell(movieInfo4gamma, g); % drift corrected

refine_res = det_maps;

if q.removeSamllRegion
    % remove tiny cells
    [movieInfo, refine_res, threshold_res, invalid_cell_ids] = ...
        removeTinyCells(movieInfo, refine_res, threshold_res, q);
    % 2.3 update movieInfo and refine_res
    [movieInfo, refine_res, g] = movieInfo_update(movieInfo, ...
        refine_res, invalid_cell_ids, g);
end

if q.saveInterMediateRes
    movieInfoAll{1} = movieInfo;
    refine_resAll{1} = refine_res;
    threshold_resAll{1} = threshold_res;
end
res_cnt = 1;
% build a map with label the real id
% incre_cnt = 0;
% for i=1:numel(refine_res)
%     idx = find(refine_res{i}>0);
%     refine_res{i}(idx) = refine_res{i}(idx) + incre_cnt;
%     incre_cnt = max(refine_res{i}(idx));
% end



%% debug
% [~, gt_tracks_ov] = loc2id(g.gt_mat, refine_res);
% [~, ratio_mat, track_id_head, ~] = validate_res(movieInfo, gt_tracks_ov);
% bestOvRatio = nanmax(ratio_mat,[],2);
% if ~isempty(find(isnan(bestOvRatio),1))
%     keyboard;
% end
% it = 2;
% movieInfo = movieInfoAll{it};
% refine_res = refine_resAll{it};
% threshold_res = threshold_resAll{it};

% add an invalid gap map
movieInfo.validGapMaps = cell(numel(refine_res),1);
for i=1:numel(refine_res)
    movieInfo.validGapMaps{i} = true(size(refine_res{i}));
end
%% iterative update transition cost start
% the initial result from min-cost flow
loopCnt = 1;
% [movieInfo, refine_res, threshold_res, g] = ...
%         missing_cell_module(movieInfo, refine_res, threshold_res, ...
%         embryo_vid, eigMaps, varMaps, g, q);
while 1
    % debug;
    %     ii = 2;
    %     id = refine_resAll{ii}{3}(1326231) + sum(movieInfoAll{ii}.n_perframe(1:2));
    %     disp([refine_resAll{ii}{3}(1326231), ...
    %         unique(refine_resAll{ii}{3}(movieInfoAll{ii}.voxIdx{id}))]);
    %     id = refine_res{3}(1326231) + sum(movieInfo.n_perframe(1:2));
    %     disp([refine_res{3}(1326231), ...
    %         unique(refine_res{3}(movieInfo.voxIdx{id}))]);
    tic;
    for dummy = 1:3
        %% save intermediate results
        res_cnt = res_cnt + 1;
        if q.saveInterMediateRes
            refine_resAll{res_cnt} = refine_res;
            threshold_resAll{res_cnt} = threshold_res;
        end
        %% module one: split/merge (optional: remove small regions)
        [movieInfo, refine_res, threshold_res, g, movieInfo_noJump] = ...
            split_merge_module(movieInfo, refine_res, threshold_res, ...
            embryo_vid, eigMaps, g, q);

        %% save intermediate results
        if q.saveInterMediateRes
            movieInfoAll{res_cnt} = movieInfo_noJump;
        end

    end
    toc;
    %% module two: missing cells (optional: remove short tracks)
    [movieInfo, refine_res, threshold_res, g] = ...
        missing_cell_module(movieInfo, refine_res, threshold_res, ...
        embryo_vid, eigMaps, varMaps, g, q);
    loopCnt = loopCnt + 1;
    if loopCnt > g.maxIter
        break;
    end
end

% finally use one to one linking
[~, g, dat_in] = trackGraphBuilder_cell(movieInfo, g);
movieInfo = mccTracker(dat_in, movieInfo, g);
% !!! LAST: test if broken tracks can be merged
movieInfo = mergeBrokenTracksV2(movieInfo, g);
res_cnt = res_cnt + 1;
if q.saveInterMediateRes
    refine_resAll{res_cnt} = refine_res;
    threshold_resAll{res_cnt} = threshold_res;
    movieInfoAll{res_cnt} = movieInfo;
end
% measure accuracy as side evidence
if isfield(g, 'gt_mat')
    if ~isempty(g.gt_mat)
        accuracy_measure(movieInfo, refine_res, g.gt_mat);
    end
end
