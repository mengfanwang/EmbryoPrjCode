function [movieInfo, movieInfoAll, dat_in_append] = mcfTracking_cell_v1(det_maps, g)
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
%% information initialization
[movieInfo,det_YXZ, voxIdxList] = tree2tracks_cell(det_maps, false);% false means we do not use existing tracking results

movieInfo.orgCoord = [movieInfo.xCoord, movieInfo.yCoord, movieInfo.zCoord];

movieInfo.Ci = zeros(g.particleNum,1)+g.observationCost; % use constant as observation cost
% initial transition cost
%[neibIdx, Cij, preNei, ovSize, preOvSize, Cji] = transitCostInitial_cell(det_maps, det_YXZ,voxIdxList, movieInfo,g);
movieInfo = transitCostInitial_cell(det_maps, det_YXZ,voxIdxList, movieInfo,g);

%% iterative update transition cost start
% the initial result from min-cost flow
movieInfoAll = cell(g.maxIter,1);
loopCnt = 1;
while 1
%     disp(loopCnt);
    if loopCnt <= 0 % result of first(or &second) iteration are initialization
        g.c_en = g.initEnter;% cost of appearance and disappearance in 1st/second run are 100
    else
        g.c_en = g.realEnter;% cost of appearance and disappearance in the scene
    end
    g.c_ex = g.c_en;
    g.observationCost = -(g.c_en+g.c_ex);
    
    % build graph using current transition cost
    [~, g, dat_in] = trackGraphBuilder_cell(movieInfo, g);
    movieInfo4gamma = mccTracker(dat_in, movieInfo, g);
    % update jump Ratio==> punishment to jump
    g.jumpCost = movieInfo4gamma.jumpRatio;
    
    movieInfo = transitCostUpt_cell(movieInfo4gamma, g);
    
    % adding extra edges to high-overlapped 
    [dat_in_append, movieInfo] = highOvEdges(movieInfo, g.side_flag);
    [~, g, dat_in] = trackGraphBuilder_cell(movieInfo, g);
    if ~isempty(dat_in_append)
        dat_in = cat(1, dat_in, dat_in_append);
    end
    % min-cost circulation for min-cost flow
    movieInfo = mccTracker(dat_in, movieInfo, g);
    
    movieInfoAll{loopCnt} = movieInfo;

    loopCnt = loopCnt+1;
    if loopCnt>g.maxIter
        break;
    end
    
end