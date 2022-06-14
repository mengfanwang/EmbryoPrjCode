function [movieInfo, movieInfoAll] = mcfTracking_cell(det_maps)
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
%% intial graph parameters
g = graphPara(length(movieInfo.xCoord));

movieInfo.Ci = zeros(g.particleNum,1)+g.observationCost; % use constant as observation cost
% initial transition cost
[neibIdx, Cij,edgeJumpCoeff] = transitCostInitial_cell(det_maps, det_YXZ,voxIdxList, movieInfo.frames,g);
movieInfo.nei = neibIdx;
movieInfo.Cij = Cij;
movieInfo.edgeJumpCoeff = edgeJumpCoeff;

%% iterative update transition cost start
% the initial result from min-cost flow
movieInfoAll = cell(g.maxIter,1);

loopCnt = 1;
while 1
    disp(loopCnt);
    if loopCnt <=2 % result of first(or &second) iteration are initialization
        g.c_en = g.initEnter;% cost of appearance and disappearance in 1st/second run are 100
    else
        g.c_en = g.realEnter;% cost of appearance and disappearance in the scene
    end
    g.c_ex = g.c_en;
    g.observationCost = -(g.c_en+g.c_ex);
    
    % build graph using current transition cost
    [~, g, dat_in] = trackGraphBuilder(movieInfo, g);
    
    % min-cost circulation for min-cost flow
    movieInfo = mccTracker(dat_in, movieInfo, g);
    movieInfoAll{loopCnt} = movieInfo;
    
    % update jump Ratio==> punishment to jump
    g.jumpCost = movieInfo.jumpRatio;

    if loopCnt>g.maxIter
        break;
    end
    % there can be drift in the data, we should correct it from tracking results
%     [movieInfo, dfEst] = driftFromTracks(movieInfo,g);
%     if dfEst == 0 % no enough tracks
%         break;
%     end
    % update transition cost based on current trajectories
    [movieInfo, zz] = transitCostGaussinBidirection_cell(movieInfo, g);

    loopCnt = loopCnt+1;
end