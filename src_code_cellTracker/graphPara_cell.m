function g = graphPara_cell(particleNum)
%generate the needed parameters for cell/region tracking
%the parameters include some basic information about detected cells
% contact ccwang@vt.edu

g.cycle_track = true; % true: circulation framework to solve tracking problem
g.stdCombined = 1; % we are cal difference between two locations (transitCostInitial_cell line 58)
g.maxDistXYZ = [50,50,5]; % max moving distance for a cell in adjacent frames
g.k = 3; % maximum number of jump allowed
g.particleNum = particleNum; % cell number
g.transitionFactor = 1;% the weight of transition cost
g.validPre = 4; % check 2 previous point to get the mean
g.validPost = 4; % check 2 post point to get the mean
%g.maxEdgeNum = 4; at most check maxEdgeNum following or previous edges
g.timeJump = false; % Consier the jump x5 = x4+v*(t5-t4) , t5-t4 may not be 1 
g.initEnter = 100; % initial enter/exit cost, force every node link to its neighbor
g.realEnter = chi2inv(1-0.01/g.particleNum, 1)/2; % 12.3546 is from chi2inv(1-0.01/particleNum) / 2
g.c_en = g.realEnter;% cost of appearance and disappearance in the scene
g.c_ex = g.c_en;
g.observationCost = -(g.c_en+g.c_ex)+0.00001; % make sure detections are all included
g.jumpCost = [];%abs(g.observationCost/g.k); % how much we should punish jump frames
g.varEstMethod = 'independent'; % median and independent
g.costCalMethod='chi1Square'; % chi1Square:use 1df chi-squre, fisher: 2df chi-square,zscore: use z-score
g.trackLength4var = 5; % tracks with smaller length will not be used to cal variance 
g.truncatedGaussian = 0.2; % remove 20% of extreme edges in variance estiamtion (truncated gaussian)
g.varPrior = 100;% use gamma distribution as variance prior, use longest 100 tracks to estimate gamma parameters
g.priorType = 'Gauss'; % prior distribution: gamma/weiAvg/Gauss/ scaled inverse chi-square  (sInvX2)
g.directionWise = false;% cal pvalue direction-wise and use fisher method (or chi-square) to combine 
g.dependencyFactor = 1.4;% the particles are not independent, based on their correlation, we add a correction factor for variance estimation.
g.splitMergeHandle = 'noJumpAll';% noJumpAll: all linking has no jump; noJump: no jump for spliting and merge; none: consider jump
g.maxIter = 7; % maximum iteration number
g.removeSingleRegion = true; % if a cell does not participate any trace, remove it.
g.detect_missing_head_tail = true; % add missing cell, do we consider head/tail node
g.applyDrift2allCoordinate = false; % correct drifting and change all y-, x- and z- coordinate
g.use_translation_as_drift = false; % use the translation as drfit.

g.blindmergeNewCell = false; % for a new detected region, blindly merge it to existing one if there is touch
g.simple_merge = true;% for two regions, if we decide to merge them, simply add their pixels together
g.par_kid_consistency_check = true; % when split a region in regionRefresh.m, check first if its two parents and two kids are consistent
g.reSplitTest = true; % for a new detected cell, if it is adjacent to another cell, re-split these two using their union voxels
g.stableNodeTest =true;
g.useOldArcProps = false;
g.considerBrokenCellOnly = true; % for linking allowing split/merge, does not consider nodes that has another good linkage already
g.addCellMissingPart = false; % if a cell missed part, we want to detect it, otherwise we can remove seeds that highly overlapped with an existing cell
g.splitMergeCost = true;% if cost of a+b->c is 20, then cost of a->c and b->c are set to 10 if true; otherwise both 20
end