function g = graphPara_cell(refine_res, frameNum)

cell_num_per_frame = cellfun(@(x) max(x(:)), refine_res);

% graph paramter
g.cycle_track = true; % true: circulation framework to solve tracking problem
g.stdCombined = 1;
g.maxDistXYZ = [25,25,5];
g.k = 3; % number of edges used before/after current particle pair
g.particleNum = sum(cell_num_per_frame);
g.transitionFactor = 1;% the weight of transition cost
g.validPre = 4; % check 2 previous point to get the mean
g.validPost = 4; % check 2 post point to get the mean
%g.maxEdgeNum = 4; at most check maxEdgeNum following or previous edges
g.timeJump = false; % Consier the jump x5 = x4+v*(t5-t4) , t5-t4 is not 1 
g.initEnter = 100; % initial 
g.realEnter = chi2inv(1-0.01/frameNum, 1)/2; % 12.3546 is from chi2inv(1-0.01/particleNum) / 2
g.c_en = g.realEnter;% cost of appearance and disappearance in the scene
g.c_ex = g.c_en;
g.observationCost = -(g.c_en+g.c_ex);
g.jumpCost = [1 1 1];%abs(g.observationCost/g.k); % how much we should punish jump frames
g.varEstMethod = 'independent'; % median and independent
g.costCalMethod='chi1Square'; % chi1Square:use 1df chi-squre, fisher: 2df chi-square,zscore: use z-score
g.trackLength4var = 5;%
g.truncatedGaussian = 0.2; % remove 20% of samples
g.varPrior = 100;% use gamma distribution as variance prior, use longest 100 tracks to estimate parameters
g.priorType = 'Gauss'; % gamma/weiAvg/Gauss/ scaled inverse chi-square  (sInvX2)
g.directionWise = false;% cal pvalue directionwise and use fisher method to combine 
g.dependencyFactor = 1.4;% the particles are not independent, based on their correlation, we add a correction factor for variance estimation.
g.maxIter = 5;% maximum number of iterations for min-cost flow tracking
g.minIntensity = 20;
g.side_flag = 0;% 1 means only consider pre-side (over-merge), 2: post side (over-segment), 0: both side
end