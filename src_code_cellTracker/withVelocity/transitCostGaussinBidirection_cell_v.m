function [movieInfo, zz] = transitCostGaussinBidirection_cell(movieInfo, g)
% update transition cost using existing trajectories using both precursors
% and ancestors
%
xCoord = movieInfo.xCoord;
yCoord = movieInfo.yCoord;
zCoord = movieInfo.zCoord;
numTrajectories = numel(movieInfo.tracks);
%Dij = movieInfo.Dij;% the distances of neighbors
zz = zeros(1e7,2);% test variable
zzCnt = 1;
%% estimate variance: we have two ways to estimate the variance of all trajectory motion
if ~isfield(g,'varEstMethod')
    fprintf('The variance estimation method is not specified!\n');
    fprintf('We use default one: median.\n');
    g.varEstMethod = 'median';
end
if strcmp(g.varEstMethod,'median')
    [devXYZ, ~,validTrack] = getStdFromTracks(movieInfo,g);
    phatDist = getStdFromTracks_cell(movieInfo, g);
    % debug
    trackStd = zeros(size(devXYZ,1),3);
    for i=1:size(devXYZ,1)
        if size(devXYZ{i},1)>1
            trackStd(i,:) = std(devXYZ{i});
        end
    end
    trackStd2 = trackStd(trackStd(:,1)~=0,:);
    % end debug
    devAll = cat(1, devXYZ{:});
    devAll = bsxfun(@minus, devAll, mean(devAll,1));
    devAll = abs(devAll);
    stdXYZ = nanmedian(devAll);
    disp(stdXYZ);
    disp(validTrack);
    stdXYZ = repmat(stdXYZ, numTrajectories, 1);
elseif strcmp(g.varEstMethod,'independent')
    % way 2: esimate variance independently for each trajectory
    [~, stdXYZ,~] = getStdFromTracks(movieInfo, g);
    phatDist = getStdFromTracks_cell(movieInfo, g);
    medStd = nanmedian(stdXYZ);
    % shorter ones use median
    shortTrackId = isnan(stdXYZ(:,1));
    stdXYZ(shortTrackId,1) = medStd(1);
    stdXYZ(shortTrackId,2) = medStd(2);
    stdXYZ(shortTrackId,3) = medStd(3);
else
    error('Unknown variance estimation method!\n');
end
%% assign std to each particle
movieInfo.particleStd = zeros(length(movieInfo.frames),3);
movieInfo.particleStd(:,1) = medStd(1);
movieInfo.particleStd(:,2) = medStd(2);
movieInfo.particleStd(:,3) = medStd(3);
for i=1:length(movieInfo.xCoord)
    curTrackNum = movieInfo.particle2track(i,1);
    if ~isnan(curTrackNum)
        movieInfo.particleStd(i,:) = stdXYZ(curTrackNum,:);
    end
end
%% estimate inherent variance of v
if g.timeJump
    [~,v0Var, s] = v2v1VarCal(movieInfo, g);
    %validOne = v2v1Var(:,1)==g.validPre & v2v1Var(:,2)==g.validPre & v2v1Var(:,3)==1;
    %movieInfo.v0Var = v0Var;%v2v1Var(validOne,8:10);
    movieInfo.s = s;
    fprintf('current s=[%.2f, %.2f, %.2f]!\n',movieInfo.s(41,3)...
        ,movieInfo.s(41,5),movieInfo.s(41,7));
end
%% we update cost of all edges rather than part of them
particle2track = movieInfo.particle2track;
nei = movieInfo.nei;
newCij = cell(length(movieInfo.xCoord),1);
edgeJumpCoeff = cell(length(movieInfo.xCoord),1);
edgeVelocity = cell(length(movieInfo.xCoord),1);
E2VVarV = cell(length(movieInfo.xCoord),1);
for i=1:length(movieInfo.xCoord)
    curNode = i;
    curTrackNum = particle2track(i,1);
    if isnan(curTrackNum)
        curStd = medStd;
    else
        curStd = stdXYZ(curTrackNum,:);
    end
    neiUp = nei{i};
    newCijSingle = nan(length(neiUp),1);
    edgeJumpCoeff{i} = nan(length(neiUp),3);
    E2VVarV{i} = nan(length(neiUp),7);
    for nn=1:length(neiUp)
        % estimate only one position using both trajectories
        neiPosition = [xCoord(neiUp(nn)),yCoord(neiUp(nn)),zCoord(neiUp(nn))];
        neiTrackNum = particle2track(neiUp(nn),1);
        if isnan(neiTrackNum)
            neiStd = medStd;
        else
            neiStd = stdXYZ(neiTrackNum,:);
        end
        timeGap = movieInfo.frames(neiUp(nn))-movieInfo.frames(curNode);
        [~, uCur, curRatioSquare, velocity,jumpCoeff, tmpE2VVarV] = ...
            estimateMeanFrom2Tracks(movieInfo, curNode,neiUp(nn), g);
        if sum(isnan(jumpCoeff))>0
            edgeJumpCoeff{i}(nn,:) = timeGap;
            edgeVelocity{i}(nn,:) = [0 0 0];
        else
            edgeJumpCoeff{i}(nn,:) = jumpCoeff;
            edgeVelocity{i}(nn,:) = velocity;
            E2VVarV{i}(nn,:) = tmpE2VVarV;
        end
        % center_shift_cost
        [csc,~,p_csc] = distance2EdgeCost(neiPosition, uCur,  curStd, neiStd, curRatioSquare,g,timeGap);
        %csc = -norminv(p_csc);    
        % overlapping_cost
        curCoord = [movieInfo.xCoord(curNode) movieInfo.yCoord(curNode) movieInfo.zCoord(curNode)];
        devXYZ = uCur-curCoord;
        ov_dist = ovDistanceRegion(movieInfo.vox{curNode}+round(devXYZ),...
            movieInfo.vox{neiUp(nn)});

        p_oc = 1-gamcdf(ov_dist, phatDist(1), phatDist(2));
        oc = norminv(1-p_oc/2);
        p_both = 1-chi2cdf(oc.^2 + csc,2);
        
        edgeCost = chi2inv(1-p_both,1);%
        
        newCijSingle(nn) = edgeCost;
        if isnan(edgeCost)
            fprintf('We found NaN edge cost!!!\n');
        end
    end
    newCij{i} = newCijSingle;
end
movieInfo.Cij = newCij;
movieInfo.edgeJumpCoeff = edgeJumpCoeff;
movieInfo.edgeVelocity = edgeVelocity;
% show costs of those edges related to trajectories
for i=1:numTrajectories
    curTrack = movieInfo.tracks{i};
    curStd = stdXYZ(i,:);
    if ~isempty(find(isnan(curStd) | curStd==0, 1))
        continue;
    end
    % start update cost Cij
    for j=1:length(curTrack)-1
        neiUp = nei{curTrack(j)};
        tmpNeiPos = neiUp==curTrack(j+1);
        zz(zzCnt,:) = movieInfo.Cij{curTrack(j)}(tmpNeiPos);
        zzCnt = zzCnt+1;
    end
end

zz = zz(1:zzCnt-1,:);

fprintf('finish cost update!\n');

end