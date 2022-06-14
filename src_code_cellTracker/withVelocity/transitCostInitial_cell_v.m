function [neibIdx, Cij, edgeJumpCoeff] = transitCostInitial_cell(det_maps, dets_YXZ, voxIdxList, frames, g)
%we can learn initial variance and mean of distances between observed object 
%locations and estimated locations from data
%INPUT:
% det_maps: cells containing label maps of all frames
% dets_YXZ: cells containing center locations of regions for each frame
% voxIdxList: cells containing voxel indexes of regions for each frame
% frames: the time/frame infor of all detections 
% g: parameters needed
%OUTPUT:
% neibIdx: neighbors that linked to each detection
% Cij: the overlapping ratio based distances
% edgeJumpCoeff: ratio of jumps among all the linkings

% contact: ccwang@vt.edu, 02/26/2020
xCoord = cell(numel(dets_YXZ),1);
yCoord = cell(numel(dets_YXZ),1);
zCoord = cell(numel(dets_YXZ),1);

for i=1:numel(dets_YXZ)
    yCoord{i} = dets_YXZ{i}(:,1);
    xCoord{i} = dets_YXZ{i}(:,2);
    zCoord{i} = dets_YXZ{i}(:,3);
end
%% using all the data to initialize transition cost
neibIdx = cell(g.particleNum, 1);
CDist = cell(g.particleNum, 1);
CX = cell(g.particleNum, 1);
CY = cell(g.particleNum, 1);
CZ = cell(g.particleNum, 1);
timepts = numel(dets_YXZ);
curLen = 0;
for i=1:timepts
    disp(i);
    for j=i+1:i+g.k
        if j>timepts
            break;
        end
%         DX = abs(createDistanceMatrix(xCoord{i},xCoord{j}));
%         DY = abs(createDistanceMatrix(yCoord{i},yCoord{j}));
%         DZ = abs(createDistanceMatrix(zCoord{i},zCoord{j}));
        DX = distanceMat2Mex(xCoord{i},xCoord{j});
        DY = distanceMat2Mex(yCoord{i},yCoord{j});
        DZ = distanceMat2Mex(zCoord{i},zCoord{j});
        if isfield(g, 'translation')
            [neighbors, distances] = ovDistanceMap(det_maps{i}, det_maps{j}, ...
                voxIdxList{i}, voxIdxList{j},...
                g.translation(j,:)-g.translation(i,:));
        else
            [neighbors, distances] = ovDistanceMap(det_maps{i}, det_maps{j}, ...
                voxIdxList{i}, voxIdxList{j});
        end
        bsIdxLen = sum(cellfun(@numel, voxIdxList(1:j-1)));
        for nn = 1:numel(neighbors)
            tmpIdx = neighbors{nn};
            if ~isempty(tmpIdx)
                neibIdx{curLen+nn} = cat(1,neibIdx{curLen+nn}, bsIdxLen+tmpIdx);
                CDist{curLen+nn} = cat(1,CDist{curLen+nn}, distances{nn});
                CX{curLen+nn} = cat(2,CX{curLen+nn}, DX(nn,tmpIdx));
                CY{curLen+nn} = cat(2,CY{curLen+nn}, DY(nn,tmpIdx));
                CZ{curLen+nn} = cat(2,CZ{curLen+nn}, DZ(nn,tmpIdx));
            end
        end
    end
    curLen = curLen+numel(voxIdxList{i});
end
validPts = ~cellfun(@isempty,CDist);
nearestDist = cellfun(@min,CDist(validPts));
nearestDX = cellfun(@min,CX(validPts));% all positive, should add sign?
nearestDY = cellfun(@min,CY(validPts));
nearestDZ = cellfun(@min,CZ(validPts));
if strcmp(g.varEstMethod,'median')
    stdFromAllPts = [nanmedian(abs(nearestDX-nanmean(nearestDX))), ...
        nanmedian(abs(nearestDY-nanmean(nearestDY))), ...
        nanmedian(abs(nearestDZ-nanmean(nearestDZ)))];
else
    stdFromAllPts = [std(nearestDX) std(nearestDY) std(nearestDZ)];
end
% fit the gamma distribution using current overlapping distances
phat = gamfit(nearestDist);

if g.stdCombined % because we are cal the diff among two locations
    stdFromAllPts = stdFromAllPts/sqrt(2);
end
Cij = cell(g.particleNum, 1);
edgeJumpCoeff = cell(g.particleNum, 1);
for i=1:numel(CDist)
    if isempty(CDist{i})
        continue;
    end
%     % start calculate cost Cij
    neipos = [CX{i}',CY{i}',CZ{i}'];
    curpos = neipos*0;
    % center_shift_cost
    [csc,~,p_csc] = distance2EdgeCost(neipos,curpos, stdFromAllPts, stdFromAllPts, ones(3,2),g);
    %csc = -norminv(p_csc);
    % overlapping_cost
    p_oc = 1-gamcdf(CDist{i}, phat(1), phat(2));
    oc = norminv(1-p_oc/2);
    p_both = 1-chi2cdf(oc.^2 + csc,2);
    
    Cij{i} = chi2inv(1-p_both,1);%
    if g.timeJump
        edgeJumpCoeff{i} = zeros(size(neipos,1), 3);
        for j=1:size(neipos,1)
            edgeJumpCoeff{i}(j,:) = frames(neibIdx{i}(j))-frames(i);
        end
    else
        edgeJumpCoeff{i} = ones(size(neipos,1), 3);
    end
end
fprintf('finish trajectory intialization with purely distance!\n');

end