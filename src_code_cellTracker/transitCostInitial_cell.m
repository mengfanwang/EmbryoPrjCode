function movieInfo = transitCostInitial_cell(det_maps, dets_YXZ, ...
    voxIdxList, movieInfo, g)
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

%% using all the data to initialize transition cost
% if ~exist('intial_distance.mat','file')
fprintf('start build intial distance and neighbor cells\n');
neibIdx = cell(g.particleNum, 1);
CDist = cell(g.particleNum, 1);
CDist_dirwise = cell(g.particleNum, 1);

timepts = numel(dets_YXZ);
curLen = 0;
for i=1:timepts
    %disp(i);
    for j=i+1:i+g.k
        if j>timepts
            break;
        end
        if isfield(g, 'translation')
            [neighbors, distances, distances_dirWise] = ovDistanceMap(...
                det_maps{i}, det_maps{j}, voxIdxList{i}, voxIdxList{j},...
                g.translation(j,:)-g.translation(i,:));
        else
            [neighbors, distances, distances_dirWise] = ovDistanceMap(...
                det_maps{i}, det_maps{j}, voxIdxList{i}, voxIdxList{j});
        end
        bsIdxLen = sum(movieInfo.n_perframe(1:j-1));
        for nn = 1:numel(neighbors)
            tmpIdx = neighbors{nn};
            if ~isempty(tmpIdx)
                neibIdx{curLen+nn} = cat(1,neibIdx{curLen+nn}, bsIdxLen+tmpIdx);
                CDist{curLen+nn} = cat(1,CDist{curLen+nn}, distances{nn});
                CDist_dirwise{curLen+nn} = cat(1,CDist_dirwise{curLen+nn}, ...
                    distances_dirWise{nn});
            end
        end
    end
    curLen = curLen+numel(voxIdxList{i});
end
%     save('intial_distance.mat', 'neibIdx','CDist');
% else
%     load('intial_distance.mat', 'neibIdx','CDist');
% end
validPts = ~cellfun(@isempty,CDist);
nearestDist = cellfun(@min,CDist(validPts));
% fit the gamma distribution using current overlapping distances
%phat = gamfit(nearestDist);
phat = fitTruncGamma(nearestDist);
movieInfo.ovGamma = phat;
Cij = cell(g.particleNum, 1);
CDist_i2j = CDist_dirwise; % !!! CDist_i2j: save both in two columns
for i=1:numel(CDist)
    if isempty(CDist{i})
        continue;
    end
    % overlapping_cost
    Cij{i} = overlap2cost(CDist{i}, phat);
    %     p_oc = 1-gamcdf(CDist{i}, phat(1), phat(2));
    %     oc = norminv(1-p_oc/2);
    %     Cij{i} = oc.^2;%
end

% cal the overlapping ratio forward and backward
ovSize = cell(numel(neibIdx),1);
for i=1:numel(neibIdx)
    ovSize{i} = zeros(length(neibIdx{i}),1);
    cur_voxIdxes = movieInfo.voxIdx{i};
    for j=1:length(ovSize{i})
        ovSize{i}(j) = length(intersect(cur_voxIdxes, ...
            movieInfo.voxIdx{neibIdx{i}(j)}));
    end
end
preNei = cell(numel(neibIdx),1);
preOvSize = cell(numel(neibIdx),1);
Cji = cell(numel(neibIdx),1);
CDist_j2i = cell(numel(neibIdx), 1);
for i=1:numel(neibIdx)
    for j=1:length(neibIdx{i})
        preNei{neibIdx{i}(j)} = cat(1, preNei{neibIdx{i}(j)}, i);
        preOvSize{neibIdx{i}(j)} = cat(1, preOvSize{neibIdx{i}(j)}, ovSize{i}(j));
        Cji{neibIdx{i}(j)} = cat(1, Cji{neibIdx{i}(j)}, Cij{i}(j));
        
        CDist_j2i{neibIdx{i}(j)} = cat(1, CDist_j2i{neibIdx{i}(j)}, ...
            CDist_i2j{i}(j,2));
    end
end

movieInfo.CDist = CDist;
movieInfo.CDist_i2j = CDist_i2j;% !!! CDist_i2j: save both i2j and j2i in two columns
movieInfo.CDist_j2i = CDist_j2i;

movieInfo.nei = neibIdx;
movieInfo.preNei = preNei;
movieInfo.ovSize = ovSize;
movieInfo.preOvSize = preOvSize;
movieInfo.Cij = Cij;
movieInfo.Cji = Cji;

fprintf('finish trajectory intialization with purely distance!\n');

end