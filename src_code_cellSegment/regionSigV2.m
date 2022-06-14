function [zscore,z_debias, sz_nei] = regionSigV2(regMap, vidComp, fmap, validNeiMap, OrSt)
%GETORDERSTAT Compute order statistics
% d: transformed data. smask: synapse mask.
% bmask: background mask, MUST be connected. imgVar: image variance

[fgBgCells, fgBgLocCells] = localizeFgBgPairs(vidComp, regMap, ...
    validNeiMap, 3, 20, 1);
sz = cellfun(@length, fgBgCells);
sz_nei = sum(sz(:,2));
% we need to estimate variance sample wise
if OrSt.fgVarSamplewiseEst
    if strcmp(OrSt.imProcMethod,'stb')
        varMap = OrSt.stbVarCropMap;
    else
        varMap = OrSt.NoStbVarCropMap;
    end
end
zscore = nan(size(fgBgCells,1), 1);
z_debias = nan(size(fgBgCells,1), 1);
for i=1:size(fgBgCells,1)
    fg_locs = fgBgLocCells{i,1};
    nei_locs = fgBgLocCells{i,2};
    fg = fgBgCells{i,1};
    fg_neighbors = fgBgCells{i,2};
    % approximation of mu and sigma
    M = max(length(fg), 10);
    N = max(length(fg_neighbors), 10);
    
    if M>10*N
        zscore(i) = nan;
        z_debias(i) = nan;
        sz_nei = nan;
    else
        if OrSt.fgVarSamplewiseEst
            OrSt.curStbVar = nanmean(varMap(cat(1, fg_locs, nei_locs))); %
            %disp(sqrt(OrSt.curStbVar));
        end
        nanVec = [];
        [zscore(i),z_debias(i)] = ordStats(fg, fg_neighbors, nanVec, OrSt);
    end
    %
    % disp([mu sigma]);
    % zscore2 = (sum_st - mu*sqrt(OrSt.NoStbVar))...
    %     / (sigma*sqrt(OrSt.NoStbVar));
end

zscore = nanmean(zscore);
z_debias = nanmean(z_debias);