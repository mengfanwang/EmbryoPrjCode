function [same_flag, cost] = verifySameCell(curVox, neiVox, prePostVox, phat, threshold)
% tell is two regions are belonging to the same cell or cell cluster
% curVox, neiVox: regions to be test
% phat: gamma distribution parameters
% threshold: for telling the significance of the cost

% 
if isempty(prePostVox{1}) && isempty(prePostVox{2})
    shift_xyz = mean(curVox)-mean(neiVox);
else
    if isempty(prePostVox{2})
        shift_xyz = mean(prePostVox{1}) - mean(neiVox);
    elseif isempty(prePostVox{1})
        shift_xyz = mean(prePostVox{2}) - mean(neiVox);
    else
        shift_xyz = (mean(prePostVox{1}) + mean(prePostVox{2}))/2 ...
        - mean(neiVox);
    end
    
end


distance = ovDistanceRegion(curVox, neiVox + shift_xyz);

cost = overlap2cost(distance, phat);
same_flag = cost <= threshold;
end