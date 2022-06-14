function [one_cell_flag, oneCellpValue, mergeVSseg] = cellNumJudge(...
    left_cells, right_cells, movieInfo, refine_res, vidMap, g, ...
    eigMaps, q)
% left_cells: left tree with leaf ids in the same level in the same cell
% right_cells: right tree with leaf ids in the same level in the same cell

% March. 2020
left_cnt = cellfun(@length, left_cells); % both are start from root node
right_cnt = cellfun(@length, right_cells);
pVal_thres = 0.05;
%% way 1: the two trees lead to one conclusion
if false
    mergeVSseg = [0 0];
    mergeVSseg(1) = sum(left_cnt==1) + sum(right_cnt==1); % element 1 means one cell, 2 means >1 cells
    mergeVSseg(1) = mergeVSseg(1)-1;% root node was counted twice
    mergeVSseg(2) = sum(left_cnt>1) + sum(right_cnt>1);
    
    oneCellpValue = binocdf(mergeVSseg(1), sum(mergeVSseg), 0.5);
    

    if oneCellpValue < pVal_thres % probability of one cell in this region is low
        one_cell_flag = 0;
        return;
    end
    if 1-oneCellpValue < pVal_thres % probability of multi-cell in this region is low
        one_cell_flag = 1;
        return;
    end
end
%% way 2: the neighbors of the root node lead to consistent conclusion
% no promising conclusion can be made, we use cumulative evidence
% pVal_thres = 0.01;
mergeVSseg = [1 0]; % 1-vs-more
oneCellpValue = 0.5; % initialize as non-significant
for i=1:min(numel(left_cnt)-1, numel(right_cnt)-1)
    if left_cnt(i+1)>1
        mergeVSseg(2) = mergeVSseg(2) + 1;
    else
        mergeVSseg(1) = mergeVSseg(1) + 1;
    end
    
    if right_cnt(i+1)>1
        mergeVSseg(2) = mergeVSseg(2) + 1;
    else
        mergeVSseg(1) = mergeVSseg(1) + 1;
    end
    oneCellpValue = binocdf(mergeVSseg(1)-1, sum(mergeVSseg), 0.5);
    if oneCellpValue < pVal_thres || (1-oneCellpValue) < pVal_thres
        break;
    end
end
if oneCellpValue < pVal_thres % probability of one cell in this region is low
    one_cell_flag = 0;
    return;
end
if 1-oneCellpValue < pVal_thres % probability of multi-cell in this region is low
    one_cell_flag = 1;
    return;
end

%% way 3: at least one side can give out a conclusion
mergeVSseg = [1 0]; % 1-vs-more
leftpValue = 0.5; % initialize as non-significant
for i=1:numel(left_cnt)-1
    if left_cnt(i+1)>1
        mergeVSseg(2) = mergeVSseg(2) + 1;
    else
        mergeVSseg(1) = mergeVSseg(1) + 1;
    end
    leftpValue = binocdf(mergeVSseg(1)-1, sum(mergeVSseg), 0.5);
    if leftpValue < pVal_thres || (1-leftpValue) < pVal_thres
        break;
    end
end
mergeVSseg = [1 0]; % 1-vs-more
rightpValue = 0.5; % initialize as non-significant
for i=1:numel(right_cnt)-1
    if right_cnt(i+1)>1
        mergeVSseg(2) = mergeVSseg(2) + 1;
    else
        mergeVSseg(1) = mergeVSseg(1) + 1;
    end
    rightpValue = binocdf(mergeVSseg(1)-1, sum(mergeVSseg), 0.5);
    if rightpValue < pVal_thres || (1-rightpValue) < pVal_thres
        break;
    end
end
% if we have conclusion on at least one side, we did experiment
% if min(leftpValue, rightpValue) < pVal_thres ...
%         || max(leftpValue, rightpValue) > 1-pVal_thres
if true
    root_node = right_cells{1};
    bases = movieInfo.kids{root_node};
    if length(bases) ~= 2
        bases = movieInfo.parents{root_node};
    end
    if length(bases) ~= 2
        one_cell_flag = nan;
    else % do further test
        root_voxIdx = movieInfo.voxIdx{root_node};
        root_frame = movieInfo.frames(root_node);
        base_voxIdx = movieInfo.voxIdx(bases);
%         base_frame = movieInfo.frames(bases);
%         smooth_sc = [5 3];
%         cost_design = [1, 2];
%         connect = 6;        
%         shift = [3 3 1];% no need to enlarge the region
        regVoxIdx = bisectRegion(root_voxIdx, root_frame, base_voxIdx, ...
            vidMap, refine_res, eigMaps, q);
        
        frames = [movieInfo.frames(root_node), movieInfo.frames(bases(1))];
        [maxCost1, ~] = voxIdx2cost(regVoxIdx{1}, base_voxIdx{1}, frames, movieInfo, refine_res);
        frames = [movieInfo.frames(root_node), movieInfo.frames(bases(2))];
        [maxCost2, ~] = voxIdx2cost(regVoxIdx{2}, base_voxIdx{2}, frames, movieInfo, refine_res);
        
        if max(maxCost1, maxCost2) < abs(g.observationCost)
            one_cell_flag = 0;
            %fprintf('the hard case was split\n');
        else
            one_cell_flag = 1;
            %fprintf('the hard case was merge\n');
        end
    end
else
    one_cell_flag = nan;
    %fprintf('the hard case was left alone\n');
end
