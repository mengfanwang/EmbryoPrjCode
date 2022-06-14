function [movieInfo,refine_res, thresholdMaps, upt_ids] = ...
    testMissingCellGotFromOneSeed(comMaps, pseudo_seed_label, ...
    parent_kid_vec, missed_frame, movieInfo,refine_res, thresholdMaps,...
    embryo_vid,eigMaps, g, q)
% if we detected a cell with a given seed, we need to test its relationship
% with its surrounding cells. Should we merge them?

upt_ids = [];
idsAndLocs = [];
if q.updateCellsAdjMissingCell
    idsAndLocs = extractCellLocsInMovieInfo(comMaps.idComp, ...
        movieInfo, missed_frame);
    if isempty(idsAndLocs) || isempty(find(idsAndLocs(:,1)==pseudo_seed_label,1))
        error('No added cell at all!');
    else
        idsAndLocs(idsAndLocs(:,1)==pseudo_seed_label,:) = [];
    end
    [movieInfo, refine_res, thresholdMaps, voxIdxCells] = addVox2adjCells(...
        comMaps, idsAndLocs, movieInfo, refine_res, thresholdMaps);
end
% if the new region adjacent to an existing cell
% should we consider z-direction???
cur_nei_regions = comMaps.idComp(...
    imdilate(comMaps.regComp, strel('disk', 1)) & ~comMaps.regComp);

cur_regions_cnt = frequency_cnt(cur_nei_regions(cur_nei_regions>0 ...
    & cur_nei_regions~=pseudo_seed_label));
if ~isempty(find(cur_regions_cnt(:,1)~=0,1)) % there indeed an adjacent cell, merge them
    [~, od] = max(cur_regions_cnt(:,2));
    adj_cell = cur_regions_cnt(od,1);
    if adj_cell <= movieInfo.n_perframe(missed_frame)
        adj_cell_id = sum(movieInfo.n_perframe(1:missed_frame-1)) ...
            + adj_cell;
    else
        adj_cell_id = adj_cell;% this is a new region whose label is
        % just its id in movieInfo
    end
    if isempty(movieInfo.voxIdx{adj_cell_id})
        adj_cell_id = [];
        
    elseif refine_res{missed_frame}(movieInfo.voxIdx{adj_cell_id}(1))...
            ~= adj_cell
        error('The adjacent cell is wrong!');
    end
else % there no adjacent cell, then create a new cell
    adj_cell_id = [];
end
append_idx = comMaps.linerInd(comMaps.fmapComp);

if g.blindmergeNewCell%!!! if we blindly merge
    if isempty(adj_cell_id)
        flag = [1 1];
    else
        flag = [2 1];
    end
else
    [flag, newSplitIdx] = newlyAddedCellValidTest(append_idx, missed_frame, ...
        adj_cell_id, parent_kid_vec, movieInfo, refine_res, embryo_vid,...
        eigMaps, comMaps, g, q);
end

if flag(2) ~= 0 && ...
        (~q.removeSamllRegion || length(append_idx) >= q.minSize)
    % flag(1) = [2 3 5] ==> merge, for 3, we can split now
    % flag(1) = [1 4] ==> split
    if ismember(flag(1), [2 5])
        movieInfo.voxIdx{adj_cell_id} = cat(1, ...
            movieInfo.voxIdx{adj_cell_id}, append_idx);
        refine_res{missed_frame}(append_idx) = adj_cell;
        upt_ids = cat(1, upt_ids, adj_cell_id);
    else
        if flag(1) == 3 % update current node/newly added node
            if isempty(newSplitIdx)
                error('newSplitIdx not assigned!');
            end
            append_idx = newSplitIdx{2};
            old_idx = movieInfo.voxIdx{adj_cell_id};
            refine_res{missed_frame}(old_idx) = 0;
            refine_res{missed_frame}(newSplitIdx{1}) = adj_cell;
            old_thres = thresholdMaps{missed_frame}(old_idx(1));
            if old_thres == 0
                error('threshold not assigned!');
            end
            thresholdMaps{missed_frame}(old_idx) = 0;
            thresholdMaps{missed_frame}(newSplitIdx{1}) = old_thres;
            movieInfo.voxIdx{adj_cell_id} = newSplitIdx{1};
            upt_ids = cat(1, upt_ids, adj_cell_id);
        elseif flag(1) == 6 % update the index of the newly added node
            if isempty(newSplitIdx)
                error('newSplitIdx not assigned!');
            end
            append_idx = newSplitIdx;
        end
        % case 1 and 4, we can directly add the new cell
        movieInfo = addNewCell2MovieInfo(movieInfo, ...
            append_idx, parent_kid_vec, missed_frame);
        refine_res{missed_frame}(append_idx) = pseudo_seed_label;
        upt_ids = cat(1, upt_ids, pseudo_seed_label);
    end
    thresholdMaps{missed_frame}(append_idx) = ...
        comMaps.pickedThreshold;
    if ~isempty(idsAndLocs)
        upt_ids = cat(1, upt_ids, idsAndLocs(:,2));
    end
else % the new cell is invalid, continue to next, so roll back
    for j=1:size(idsAndLocs, 1)
        if ~isempty(voxIdxCells{j,2})
            movieInfo.voxIdx{idsAndLocs(j,2)} = voxIdxCells{j,1};
            refine_res{missed_frame}(voxIdxCells{j,2}) = 0;
            thresholdMaps{missed_frame}(voxIdxCells{j,2}) = 0;
        end
    end
    %         if ~isempty(rm_cell_id_idx_p_k) % a cell is removed, it also should roll back
    %             movieInfo.voxIdx(rm_cell_id_idx_p_k{1}(:,1)) = rm_cell_id_idx_p_k{2};
    %             movieInfo.parents(rm_cell_id_idx_p_k{1}(:,1)) = rm_cell_id_idx_p_k{3};
    %             movieInfo.kids(rm_cell_id_idx_p_k{1}(:,1)) = rm_cell_id_idx_p_k{4};
    %             for j=1:size(rm_cell_id_idx_p_k{1},1)
    %                 refine_res{missed_frame}(rm_cell_id_idx_p_k{2}{j}) = ...
    %                     rm_cell_id_idx_p_k{1}(j,2);
    %                 thresholdMaps{missed_frame}(rm_cell_id_idx_p_k{2}{j}) = ...
    %                     rm_cell_id_idx_p_k{1}(j,3);
    %             end
    %             if length(upt_ids) == size(rm_cell_id_idx_p_k{1},1)
    %                 upt_ids = []; % no need to update any more
    %             else
    %                 error('some other ids for checking, need carefully re-check');
    %             end
    %         end
end
end