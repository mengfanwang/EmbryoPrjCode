function [refine_res, movieInfo, check_kids] = mergeBrokenTracks(movieInfo, g, refine_res, embryo_vid)
% there are some broken trajectories because region size/location 
% difference, we should link them
% Test the end node of each track (with heads of other tracks)
% Both cases, the node has a unique best choice to link

% V1 (this): we re-grow the two regions at adjacnt frames
% V2: we use distance rather than overlap to link them

num_frames = numel(refine_res);
upt_voxIdx = cell(numel(movieInfo.xCoord), 1);
iou_score = zeros(numel(movieInfo.xCoord), 1);
parentKidVec = nan(numel(movieInfo.xCoord)*4, 2);
pkCnt = 0;
grCnt = 0;
fprintf('Start to processing broken tracks\n');
for i=1:numel(movieInfo.xCoord)
    if mod(i, 100) == 0
        disp(i);
    end
    if movieInfo.frames(i) == num_frames
        continue;
    end
    cur_track_id = movieInfo.particle2track(i,1);
    if ~isnan(cur_track_id)
        cur_track = movieInfo.tracks{cur_track_id};
%         cur_locs = movieInfo.particle2track(i,2);
        if movieInfo.frames(i) ~= movieInfo.frames(cur_track(end))
            continue;
        end
    end
    % find the head of another track (or isolated cell) in next frame
    nei = movieInfo.nei{i};
    nei = nei(movieInfo.frames(nei) == (movieInfo.frames(i)+1));% may not necessary
    locsInTracks = movieInfo.particle2track(nei, 1);
    valid_nei = false(length(nei), 1);
    for j=1:length(nei)
        if isnan(locsInTracks(j))
            valid_nei(j) = true;
        else
            if isempty(movieInfo.parents{nei(j)})
                valid_nei(j) = true;
            end
%             node1st = movieInfo.tracks{locsInTracks(j)}(1);
%             frame1st = movieInfo.frames(node1st);
%             if movieInfo.frames(nei(j)) == frame1st
%                 valid_nei(j) = true;
%             end
        end
    end
    candidateKids = nei(valid_nei);
    if isempty(candidateKids)
        continue;
    end
    % if the head of another track and current node form the best pair
    kidsFrames = movieInfo.frames(candidateKids);
    frs = unique(kidsFrames);% currently, we have test only one frame
    for j=1:length(frs)
        valid_kids = candidateKids(kidsFrames==frs(j));
        flag = false(length(valid_kids), 1);
        for k = 1:length(valid_kids)
            [max_parent, ~] = bestOvNei(valid_kids(k), movieInfo, ...
                movieInfo.frames(i), 1);
            if max_parent == i
                flag(k) = true;
            end
        end
        valid_kids = valid_kids(flag);
        % if want add best kids to test array, uncomment these 4 lines
        %         [max_kid, ~] = bestOvNei(i, movieInfo, frs(j), 0);
        %         if isempty(find(valid_kids == max_kid, 1))
        %             valid_kids = cat(1, valid_kids, max_kid);
        %         end
        parent_id = i;
        for k=1:length(valid_kids)
            kid_id = valid_kids(k);
            [p_k_vec, minmaxIds, minIdVoxIdx, maxIdVoxIdx] = ...
                segmentErrorTest(parent_id, kid_id, g, movieInfo, ...
                refine_res, embryo_vid);
            minId = minmaxIds(1);
            maxId = minmaxIds(2);
            if ~isempty(p_k_vec)
                pkCnt = pkCnt + 1;
                parentKidVec(pkCnt,:) = p_k_vec;
            elseif ~isempty(minIdVoxIdx)
%                 if ~isempty(upt_voxIdx{minId})
%                     keyboard;
%                 end
                
                grCnt = grCnt + 1;
                upt_voxIdx{minId} = minIdVoxIdx;
                upt_voxIdx{maxId} = maxIdVoxIdx;
                if minId == parent_id
                    in_l = length(intersect(minIdVoxIdx, movieInfo.voxIdx{kid_id}));
                    un_l = length(minIdVoxIdx) + length(movieInfo.voxIdx{kid_id})-in_l;
                else
                    in_l = length(intersect(minIdVoxIdx, movieInfo.voxIdx{parent_id}));
                    un_l = length(minIdVoxIdx) + length(movieInfo.voxIdx{parent_id})-in_l;
                end
                iou_score(minId) = in_l/un_l;
                %iou_score(minId) = in_l/un_l;
            end
        end
    end
end

% merge regions
parentKidVec = parentKidVec(1:pkCnt,:);
[p_ids, od] = sort(parentKidVec(:,1));
parentKidVec = parentKidVec(od,:);
endPts = find(p_ids(2:end) - p_ids(1:end-1));
stPts = [1; endPts+1];
endPts = [endPts; length(p_ids)];
check_kids = [];
for i=1:length(stPts)
    kids = parentKidVec(stPts(i):endPts(i), 2);
    processed = find(iou_score(kids)>0,1);
    if ~isempty(processed)
        continue;
    end
    parent = parentKidVec(stPts(i), 1);
    [voxIdx_newKid, ~] = growHighOvReg(parent, kids, movieInfo, refine_res, embryo_vid);
    in_l = length(intersect(voxIdx_newKid, movieInfo.voxIdx{parent}));
    un_l = length(minIdVoxIdx) + length(movieInfo.voxIdx{parent}) - in_l;
    iou = in_l/un_l;
%     processed = find(cellfun(@length, upt_voxIdx(kids))>0);
%     if ~isempty(processed)
%         if ~isempty(find(iou_score(processed) > iou, 1))
%             continue;
%         else
%             
%         end
%     end
    iou_score(kids) = iou;
    grCnt = grCnt + 1;
    upt_voxIdx{kids(1)} = voxIdx_newKid;
    if length(kids) > 1
        check_kids = cat(1, check_kids, kids(2:end));
    end
end
% replace updated regions
for i=1:length(check_kids)
    movieInfo.voxIdx{check_kids(i)} = [];
    % there can be regions both exist in parentKidVec(:,2) (kids) and 
    % upt_voxIdx. This should be avoied, cause we handle one condition each time.
    upt_voxIdx{check_kids(i)} = [];
end
valid_ids = find(cellfun(@length, upt_voxIdx) > 0);
for i=1:length(valid_ids)
    frame = movieInfo.frames(valid_ids(i));
    realId = refine_res{frame}(movieInfo.voxIdx{valid_ids(i)}(1));
    if realId == 0
        keyboard;
    end
    refine_res{frame}(movieInfo.voxIdx{valid_ids(i)}) = 0;
    refine_res{frame}(upt_voxIdx{valid_ids(i)}) = realId;
end
movieInfo.voxIdx(valid_ids) = upt_voxIdx(valid_ids);

check_kids = cat(1, check_kids, valid_ids);

fprintf('Done with broken tracks with %d growed regions\n', length(valid_ids));
end