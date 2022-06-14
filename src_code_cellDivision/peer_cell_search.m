function [parent, peer_node] = peer_cell_search(current_candidate_id, ...
    division_candidate_ids, movieInfo, refine_res, embryo_vid, g)
% We search the sourrouding area of the cell to see if a peer cell exist.
% Such peer cell should be either a candidate too or has a meaningful
% parent.
% We move these two cells together to their center, they jointly should be
% a good region to link to their parent. Their common parent should
% satisify that
% 1) it is a parent of one of the parent or end of a trace
% 2) the cost of linking it to the joint kids are better than previous
% linking cost
    

current_candidate_xyz = movieInfo.vox{current_candidate_id};

candidate_f = movieInfo.frames(current_candidate_id);

id_map = crop3D(refine_res{candidate_f}, movieInfo.voxIdx{current_candidate_id},...
    g.maxDistXYZ);

reg_map = id_map==refine_res{candidate_f}(movieInfo.voxIdx{current_candidate_id}(1));

candidate_peer_real_ids = unique(id_map(~reg_map & id_map > 0));

adj_peer_real_ids = unique(id_map(imdilate(reg_map,strel('sphere', 1)) & ...
    ~reg_map & id_map > 0));

if candidate_f > 1
    candidate_peers_id = sum(movieInfo.n_perframe(1:candidate_f-1)) + ...
        candidate_peer_real_ids;
    adj_peer_ids = sum(movieInfo.n_perframe(1:candidate_f-1)) + ...
        adj_peer_real_ids;
else
    candidate_peers_id = candidate_peer_real_ids;
    adj_peer_ids = adj_peer_real_ids;
end




priors = ismember(candidate_peers_id, division_candidate_ids);
prior_ids = candidate_peers_id(priors);
candidate_peers_id(priors) = [];
candidate_peers_id = cat(1, prior_ids, candidate_peers_id);

currrent_candidate_center_xyz = mean(current_candidate_xyz,1);
bound_xyz = size(refine_res{1});
bound_xyz = bound_xyz([2 1 3]);
parent = nan;
peer_node = nan;
se = strel('sphere', 3);
nhood = se.Neighborhood(:,:,3:5);
for i=1:length(candidate_peers_id)
    if ~isempty(find(adj_peer_ids == candidate_peers_id(i),1))
        adj_flag = true;
        k1 = movieInfo.kids{current_candidate_id};
        k2 = movieInfo.kids{candidate_peers_id(i)};
        split_cnt = 0;
        while true
            if length(k1)~=1 || length(k2)~=1
                break;
            end
            if movieInfo.frames(k1) < movieInfo.frames(k2)
                k1 = movieInfo.kids{k1};
                continue;
            end
            if movieInfo.frames(k1) > movieInfo.frames(k2)
                k2 = movieInfo.kids{k2};
                continue;
            end

            adj_flag = regionAdjacent(movieInfo, refine_res, k1, k2, nhood);
            if ~adj_flag
                if split_cnt == 0
                    start_split_pair = [k1, k2];
                end
                split_cnt = split_cnt + 1;
                if split_cnt == 5
                    break;
                end
            elseif split_cnt > 0
                break;
            end
            k1 = movieInfo.kids{k1};
            k2 = movieInfo.kids{k2};
        end
        if adj_flag
            continue;
        end
    else
        start_split_pair = [current_candidate_id, candidate_peers_id(i)];
    end
    sz_ratio = length(movieInfo.voxIdx{start_split_pair(1)}) ...
        / length(movieInfo.voxIdx{start_split_pair(2)});
    if sz_ratio < 0.5 || sz_ratio > 1.5
        continue;
    end
    peer_center_xyz = mean(movieInfo.vox{candidate_peers_id(i)}, 1);
    joint_center_xyz = (currrent_candidate_center_xyz + peer_center_xyz) ./2;

    if ~isempty(movieInfo.parents{candidate_peers_id(i)}) 
        % there should be only one parent
        test_parents = movieInfo.parents{candidate_peers_id(i)}; 
    else
        test_parents = search_tail_from_center(joint_center_xyz, candidate_f, ...
            movieInfo, refine_res, g);
    end
    if isempty(test_parents)
        continue;
    end
    if ~isempty(find(adj_peer_ids == candidate_peers_id(i),1))
        curr_can_shift_xyz = current_candidate_xyz;
        peer_shift_xyz = movieInfo.vox{candidate_peers_id(i)};
    else
        curr_can_shift_xyz = shift2touch_center(joint_center_xyz, ...
            current_candidate_xyz, bound_xyz);
        peer_shift_xyz = shift2touch_center(joint_center_xyz, ...
            movieInfo.vox{candidate_peers_id(i)}, bound_xyz);
    end
    valid_parent = nan;
    for j = 1:length(test_parents) % sort already by time (reversed)
        link_cost = voxIdx2cost(cat(1, curr_can_shift_xyz, peer_shift_xyz), ...
             movieInfo.vox{test_parents(j)}, ...
            [candidate_f movieInfo.frames(test_parents(j))], ...
            movieInfo, size(refine_res{1}), movieInfo.jumpRatio);
        max_cost = inf;
        if length(movieInfo.kids{test_parents(j)}) == 1% should be only one kid
            max_cost = movieInfo.Cij{test_parents(j)}(...
                movieInfo.nei{test_parents(j)} == movieInfo.kids{test_parents(j)});
        elseif length(movieInfo.kids{test_parents(j)}) > 1
            max_cost = -inf;
        end
        
        if link_cost < min(abs(g.c_ex), max_cost)
            disp(link_cost);
            valid_parent = test_parents(j);
            break;
        end
    end
    if ~isnan(valid_parent)
        parent = valid_parent;
        peer_node = candidate_peers_id(i);
        break;
    end
end