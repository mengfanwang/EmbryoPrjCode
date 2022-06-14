function [combined_tracks, trajectories, parents, kids] = recover_paths_ssp_cell(res_G, numTracks, sink_node)
% parents: parent nodes of the current nodes
% kids: kid nodes of the current node

% recover the tracks
trajectories = cell(numTracks,1);
%stNodes = res_G.Edges.EndNodes(res_G.Edges.EndNodes(:,1)==sink_node,2);
[~, stNodes] = outedges(res_G,sink_node);
for i=1:numTracks
    %fprintf('processing track id %d out of %d tracks\n', i, numTracks);
    tail = stNodes(i);
%     if floor(tail/2) == 5145
%         keyboard;
%     end
%     res_G = flipedge(res_G, sink_node, tail);
    trajectories{i} = tail;
    %head = res_G.Edges.EndNodes(res_G.Edges.EndNodes(:,1) == tail, 2);
    [~, head] = outedges(res_G,tail);
    if mod(tail, 2) > 0
        head = head(head==1 | head == tail-1);
    else
        head(head==sink_node) = [];
    end
%     head(head>tail) = [];
    head = head(1);
    while head > 1
        trajectories{i} = [head, trajectories{i}];
%         if tail==2576 && head==2325
%             keyboard;
%         end
        %res_G = flipedge(res_G, tail, head);
        res_G = rmedge(res_G, tail, head);
        tail = head;
        %head = res_G.Edges.EndNodes(res_G.Edges.EndNodes(:,1) == tail, 2);
        [~, head]= outedges(res_G,tail);
        if mod(tail, 2) > 0
            head = head(head==1 | head == tail-1);
        else
            head(head==sink_node) = [];
        end
        head = head(1);
%         res_G = flipedge(res_G, tail, head);

    end
    if head == 1
        %res_G = flipedge(res_G, tail, head);
        res_G = rmedge(res_G, tail, head);
    end
%     if isempty(trajectories{i})
%         res_G = flipedge(res_G, tail, head);
%     end
    trajectories{i} = unique(floor(trajectories{i}/2), 'stable');
%     if sum(trajectories{i} == 1810)>0
%         keyboard;
%     end
end

cell_num = floor(sink_node/2) - 1;
parents = cell(cell_num,1);
kids = cell(cell_num,1);
cell2track = zeros(cell_num,1);
track_cnt = 0;
combined_tracks = cell(numTracks, 1);
for i=1:numTracks
    % find the parent nodes and kid nodes for each cell
    for j=1:numel(trajectories{i})-1
        kids{trajectories{i}(j)} = cat(1, kids{trajectories{i}(j)}, trajectories{i}(j+1));
        parents{trajectories{i}(j+1)} = cat(1, parents{trajectories{i}(j+1)}, trajectories{i}(j));
    end
    % merge overlapped tracks
    pre_track_ids = cell2track(trajectories{i});
    pre_track_ids = unique(pre_track_ids(pre_track_ids~=0));
%     if find(trajectories{i}>=1606 & trajectories{i}<=1608,1)
%         keyboard;
%     end
    if isempty(pre_track_ids) % new track
        track_cnt = track_cnt + 1;
        cell2track(trajectories{i}) = track_cnt;
        combined_tracks{track_cnt} = trajectories{i};
    elseif length(pre_track_ids) == 1
        cell2track(trajectories{i}) = pre_track_ids;
        combined_tracks{pre_track_ids} = cat(2, combined_tracks{pre_track_ids}, ...
            trajectories{i});
    else
        for j=2:length(pre_track_ids)
            combined_tracks{pre_track_ids(1)} = cat(2, combined_tracks{pre_track_ids(1)}, ...
            combined_tracks{pre_track_ids(j)});
            cell2track(combined_tracks{pre_track_ids(j)}) = pre_track_ids(1);
            combined_tracks{pre_track_ids(j)} = [];
        end
        combined_tracks{pre_track_ids(1)} = cat(2, combined_tracks{pre_track_ids(1)}, ...
            trajectories{i});
        cell2track(trajectories{i}) = pre_track_ids(1);
    end
    
end

% sanity check of parents and kids
for i=1:cell_num
    st_p = sort(parents{i});
    st_k = sort(kids{i});
    if length(st_p)>1
        if ~isempty(find(st_p(2:end)-st_p(1:end-1)==0, 1))
            keyboard;
        end
    end
    if length(st_k)>1
        if ~isempty(find(st_k(2:end)-st_k(1:end-1)==0, 1))
            keyboard;
        end
    end
end
combined_tracks = combined_tracks(1:track_cnt);
tr_len = cellfun(@length, combined_tracks);
combined_tracks = combined_tracks(tr_len>0);
for i=1:numel(combined_tracks)
    combined_tracks{i} = unique(sort(combined_tracks{i}, 'ascend'));
%     if ~isempty(find(combined_tracks{i}(2:end) - ...
%             combined_tracks{i}(1:end-1)==0,1))
%         keyboard;
%     end
end

fprintf('after merging, totally get %d tracks\n', numel(combined_tracks));
end