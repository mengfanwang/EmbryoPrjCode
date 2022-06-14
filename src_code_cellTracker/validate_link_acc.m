function accuracy = validate_link_acc(movieInfo, label_maps, gt_tracks, movieInfo_TGMM)

num_frames = numel(movieInfo.validGapMaps);
%if isempty(label_maps) % generally this would
    % build label map first
    label_maps = cell(num_frames, 1);
    for i=1:num_frames
        label_maps{i} = uint16(zeros(size(movieInfo.validGapMaps{1})));
    end
    for i=1:numel(movieInfo.voxIdx)
        label_maps{movieInfo.frames(i)}(movieInfo.voxIdx{i}) = i;
    end
%end
accuracy = zeros(numel(gt_tracks), 5); % correct one, track length, ratio
if nargin == 3 % measure our tracking accuracy
    %used_id = [];
    for i=1:numel(gt_tracks)
        cur_track = gt_tracks{i};
        track_length = size(cur_track,1);
        p_id = [];
        for j=1:track_length - 1
            if isempty(p_id)
                pt = cur_track(j,:);
                p_id = point2detection_id(pt,label_maps);
            end
            kt = cur_track(j+1,:);
            k_id = point2detection_id(kt,label_maps);
            if ~isnan(p_id) && ~isnan(k_id)
                if ~isempty(find(movieInfo.kids{p_id}==k_id,1))
                    accuracy(i,1) = accuracy(i,1) + 1;
                    %used_id = cat(1, used_id, p_id, k_id);
                else
                    accuracy(i,4) = accuracy(i,4) + 1;
                end
            else
                accuracy(i,5) = accuracy(i,5) + 1;
            end
            p_id = k_id;
        end
        accuracy(i,2) = track_length-1;
    end
    accuracy(:,3) = accuracy(:,1) ./ accuracy(:,2);
else
    % build relationship between movieInfo_TGMM and movieInfor
    map2TGMM = cell(numel(movieInfo.voxIdx), 1);
    for i=1:size(movieInfo_TGMM.orgCoord, 1)
        cur_point = movieInfo_TGMM.orgCoord(i,:);
        cur_point = [movieInfo_TGMM.frames(i), ...
            cur_point(2), cur_point(1), cur_point(3)];
        [cor_id, cor_ids] = point2detection_id(cur_point,label_maps);
        if ~isnan(cor_id)
            val_ids = unique(cor_ids);
            for j=1:length(val_ids)
                map2TGMM{val_ids(j)} = cat(1, map2TGMM{val_ids(j)},...
                    i);
            end
        end
    end
    for i=1:numel(gt_tracks)
        cur_track = gt_tracks{i};
        track_length = size(cur_track,1);
        p_id = [];
        for j=1:track_length - 1
            if isempty(p_id)
                pt = cur_track(j,:);
                p_id = point2detection_id(pt,label_maps);
            end
            kt = cur_track(j+1,:);
            k_id = point2detection_id(kt,label_maps);
            if ~isnan(p_id) && ~isnan(k_id)
                TGMM_pids = map2TGMM{p_id};
                TGMM_kids = map2TGMM{k_id};
                if ~isempty(TGMM_pids) && ~isempty(TGMM_kids)
                    all_kids = cat(1, movieInfo_TGMM.kids{TGMM_pids});
                    if ~isempty(intersect(all_kids, TGMM_kids))
                        accuracy(i,1) = accuracy(i,1) + 1;
                        %used_id = cat(1, used_id, p_id, k_id);
                    else
                        accuracy(i,4) = accuracy(i,4) + 1;
                    end
                else
                    accuracy(i,5) = accuracy(i,5) + 1;
                end
%             else % if I have detection error, directly view TGMM as correct
%                 accuracy(i,1) = accuracy(i,1) + 1;
            end
            p_id = k_id;
        end
        accuracy(i,2) = track_length-1;
    end
    accuracy(:,3) = accuracy(:,1) ./ accuracy(:,2);
    
end

end
