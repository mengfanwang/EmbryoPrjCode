function movieInfo = tree2tracks(detections, track_flag)
% change the tracking results from saving in a tree structure to isolated
% form, where each track are saved in a single cell
% INPUT:
% detections: t-by-1 cells, each cell corresponds to detections in one
% frame. For each cell, if it is n*3, each row is the yxz coordicates of
% one detection. If it is n*4, the last column is the parent id in the
% previous frame, which means tracking results is attained.
% track_flag: do we need to re-format the trajectories? Default is true(yes).
% OUTPUT:
% movieInfo: a struct consist of all the infor needed for our own tracking
% framework such as: y,x,z coordinate, frame stamps. In movieInfo, all the
% detections are in one matrix rather than cells as in the input
% "detections".

if nargin < 2
    track_flag = true;
end
movieInfo=struct('xCoord',[],'yCoord',[], 'zCoord', [], ...
    'frames', [],  'tracks',[], 'particle2track', []);
% coordinates
track_num = 0;
for i=1:numel(detections)
    detections{i} = cat(2, detections{i}, i*ones(size(detections{i},1),1));
    track_num = track_num + sum(isnan(detections{i}(:,4)));
end
det_all = cat(1, detections{:});
movieInfo.yCoord = det_all(:,1);
movieInfo.xCoord = det_all(:,2);
movieInfo.zCoord = det_all(:,3);
movieInfo.frames = det_all(:,end); % time
if size(detections{1},2) > 4 && track_flag % there is track infor
    % tracks
    tracks = cell(track_num, 1);
    det_cnt = 0;
    cnt = 0;
    for f=1:numel(detections)
        det_ids = (1:size(detections{f},1))' + det_cnt;
        % 4:parent, 5:frames, 6:track_id, 7: det_id
        detections{f} = cat(2, detections{f}, nan(size(detections{f},1), 1), det_ids);
        
        for i=1:size(detections{f},1)
            cur_pt = detections{f}(i,:);
            if isnan(cur_pt(4)) % new trajectory
                cnt = cnt + 1;
                detections{f}(i,6) = cnt;
                tracks{cnt} = [det_ids(i), nan];
            else % existing trajectory
                parent_pt = detections{f-1}(cur_pt(4), :);
                tr_id = parent_pt(6);
                detections{f}(i,6) = tr_id;
                tracks{tr_id} = cat(1, tracks{tr_id}, [det_ids(i), parent_pt(7)]);
            end
        end
        det_cnt = det_cnt + size(detections{f},1);
    end
    
    if cnt ~= track_num
        error('Wrong number of trajectories!');
    end
    movieInfo.tracks = tracks;
    det_all = cat(1, detections{:});
    movieInfo.particle2track = det_all(:,6);
end