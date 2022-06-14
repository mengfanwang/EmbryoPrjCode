function movieInfo = mergeBrokenTracksV2(movieInfo, g)
% there are some broken trajectories because region size/location 
% difference, we should link them

% V1: we re-grow the two regions at adjacnt frames
% V2 (this): 1: check jump 2: use distance rather than overlap to link
% broken tracks

movieInfo = relinkJumpOveredCell(movieInfo, g);
%% use distance rather than overlap to link broken tracks
% build distance distribution
intrack_distance = cell(numel(movieInfo.parents),1);
for i=1:numel(movieInfo.parents)
    if ~isempty(movieInfo.parents{i})% does not consider drift?
        xd = movieInfo.xCoord(movieInfo.parents{i})-movieInfo.xCoord(i);
        yd = movieInfo.yCoord(movieInfo.parents{i})-movieInfo.yCoord(i);
        zd = movieInfo.zCoord(movieInfo.parents{i})-movieInfo.zCoord(i);
        intrack_distance{i} = sqrt(xd.^2 + yd.^2 + zd.^2);
    end
end
intrack_distance = cat(1, intrack_distance{:});
phat = fitTruncGamma(intrack_distance);
% find track head and tails
track_heads = find(cellfun(@length, movieInfo.parents)==0);
track_heads(track_heads>length(movieInfo.xCoord)) = [];
track_tails = find(cellfun(@length, movieInfo.kids)==0);
track_tails(track_tails>length(movieInfo.xCoord)) = [];

valid_ids = union(track_heads, track_tails);
% these elements are enough for the tracklet based linking
movieInfo_trackLet=struct('xCoord',[],'yCoord',[], 'zCoord', [], ...
    'frames', [], 'CDist', [], 'Ci', [], 'nei',[], 'Cij',[], 'ovGamma', []);

movieInfo_trackLet.xCoord = movieInfo.xCoord(valid_ids);
movieInfo_trackLet.yCoord = movieInfo.yCoord(valid_ids);
movieInfo_trackLet.zCoord = movieInfo.zCoord(valid_ids);
movieInfo_trackLet.frames = movieInfo.frames(valid_ids);
movieInfo_trackLet.CDist = movieInfo.CDist(valid_ids);
movieInfo_trackLet.Ci = movieInfo.Ci(valid_ids);
movieInfo_trackLet.ovGamma = phat;

map = nan(length(movieInfo.xCoord), 1);
map(valid_ids) = 1:length(valid_ids);

movieInfo_trackLet.nei = cell(length(valid_ids), 1);
% only tails of tracks is valid for further linking (to heads)
movieInfo_trackLet.nei(map(track_tails)) = movieInfo.nei(track_tails);
% remove kids that not a head of a track
valid_kids_map = false(length(movieInfo.xCoord), 1);
valid_kids_map(track_heads) = true;
movieInfo_trackLet.nei = cellfun(@(x) x(valid_kids_map(x)), ...
    movieInfo_trackLet.nei,'UniformOutput', false);

movieInfo_trackLet.Cij = cell(length(valid_ids), 1);
for i=1:numel(movieInfo_trackLet.nei)
    if ~isempty(movieInfo_trackLet.nei{i})
        flag = false(length(movieInfo_trackLet.nei{i}),1);
        for k = 1:length(movieInfo_trackLet.nei{i})
            flag(k) = bestOvPair(movieInfo, valid_ids(i),...
                movieInfo_trackLet.nei{i}(k));
        end
        movieInfo_trackLet.nei{i} = movieInfo_trackLet.nei{i}(flag);
        if isempty(movieInfo_trackLet.nei{i})
           continue;
        end
        cur_nei = map(movieInfo_trackLet.nei{i});% change to new id
        cur_nei = cur_nei(~isnan(cur_nei));
        movieInfo_trackLet.nei{i} = cur_nei;
        if ~isempty(cur_nei)
            jump_punish = ...
                movieInfo.jumpRatio(movieInfo_trackLet.frames(cur_nei)-...
                movieInfo_trackLet.frames(i));
            xd = movieInfo_trackLet.xCoord(i)-movieInfo_trackLet.xCoord(cur_nei);
            yd = movieInfo_trackLet.yCoord(i)-movieInfo_trackLet.yCoord(cur_nei);
            zd = movieInfo_trackLet.zCoord(i)-movieInfo_trackLet.zCoord(cur_nei);
            movieInfo_trackLet.CDist{i} = sqrt(xd.^2 + yd.^2 + zd.^2);
            movieInfo_trackLet.Cij{i} = overlap2cost(...
                movieInfo_trackLet.CDist{i}, phat, jump_punish);
        end
    end
end
% mcc tracking
[~, g, dat_in] = trackGraphBuilder_cell(movieInfo_trackLet, g);
movieInfo_trackLet = mccTracker(dat_in, movieInfo_trackLet, g, false, false);

% update movieInfo using these new tracks
% including movieInfo.tracks; pathCost; particle2track; parents; kids
track_len = cellfun(@length, movieInfo_trackLet.tracks);
movieInfo_trackLet.tracks(track_len<2) = [];
for i=1:numel(movieInfo_trackLet.tracks)
    cur_track = movieInfo_trackLet.tracks{i};
    % parents and kids
    cur_track = valid_ids(cur_track); % retrieve real ids
    for j=1:length(cur_track)-1
        movieInfo.parents{cur_track(j+1)} = cur_track(j);
        movieInfo.kids{cur_track(j)} = cur_track(j+1);
    end
end
% union tracks for merging
% merged_tracks = cell(numel(movieInfo_trackLet.tracks), 1);
loopCnt = 0;
while true
    for j=1:numel(movieInfo_trackLet.tracks)
        if isempty(movieInfo_trackLet.tracks{j})
            continue;
        end
        track_ids = movieInfo.particle2track(...
            valid_ids(movieInfo_trackLet.tracks{j}),1);
        track_ids = track_ids(~isnan(track_ids));
        if isempty(track_ids)
            continue;
        end
        for k = j+1:numel(movieInfo_trackLet.tracks)
            %             if j==20 && k==28
            %                 keyboard;
            %             end
            if isempty(movieInfo_trackLet.tracks{k})
                continue;
            end
            track_ids_new = movieInfo.particle2track(...
            valid_ids(movieInfo_trackLet.tracks{k}),1);
            track_ids_new = track_ids_new(~isnan(track_ids_new));
            if ~isempty(intersect(track_ids, track_ids_new))
                movieInfo_trackLet.tracks{j} = union(...
                    movieInfo_trackLet.tracks{j}, ...
                    movieInfo_trackLet.tracks{k});
                movieInfo_trackLet.tracks{k} = [];
                track_ids = union(track_ids, track_ids_new);
            end
        end
    end
    track_len = cellfun(@length, movieInfo_trackLet.tracks);
    loopCnt = loopCnt + 1;
    if loopCnt > 3 || sum(track_len<2)==0
        break;
    end
    movieInfo_trackLet.tracks(track_len<2) = [];
end
rm_track_ids = cell(numel(movieInfo_trackLet.tracks), 1);
for i=1:numel(movieInfo_trackLet.tracks)
    cur_track = movieInfo_trackLet.tracks{i};
    cur_track = valid_ids(cur_track); % retrieve real ids
    % tracks and pathCost (cost is not accurate since it is not important)
    track_ids = movieInfo.particle2track(cur_track,1); 
    if ~isempty(find(track_ids == 154,1))
        disp(i);
    end
    tmp_track = [];
    tmpCost = 0;
    for j=1:length(track_ids)
        if isnan(track_ids(j))
            tmp_track = cat(1, tmp_track, cur_track(j));
        else
            if isempty(find(track_ids(1:j-1)==track_ids(j), 1))
                tmp_track = cat(1, tmp_track, ...
                    movieInfo.tracks{track_ids(j)});
                tmpCost = tmpCost + movieInfo.pathCost(track_ids(j));
            end
        end
    end
    tmp_track = sort(tmp_track, 'ascend');
    track_ids = unique(track_ids);
    keep_id = find(~isnan(track_ids),1);
    if ~isempty(keep_id)
        movieInfo.tracks{track_ids(keep_id)} = tmp_track;
        movieInfo.pathCost(track_ids(keep_id)) = tmpCost;
        cur_rm_track_ids = track_ids;
        cur_rm_track_ids(keep_id) = [];
        cur_rm_track_ids(isnan(cur_rm_track_ids)) = [];
    else
        movieInfo.tracks = cat(1, movieInfo.tracks, tmp_track);
        movieInfo.pathCost = cat(1, movieInfo.pathCost, tmpCost);
        cur_rm_track_ids = [];
    end
    
    if ~isempty(cur_rm_track_ids)
        rm_track_ids{i} = cur_rm_track_ids;
    end
end
rm_track_ids = unique(cat(1, rm_track_ids{:}));
movieInfo.tracks(rm_track_ids) = [];
movieInfo.pathCost(rm_track_ids) = [];

% particle2track
% 1: the track it belongs to; 2: its position
particle2track = nan(length(movieInfo.xCoord), 3);
for j=1:numel(movieInfo.tracks)
    particle2track(movieInfo.tracks{j},1:2) = ...
        [j+zeros(length(movieInfo.tracks{j}),1),...
        [1:length(movieInfo.tracks{j})]'];
end
movieInfo.particle2track = particle2track;

%% we should update the regions that linked by this function


