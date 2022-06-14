function phatDist = getStdFromTracks_cell(movieInfo, g)
% get the standard deviation from existing tracks
numTrajectories = numel(movieInfo.tracks);
validTrack = 0;
devDist = cell(numTrajectories,1);

for i=1:numTrajectories
    curTrack = movieInfo.tracks{i};
    if length(curTrack)<g.trackLength4var % only use those trajectories >=g.trackLength4var
        continue;
    end
    validTrack = validTrack+1;
    % way 1 of estimate variance
    cur_dist = zeros(length(curTrack)-1, 1);
    
    for j=1:length(curTrack)-1
        cur_dist(j) = movieInfo.CDist{curTrack(j)}...
            ((movieInfo.nei{curTrack(j)}==curTrack(j+1)));
    end
    % way 2: consider edges associated with over-segmented regions all as correct linking
%     cur_dist = zeros(2*length(curTrack), 1);
%     cd_cnt = 0;
%     cur_frames = movieInfo.frames(curTrack);
%     uni_frs = sort(unique(cur_frames));
%     for j=1:length(uni_frs)-1
%         cur_cadidates = curTrack(curTrack==uni_frs(j));
%         nei_cadidates = curTrack(curTrack==uni_frs(j+1));
%         for cc=1:length(cur_cadidates)
%             tmpDist = movieInfo.CDist{cur_cadidates(cc)};
%             for nn = 1:length(nei_cadidates)
%                 loc = find(movieInfo.nei{cur_cadidates(cc)}...
%                     ==nei_cadidates(nn));
%                 if ~isempty(loc)
%                     cd_cnt = cd_cnt + 1;
%                     cur_dist(cd_cnt) = tmpDist(loc);
%                 else
%                     keyboard;
%                 end
%             end
%         end
%     end
%   cur_dist = cur_dist(1:cd_cnt);
    devDist{i} = cur_dist;
end
all_dev = cat(1, devDist{:});
all_dev(isnan(all_dev(:,1)) | isinf(all_dev(:,1)),:) = [];
%phatDist = gamfit(all_dev);
phatDist = fitTruncGamma(all_dev);

end