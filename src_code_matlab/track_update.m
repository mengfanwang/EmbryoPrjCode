function movieInfo = track_update(movieInfo, refine_res)
% after merging regions, the tracks also no longer valid, we will re-label
% the tracks based on new label maps

numTracks = numel(movieInfo.tracks);
region_cnt = cellfun(@(x) max(x(:)), refine_res);
new_tracks = cell(numTracks, 1);
for i=1:numTracks
    cur_track = movieInfo.tracks{i};
    cur_voxIdx = movieInfo.voxIdx(cur_track);
    region1stEles = cellfun(@(x) x(1), cur_voxIdx);
    cur_frames = movieInfo.frames(cur_track);
    newIds = zeros(length(cur_track),1);
    for j=1:length(cur_frames)
        newIds(j) = refine_res{cur_frames(j)}(region1stEles(j)) + ...
            sum(region_cnt(1:cur_frames(j)-1));
    end
    new_tracks{i} = unique(newIds);
%     cur_frames = movieInfo.frames(new_tracks{i});
%     uni_frs = unique(cur_frames);
%     for j=1:length(uni_frs)
%         nodes_cur_fr = new_tracks{i}(cur_frames==uni_frs(j));
%         if length(nodes_cur_fr)>1
%             movieInfo.vox{nodes_cur_fr(1)} = cat(1, movieInfo.vox{nodes_cur_fr});
%             movieInfo.voxIdx{nodes_cur_fr(1)} = cat(1, movieInfo.voxIdx{nodes_cur_fr});
%             new_tracks{i}(cur_frames==uni_frs(j)) = nodes_cur_fr(1);
%         end
%     end
%     new_tracks{i} = unique(new_tracks{i});
end

movieInfo = tree2tracks_cell(refine_res, false);

movieInfo.tracks = new_tracks;

particle2track = nan(sum(region_cnt), 2); % the track it belongs to and position
for i=1:numTracks
    particle2track(new_tracks{i},:) = [i+zeros(length(new_tracks{i}),1),...
        [1:length(new_tracks{i})]'];
end
movieInfo.particle2track = particle2track;

end