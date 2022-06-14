function tracks = movieInfo2tracks(movieInfo)
tracks = cell(numel(movieInfo.tracks), 1);
for i=1:numel(movieInfo.tracks)
    cur_track = movieInfo.tracks{i};
    tracks{i} = zeros(length(cur_track), 4);
    tracks{i}(:,1) = movieInfo.frames(cur_track);
    
    for j=1:length(cur_track)
        idx = movieInfo.voxIdx{cur_track(j)};
        [yy, xx, zz] = ind2sub(size(movieInfo.validGapMaps{1}), idx);
        tracks{i}(j,2:4) = mean([yy, xx, zz]);
    end
end
end