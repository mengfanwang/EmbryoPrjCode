function [refine_res, merge_cnt] = track2seg(movieInfo, refine_res, vidMap, simple_merge)
% refine the segmentation results using tracking results
%INPUT:
% movieInfo: all information about segmentation and tracks
% refine_res: cells containing the label maps
% vidMap: the original images(cells each of which corresponds to one frame)
% simple_merge: simply label the regions for merging as the same index
%OUTPUT:

% contact: ccwang@vt.edu, 03/04/2020

% find previous neighbors
% preNei = cell(length(movieInfo.xCoord),1);
% for i=1:numel(movieInfo.nei)
%     for j=1:length(movieInfo.nei{i})
%         preNei{movieInfo.nei{i}(j)} = cat(1, preNei{movieInfo.nei{i}(j)}, i);
%     end
% end
% movieInfo.preNei = preNei;


track_len = cellfun(@length, movieInfo.tracks);
[~, od] = sort(track_len, 'descend');
merge_cnt = 0;
for tr_id=1:numel(movieInfo.tracks)
    i = od(tr_id);
    fprintf('processing %d/%d tracks, id:%d\n', tr_id, numel(movieInfo.tracks),i);
    
    % way 1: purely based on overlapping ratio: merge segments highly
    % overlapped with the same segment in the previous/following frame
    %     voxIdxCells = movieInfo.voxIdx(movieInfo.tracks{i});
    %     sz_trend = cellfun(@length, voxIdxCells);
    %     while true
    %         [~, id] = max(sz_trend);
    %         cur_reg_id = movieInfo.tracks{i}(id);
    %         preNei = movieInfo.preNei{cur_reg_id};
    %         nei = movieInfo.nei{cur_reg_id};
    %         neighbors = cat(1, preNei, nei);
    %         mergeReg = regMergeTest(neighbors, cur_reg_id, movieInfo);
    %
    %         [refine_res, movieInfo] = mergedRegGrow(mergeReg, vidMap, ...
    %             refine_res, movieInfo, simple_merge);
    %         %movieInfo.voxIdx{cur_reg_id} = [];
    %         movieInfo.tracks{i}(id) = [];
    %         voxIdxCells(id) = [];
    %         if isempty(voxIdxCells)
    %             break;
    %         end
    %         sz_trend = cellfun(@length, voxIdxCells);
    %     end
    % way 2: use the tracking results: if two segments in the same frame
    % belong to the same track, they should be merged
    cur_track = movieInfo.tracks{i};
    if length(cur_track)<2
        continue;
    end
    cur_frames = movieInfo.frames(cur_track);
    locs = find((cur_frames(1:end-1) - cur_frames(2:end))==0);
    merge_check = unique(cur_frames(locs));
    mergeReg = cell(length(merge_check), 1);
    for j = 1:length(merge_check)
        mergeReg{j} = cur_track(cur_frames==merge_check(j));
    end
    [refine_res, movieInfo] = mergedRegGrow(mergeReg, vidMap, ...
        refine_res, movieInfo, simple_merge);
    merge_cnt = merge_cnt + numel(mergeReg);
end

for i=1:numel(refine_res)
    refine_res{i} = rearrange_id(refine_res{i});
%     [~, flag] = region_sanity_check(refine_res{i});
%     if flag 
%         keyboard;
%         zzshow(label2rgb3d(refine_res{i}));
%     end
end
end