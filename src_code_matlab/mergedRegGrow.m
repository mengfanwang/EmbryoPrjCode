function [refine_res, movieInfo, newlyReg] = mergedRegGrow(...
    mergeReg, vidMap, refine_res, movieInfo, simple_merge, merge_error_prob,...
    stableNodeTest)
% region grow to merge several given seed regions
%INPUT:
% mergeReg: cells, each cell containing the regions that should be merged
% scoreMap: cells, each containing the score map for region grow
% refine_res: cells, label map of each frame
% movieInfo: information about cells and tracks
% simple_merge: simply label the regions for merging as the same index
%OUTPUT:

% ccwang@vt.edu, 03/04/2020
if nargin == 5
    merge_error_prob = 0.5 + zeros(numel(mergeReg),1);
    stableNodeTest = false;
end

connect = 4;
cost_design = [1 2];
smooth_sc = 1;
shift = [3 3 2];
%newlyReg = cell(numel(mergeReg), 1);
for i=1:numel(mergeReg)
    cur_frame = unique(movieInfo.frames(mergeReg{i}));
    if length(cur_frame) > 1
        %keyboard;
        mergeReg{i} = [];
        continue;
    end
    if stableNodeTest && ...
            ~isempty(find(~isinf(movieInfo.arc_avg_mid_std(mergeReg{i},4)),1))
        continue;
    end
    real_ids = nan(length(mergeReg{i}), 1);
    for j=1:length(real_ids)
        if ~isempty(movieInfo.voxIdx{mergeReg{i}(j)})
            real_ids(j) = refine_res{cur_frame}(movieInfo.voxIdx{mergeReg{i}(j)}(1));
            if real_ids(j) + sum(movieInfo.n_perframe(1:cur_frame-1))...
                    ~= mergeReg{i}(j)
                keyboard;
            end
        end
    end
    % avoid duplicated merge
    [real_ids, ia] = unique(real_ids);
    if ~isempty(find(real_ids==0,1))
        keyboard;
    end
    if length(real_ids)<2
        continue;
    else
        mergeReg{i} = mergeReg{i}(ia);
    end
    yxz = cat(1, movieInfo.voxIdx{mergeReg{i}});
    if simple_merge
        new_idx = yxz;
    else
        [vidComp, linerInd] = crop3D(vidMap{cur_frame}, yxz, shift);
        idComp = crop3D(refine_res{cur_frame}, yxz, shift);
        eig2dComp = principalCv2d(vidComp, idComp, smooth_sc);%
        eig2dComp(eig2dComp<0) = 0;
        scoreMap = scale_image(eig2dComp, 1e-3,1);
        
        fMap = imfill(eig2dComp~=0 ,'holes');
        fMap(idComp>0) = 0;
        idMap = zeros(size(idComp));
        for j=1:length(real_ids)
            fMap(idComp==real_ids(j)) = 1;
            idMap = idMap | idComp==real_ids(j);
        end
        newLabel = regionGrow(idMap, scoreMap, fMap, connect, cost_design);
        new_idx = linerInd(newLabel>0);
    end
    % merge will no longer use their gaps since we merge the two regions
    if ~isempty(movieInfo.validGapMaps) && merge_error_prob(i) < 0.05
        % if we have >99% confidence that the merge is correct, we fix it,
        % it will never be split in the future
        movieInfo.validGapMaps{cur_frame}(new_idx) = false;
    end
    refine_res{cur_frame}(new_idx) = real_ids(1);
    % relabel the related regions in movieInfo (which are saved in refine_res)
    % with the same indexes, thus we can avoid duplicated merge
    movieInfo.voxIdx{mergeReg{i}(1)} = new_idx;
%     [y,x,z] = ind2sub(size(refine_res{cur_frame}), new_idx);
%     movieInfo.vox{mergeReg{i}(1)} = [x, y, z];
    for j=2:length(mergeReg{i})
        movieInfo.voxIdx{mergeReg{i}(j)} = [];
%         movieInfo.vox{mergeReg{i}(j)} = [];
    end
end
newlyReg = unique(cat(2, mergeReg{:}));

end