function [seed_region, missed_frames] = extractSeedRegFromGivenCell(...
    parent_kid_vec, movieInfo, thresholdMaps, embryo_vid)
% given two footprints of the same cell (or one), get the seed in the 
% jumped frame in between
% NOTE: here we do not consider if the seed is already occupied by an
% existing cell, we will check this in the following steps
p_id = parent_kid_vec(1);
k_id = parent_kid_vec(2);
if isnan(p_id) || isnan(k_id)% starting point
    if ~isnan(p_id)
        missed_frames = movieInfo.frames(p_id)+1;
        idx = movieInfo.voxIdx{p_id};
        thres = thresholdMaps{movieInfo.frames(p_id)}(idx);
    elseif ~isnan(k_id)
        missed_frames = movieInfo.frames(k_id)-1;
        idx = movieInfo.voxIdx{k_id};
        thres = thresholdMaps{movieInfo.frames(k_id)}(idx);
    else
        error('at least one element in parent_kid_vec should be valid!');
    end
    if find(thres==0,1)
        keyboard;
    end
    if isempty(thres(thres>0)) %|| thres == 0
        seed_region = [];
        return;
    else
        thres = mode(thres(thres>0));
        vals = embryo_vid{missed_frames}(idx);
        seed_region = idx(vals>thres);
        if isempty(seed_region)
            return;
        end
    end
else
    missed_frames = movieInfo.frames(p_id)+1 : ...
        movieInfo.frames(k_id)-1;
    % find seed region: can be more complicated
    seed_region = intersect(movieInfo.voxIdx{p_id}, movieInfo.voxIdx{k_id});
    if isempty(seed_region)
        return;
    end
end

%% pick the largest region indicated by the seed region
% way 1
%     tic;
%     bw = false(size(thresholdMaps{missed_frames}));
%     bw(seed_region) = true;
%     out_bw = pickLargestReg(bw, 6);
%     seed_region1 = find(out_bw);
%     toc;
% way 2
[bw, st_pt, idx] = mapIdx2Bw(seed_region, size(thresholdMaps{1}));
[bw, rm_flag, idx] = pickLargestReg(bw, 6, idx);
if rm_flag % otherwise, no need to change seed_region
    [yy, xx, zz] = ind2sub_direct(size(bw), idx);
    seed_region = sub2ind_direct(size(thresholdMaps{1}), yy+(st_pt(1)-1), ...
        xx+(st_pt(2)-1), zz+(st_pt(3)-1));
end
end