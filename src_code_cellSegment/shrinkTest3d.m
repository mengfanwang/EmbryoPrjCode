function [outLabel, reg_split] = shrinkTest3d(inLabel, shrinkScale, comMaps, q)
% if its real cell, not matter how we shrink, it should still be a unified
% region
% INPUT:
% newLabel: the input cell detection results
% regTest: 1=>pixels of edge between two cells, 2=>pixels of cell for testing,
% 3=>pixels of the other cell for testing
% shrinkScale: the scale to shrink the data
% OUTPUT:
% outLabel: the output cell detection results
% outReg: the output edge, region label map

% contact: ccwang@vt.edu, 02/07/2020

if length(shrinkScale) == 1
    shrinkScale = [shrinkScale,shrinkScale];
end
if (shrinkScale(2)) > 0
    se1 = strel('sphere', shrinkScale(1));
    se1 = se1.Neighborhood;
    se2 = strel('sphere', shrinkScale(2));
    se2 = se2.Neighborhood;
    if shrinkScale(1) > shrinkScale(2)
        delta = shrinkScale(1) - shrinkScale(2);
        mid = se1(:,:,shrinkScale(1) + 1);
        se1 = se1 * 0;
        se1(delta+1:end-delta, delta+1:end-delta, delta+1:end-delta)...
            = se2;
        se1(:,:,shrinkScale(1) + 1) = mid;
        se = se1;
    else
        delta = shrinkScale(2) - shrinkScale(1);
        mid = se2(:,:,shrinkScale(2) + 1);
        se2 = se2 * 0;
        se2(delta+1:end-delta, delta+1:end-delta, delta+1:end-delta)...
            = se1;
        se2(:,:,shrinkScale(2) + 1) = mid;
        se = se2;
    end
else
    se = strel('sphere', shrinkScale(1));
    se = se.Neighborhood;
    mid = se(:,:,shrinkScale(1) + 1);
    se = se * 0;
    se(:,:,shrinkScale(1) + 1) = mid;
end
%se = offsetstrel('ball', shrinkScale(1), shrinkScale(2));

n_old = max(inLabel(:));
reg_num = 0;
reg_split = false;
outLabel = zeros(size(inLabel));
for i=1:n_old
    loopCnt = 1;
    tmp_label = inLabel==i;
    while 1
        tmp_label = imerode(tmp_label, se);
        loopCnt = loopCnt + 1;
        if loopCnt >= 2
            break;
        end
    end
    tmp_label = bwareaopen(tmp_label, 10, 6);
    [tmp_id_label, n] = bwlabeln(tmp_label, 6);
    fgIn = inLabel==i;
    if n<2 % no new region appears
        n = 1;
        outLabel(fgIn) = n + reg_num;
    else % new region appears
        reg_split = true;
        %fprintf('split by shrinking!\n');
        [newLabel, n] = regionGrow(tmp_id_label, ...
            comMaps.score2dMap + comMaps.score3dMap,...
            fgIn, q.growConnectInRefine, ...
            q.cost_design, false);
        outLabel(fgIn) = newLabel(fgIn) + reg_num;
    end
    
    reg_num = reg_num + n;
end

end