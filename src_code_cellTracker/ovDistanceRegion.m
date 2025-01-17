function [maxDistance, minDistance, distances, re_ratio] = ovDistanceRegion...
    (curRegVox, nextRegVox, frame_shift, ovFlag)
% cal the overlapping distance between two regions
%INPUT:
% curRegIdx: coordinates of the current region y,x,z (or x, y, z)
% nextRegIdx: coordinates of the neighboring region with order the same as curRegIdx
%!!! frame_shift: should also be the same order as curRegIdx
% ovFlag: whether using overlapping ratio as distance
%OUTPUT:
% distance: overlapping distances between this two region

% contact: ccwang@vt.edu, 02/26/2020
if nargin < 3
    frame_shift = [];
    ovFlag = false;
elseif nargin < 4
    ovFlag = false;
end

re_ratio = 0;
if ~ovFlag % after downsampling, we can use this method to cal distance
    if size(nextRegVox,1) < 2 || size(curRegVox,1) < 2 % less than 2 pixels
        distances = [100, 100];
    else
        if numel(frame_shift) ~= 3
            %% way 1: O(n^2)
            DY = distanceMat2Mex(curRegVox(:,1),nextRegVox(:,1));
            DX = distanceMat2Mex(curRegVox(:,2),nextRegVox(:,2));
            DZ = distanceMat2Mex(curRegVox(:,3),nextRegVox(:,3));
            distances = sqrt(DX.^2 + DY.^2 + DZ.^2);
            distances_c2n = min(distances, [], 2);
            distances_n2c = min(distances, [], 1);
            distances = [mean(distances_c2n), mean(distances_n2c)];
        else
            %% way 2: bwdist (O(n))
            %         st_pt = min([min(curRegVox); min(nextRegVox)]);
            %         end_pt = max([max(curRegVox); max(nextRegVox)]);
            %         bw_sz = ceil(end_pt-st_pt + 1);
            %         bw = false(bw_sz);
            %         cell1_sub = curRegVox - st_pt + 1;
            %         cell1_idx = sub2ind(bw_sz, cell1_sub(:,1), cell1_sub(:,2), cell1_sub(:,3));
            %         cell2_sub = nextRegVox - st_pt + 1;
            %         cell2_idx = sub2ind(bw_sz, cell2_sub(:,1), cell2_sub(:,2), cell2_sub(:,3));
            %         bw(cell1_idx) = true;
            %         dist2cell1 = bwdist(bw);
            %         distances_n2c_way2 = dist2cell1(cell2_idx);
            %
            %         bw(cell1_idx) = false;
            %         bw(cell2_idx) = true;
            %         dist2cell2 = bwdist(bw);
            %         distances_c2n_way2 = dist2cell2(cell1_idx);
            %         distances = [mean(distances_c2n_way2), mean(distances_n2c_way2)];
            %         redundant_sz = prod(bw_sz) - length(cell2_idx);
            %         re_ratio = redundant_sz/prod(bw_sz);
            
            %% way 3: lower envelops (O(n))
            %         st_pt = min([min(curRegVox); min(nextRegVox)]);
            %         end_pt = max([max(curRegVox); max(nextRegVox)]);
            %         bw_sz = ceil(end_pt-st_pt + 1);
            %         bw = false(bw_sz);
            %         cell1_sub = curRegVox - st_pt + 1;
            %         cell1_idx = sub2ind(bw_sz, cell1_sub(:,1), cell1_sub(:,2), cell1_sub(:,3));
            %         cell2_sub = nextRegVox - st_pt + 1;
            %         cell2_idx = sub2ind(bw_sz, cell2_sub(:,1), cell2_sub(:,2), cell2_sub(:,3));
            %         bw(cell1_idx) = true;
            %
            %         dist2cell1 = edt_3d(bw);
            %         distances_n2c_way3 = sqrt(dist2cell1(cell2_idx));
            %
            %         bw(cell1_idx) = false;
            %         bw(cell2_idx) = true;
            %         dist2cell2 = edt_3d(bw);
            %         distances_c2n_way3 = sqrt(dist2cell2(cell1_idx));
            
            %% way 4: lower envelop with cropped region
            %         st_pt = min(curRegVox);
            %         end_pt = max(curRegVox);
            %         bw_sz = ceil(end_pt-st_pt + 1); % must be integer
            %         ref_cell = false(bw_sz);
            %         cell1_sub = curRegVox - st_pt + 1;
            %         cell1_idx = sub2ind(bw_sz, cell1_sub(:,1), cell1_sub(:,2), cell1_sub(:,3));
            %         ref_cell(cell1_idx) = true;
            %         boxSize = prod(bw_sz);
            %
            %         st_pt2 = min(nextRegVox);
            %         end_pt2 = max(nextRegVox);
            %         bw_sz = ceil(end_pt2-st_pt2 + 1);
            %         mov_cell = false(bw_sz);
            %         cell2_sub = nextRegVox - st_pt2 + 1;
            %         cell2_idx = sub2ind(bw_sz, cell2_sub(:,1), cell2_sub(:,2), cell2_sub(:,3));
            %         mov_cell(cell2_idx) = true;
            %         boxSize = boxSize + prod(bw_sz);
            %
            %         mov_shift = st_pt2 - st_pt - frame_shift;
            %         dist2cell1 = edt_3d_regFixed(ref_cell, mov_cell, mov_shift);
            %         distances_n2c_way4 = sqrt(dist2cell1(cell2_idx));
            %
            %         dist2cell2 = edt_3d_regFixed(mov_cell, ref_cell, -mov_shift);
            %         distances_c2n_way4 = sqrt(dist2cell2(cell1_idx));
            %% way 5 fully c++ based way4
            st_pt = min(curRegVox, [], 1);
            end_pt = max(curRegVox, [], 1);
            bw_sz = ceil(end_pt-st_pt + 1); % must be integer
            ref_cell = false(bw_sz);
            cell1_sub = curRegVox - st_pt + 1;
            cell1_idx = sub2ind(bw_sz, cell1_sub(:,1), cell1_sub(:,2), cell1_sub(:,3));
            ref_cell(cell1_idx) = true;
            boxSize = prod(bw_sz);
            
            st_pt2 = min(nextRegVox);
            end_pt2 = max(nextRegVox);
            bw_sz = ceil(end_pt2-st_pt2 + 1);
            mov_cell = false(bw_sz);
            cell2_sub = nextRegVox - st_pt2 + 1;
            cell2_idx = sub2ind(bw_sz, cell2_sub(:,1), cell2_sub(:,2), cell2_sub(:,3));
            mov_cell(cell2_idx) = true;
            boxSize = boxSize + prod(bw_sz);
            
            mov_shift = single(st_pt2 - st_pt - frame_shift);
            dist2cell1 = edt_3dMex(ref_cell, mov_cell, mov_shift);
            distances_n2c_way5 = dist2cell1(cell2_idx);%sqrt(dist2cell1(cell2_idx));
            dist2cell2 = edt_3dMex(mov_cell, ref_cell, -mov_shift);
            distances_c2n_way5 = dist2cell2(cell1_idx);%sqrt(dist2cell2(cell1_idx));
            if false%q.sqrtDistance
                distances_n2c_way5 = sqrt(distances_n2c_way5);
                distances_c2n_way5 = sqrt(distances_c2n_way5);
            end
            distances = [mean(distances_c2n_way5), mean(distances_n2c_way5)]; % i2j and j2i
            
            redundant_sz = boxSize - length(cell2_idx);
            re_ratio = redundant_sz/boxSize;
            
            %         if sum(abs([distances_c2n(:)-distances_c2n_way5(:)]))>0.01 || ...
            %                 sum(abs([distances_n2c(:)-distances_n2c_way5(:)]))>0.01
            %             disp(sum(abs([distances_c2n-distances_c2n_way5])));
            %         end
        end
    end
else % purely based on overlapping ratio
    C = intersect(curRegVox, nextRegVox,'rows');
    distances = size(C,1) / (size(curRegVox,1) + size(nextRegVox,1) - size(C,1));
    distances = -log(distances);
end
maxDistance = max(distances);
minDistance = min(distances);
end