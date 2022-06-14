function [dat_in, src_node, sink_node] = graphCut_4d(scoreMap, fMapOrg, ...
    sMap, tMap, punish, connect, cost_design)
% build a graph for graph cut with temporal information
% contact: ccwang@vt.edu

% sMap = newLabel{2} == cur_id;
% if ~isempty(fMapOrg{1})
%     sMap = sMap & fMapOrg{1};
% end
% if ~isempty(fMapOrg{3})
%     sMap = sMap & fMapOrg{3};
% end    
% if isempty(find(sMap,1))
%     dat_in = [];
%     return;
% end
vox_num = numel(scoreMap{2});
fMapAll = fMapOrg{2};
if ~isempty(fMapOrg{1})
    fMapAll = fMapAll | fMapOrg{1};
end
if ~isempty(fMapOrg{3})
    fMapAll = fMapAll | fMapOrg{3};
end
fMap = cell(3,1);
fMap{2} = fMapAll;
if ~isempty(fMapOrg{1})
    fMap{1} = fMapAll;
end
if ~isempty(fMapOrg{3})
    fMap{3} = fMapAll;
end
% 
valid_vox = cell(3,1);
for i=1:3
    valid_vox{i} = find(fMap{i});
end
vox_len = cellfun(@length, valid_vox);
for i=2:numel(fMap)
    valid_vox{i} = valid_vox{i} + vox_num * (i-1);
end
valid_vox = cat(1, valid_vox{:});
valid_vox_num = length(valid_vox);
adjMap = cell(3,1);
for i=1:3
    adjMap{i} = zeros(size(scoreMap{2})); % map of neighbors
end

dat_in = zeros(valid_vox_num*connect*2 + vox_num*7, 3);
k_dat = 0;
% connections inside principal map or gradient map
for i=1:valid_vox_num
    cur_vox_id = valid_vox(i);
    
    if i > vox_len(1) + vox_len(2)
        increment = vox_num*2;
        cur_fr = 3;
    elseif i > vox_len(1)
        increment = vox_num;
        cur_fr = 2;
    else
        increment = 0;
        cur_fr = 1;
    end
    ind = cur_vox_id - increment;
    vid = scoreMap{cur_fr};
    p1 = vid(ind);
    [~, ~, nei_ids] = neighbours(ind, vid, connect);
    p2 = vid(nei_ids);
    if cost_design(1)==1
        costs = (2./(p1+p2)).^cost_design(2);
    elseif cost_design(1)==2
        p2(p2==0) = p1; % or inf?
        costs = (1./sqrt(p1.*p2)).^cost_design(2);
    end

    k_dat = k_dat + length(nei_ids);
    dat_in(k_dat-length(nei_ids)+1:k_dat,:) = ...
        [cur_vox_id + nei_ids*0, increment+nei_ids, costs];

    k_dat = k_dat + length(nei_ids);
    dat_in(k_dat-length(nei_ids)+1:k_dat,:) = ...
        [increment+nei_ids, cur_vox_id + nei_ids*0, costs];
    if cur_fr == 2 % add temporal links
        if ~isempty(fMap{1})
            if fMap{1}(ind) > 0
                left_id = ind;
                k_dat = k_dat + 1;
                dat_in(k_dat,:) = [left_id, cur_vox_id, punish];
                k_dat = k_dat + 1;
                dat_in(k_dat,:) = [cur_vox_id, left_id, punish];
            end
        end
        if ~isempty(fMap{3})
            if fMap{3}(ind) > 0
                right_id = ind + 2*vox_num;
                k_dat = k_dat + 1;
                dat_in(k_dat,:) = [right_id, cur_vox_id, punish];
                k_dat = k_dat + 1;
                dat_in(k_dat,:) = [cur_vox_id, right_id, punish];
            end
        end
    end
    %adjMap(nei_ids) = 1;
end
% connections to src or sink
% srcIds = cell(3,1);
% srcIds{2} = find(sMap) + vox_num;
% for i=[1 3]
%     if ~isempty(newLabel{i})
%         srcIds{i} = srcIds{2} + vox_num * (i-2);
%     end
% end
src_node = vox_num*3+1;
if ~iscell(sMap)
    srcIds = find(sMap) + vox_num;
    k_dat = k_dat + length(srcIds);
    dat_in(k_dat-length(srcIds)+1:k_dat,:) = [src_node + srcIds*0, srcIds, inf + srcIds*0];
else
    for i=1:numel(sMap)
        if ~isempty(sMap{i})
            srcIds = find(sMap{i}) + vox_num*(i-1);
            k_dat = k_dat + length(srcIds);
            dat_in(k_dat-length(srcIds)+1:k_dat,:) = [src_node + srcIds*0, srcIds, inf + srcIds*0];
        end
    end
end

sink_node = vox_num*3+2;
if ~iscell(tMap)
    sinkIds = find(tMap) + vox_num;
    k_dat = k_dat + length(sinkIds);
    dat_in(k_dat-length(sinkIds)+1:k_dat,:) = [sinkIds, sink_node + sinkIds*0, inf + sinkIds*0];
else
    for i=1:numel(tMap)
        if ~isempty(tMap{i})
            sinkIds = find(tMap{i}) + vox_num*(i-1);
            k_dat = k_dat + length(sinkIds);
            dat_in(k_dat-length(sinkIds)+1:k_dat,:) = [sinkIds, sink_node + sinkIds*0, inf + sinkIds*0];
        end
    end
end
dat_in = dat_in(1:k_dat,:);
end