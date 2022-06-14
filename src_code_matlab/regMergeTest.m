function mergeReg = regMergeTest(neighbors, cur_reg_id, movieInfo)
% test if the neighbors should be merged to current region
%INPUT:
% neighbors: ids of neighboring regions
% cur_reg_id: id of tested region
% movieInfo:
%OUTPUT:
% mergeReg: cells, each of which containing the ids should be merged
% 
% contact: ccwang@vt.edu, 03/04/2020

cur_voxIdx = movieInfo.voxIdx{cur_reg_id};
frs = movieInfo.frames(neighbors);
uni_fr = unique(frs);
% cur_tr_id = movieInfo.particle2track(cur_reg_id, 1);
% cur_track = movieInfo.tracks{cur_tr_id};
mergeReg = cell(length(uni_fr), 1);
for i=1:length(uni_fr)
    cur_neis = neighbors(frs == uni_fr(i));
    voxIdxes = movieInfo.voxIdx(cur_neis);
    if length(voxIdxes)>1
        for j=1:numel(voxIdxes)
            intersect_reg = intersect(voxIdxes{j}, cur_voxIdx);
            if length(intersect_reg) > length(voxIdxes{j})/2
                mergeReg{i} = cat(1, mergeReg{i}, cur_neis(j));
            end
        end
    end
end
mrl = cellfun(@length, mergeReg);
mergeReg(mrl<2) = []; % if only one region, no need to merge
end