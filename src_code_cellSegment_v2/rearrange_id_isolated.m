function [newIdMap, cnt] = rearrange_id_isolated(idMap, minSz)
% INPUT:
% idMap: the original id map, the ids may not be consecutive
% OUTPUT:
% newIdMap: new id map, the ids are consecutive, and each isolated region
% has a unique id. too small region will be removed
    com = bwconncomp(idMap>0);
    newIdMap = zeros(size(idMap), class(idMap));
    cnt = 0;
    for ii = 1:com.NumObjects
        if length(com.PixelIdxList{ii}) >= minSz
            cnt = cnt + 1;
            newIdMap(com.PixelIdxList{ii}) = cnt;
        end
    end

end