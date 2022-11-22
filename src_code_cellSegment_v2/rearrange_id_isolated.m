function [newIdMap, cnt] = rearrange_id_isolated(idMap, minSz)
% INPUT:
% idMap: the original id map, the ids may not be consecutive
% OUTPUT:
% newIdMap: new id map, the ids are consecutive, and each isolated region
% has a unique id. too small region will be removed

% bug: should consider different id
%     com = bwconncomp(idMap>0);
%     newIdMap = zeros(size(idMap), class(idMap));
%     cnt = 0;
%     for ii = 1:com.NumObjects
%         if length(com.PixelIdxList{ii}) >= minSz
%             cnt = cnt + 1;
%             newIdMap(com.PixelIdxList{ii}) = cnt;
%         end
%     end

    newIdMap = zeros(size(idMap));
    s = regionprops3(idMap, {'VoxelIdxList'});
    cnt = 0; 

    for ii=1:numel(s.VoxelIdxList)  
        fprintf("%d/%d\n", ii, numel(s.VoxelIdxList))
        if ~isempty(s.VoxelIdxList{ii})
            idMap_temp = zeros(size(idMap));
            idMap_temp(s.VoxelIdxList{ii}) = 1;
            com = bwconncomp(idMap_temp);
            for jj = 1:com.NumObjects
                if length(com.PixelIdxList{jj}) >= minSz
                    cnt = cnt + 1;
                    newIdMap(com.PixelIdxList{jj}) = cnt;
                end
            end
        end
    end

end