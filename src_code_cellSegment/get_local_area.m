function [comMaps, OrSt] = get_local_area(vid, vid_stb, idMap, i,...
    eig2d, eig3d, yxz, OrSt, q)
% crop out the local surrounding area for a given seed region
% i is the id labeling the seed region in idMap, yxz is its location
% indexs.
if iscell(vid)
    comMaps = cell(3,1);
    comMaps{1} = cropNeedRegions(vid{1}, vid_stb{1}, idMap{1}, nan,...
        [], [], yxz, OrSt.imProcMethod, q.shift);
    comMaps{2} = cropNeedRegions(vid{2}, vid_stb{2}, idMap{2}, i, ...
        eig2d, eig3d, yxz, OrSt.imProcMethod, q.shift);
    comMaps{3} = cropNeedRegions(vid{3}, vid_stb{3}, idMap{3}, nan, ...
        [], [], yxz, OrSt.imProcMethod, q.shift);
else
    comMaps = cropNeedRegions(vid, vid_stb, idMap, i, ...
        eig2d, eig3d, yxz, OrSt.imProcMethod, q.shift);
end
OrSt.stbVarCropMap = crop3D(OrSt.stbVarMap, yxz, q.shift);
OrSt.NoStbVarCropMap = crop3D(OrSt.NoStbVarMap, yxz, q.shift);