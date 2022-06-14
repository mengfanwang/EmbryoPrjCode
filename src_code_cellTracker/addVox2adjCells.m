function [movieInfo, refine_res, thresholdMaps, voxIdxCells] = ...
    addVox2adjCells(comMaps, idsAndLocs, movieInfo, refine_res, thresholdMaps)
% add voxes to an existing cell region
voxIdxCells = cell(size(idsAndLocs, 1), 2);
if isempty(idsAndLocs)
    return;
end
for i=1:size(idsAndLocs, 1)
    idx_all = comMaps.linerInd(comMaps.idComp==idsAndLocs(i,1));
    fr = movieInfo.frames(idsAndLocs(i,2));
    idx = idx_all(refine_res{fr}(idx_all) == 0);
    if ~isempty(idx)
        thresholdMaps{fr}(idx) = comMaps.pickedThreshold; % there shold be only one
        refine_res{fr}(idx) = idsAndLocs(i,1);
        voxIdxCells{i,1} = movieInfo.voxIdx{idsAndLocs(i,2)};
        voxIdxCells{i,2} = idx;
        movieInfo.voxIdx{idsAndLocs(i,2)} = cat(1, ...
            movieInfo.voxIdx{idsAndLocs(i,2)}, idx);
    end
end




end