function [movieInfo, movieInfo_noJump] = handleMergeSplitRegions(movieInfo, g)
% add multi-to-one or one-to-multi choices; for >1 regions that has high 
% overlapping ratios with one pre- or post- neighbor, add links for such
% region to src or sink such that one region can incident with >1 regions


% we first build a movieInfo with no jump
movieInfo_noJump = movieInfo;
% noJumpAll: all linking has no jump; 
% noJump: no jump for spliting and merge
% none: consider jump
if strcmp(g.splitMergeHandle, 'noJumpAll')
    for i=1:numel(movieInfo_noJump.nei)
        % we only need to update the following 7 field
        if ~isempty(movieInfo_noJump.nei)
            validNeiIdx = (movieInfo_noJump.frames(movieInfo_noJump.nei{i}) == ...
                movieInfo_noJump.frames(i)+1) & ~isinf(movieInfo_noJump.Cij{i});
            movieInfo_noJump.nei{i} = movieInfo_noJump.nei{i}(validNeiIdx);
            movieInfo_noJump.Cij{i} = movieInfo_noJump.Cij{i}(validNeiIdx);
            movieInfo_noJump.ovSize{i} = movieInfo_noJump.ovSize{i}(validNeiIdx);
            movieInfo_noJump.CDist{i} = movieInfo_noJump.CDist{i}(validNeiIdx);
            movieInfo_noJump.CDist_i2j{i} = movieInfo_noJump.CDist_i2j{i}(validNeiIdx,:);
        end
        if ~isempty(movieInfo_noJump.preNei)
            validPreNeiIdx = movieInfo_noJump.frames(movieInfo_noJump.preNei{i}) == ...
                movieInfo_noJump.frames(i)-1 & ~isinf(movieInfo_noJump.Cji{i});
            movieInfo_noJump.preNei{i} = movieInfo_noJump.preNei{i}(validPreNeiIdx);
            movieInfo_noJump.Cji{i} = movieInfo_noJump.Cji{i}(validPreNeiIdx);
            movieInfo_noJump.preOvSize{i} = movieInfo_noJump.preOvSize{i}(validPreNeiIdx);
            movieInfo_noJump.CDist_j2i{i} = movieInfo_noJump.CDist_j2i{i}(validPreNeiIdx);
        end
    end
end

% get the arcs about split/merge among regions
[dat_in_append, movieInfo_noJump] = highOvEdges(movieInfo_noJump, g);
% get the other arcs
[~, g, dat_in] = trackGraphBuilder_cell(movieInfo_noJump, g);
if ~isempty(dat_in_append)
    fprintf('%d added edges\n',size(dat_in_append,1));
    dat_in = cat(1, dat_in, dat_in_append);
end
% linking with new graph design and update movieInfo (not movieInfo_noJump)
movieInfo = mccTracker(dat_in, movieInfo, g);

end