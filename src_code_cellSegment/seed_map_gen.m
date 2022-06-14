function seedsMap = seed_map_gen(fg, gap3d, gap2d, min_seed_sz)
% 3d map might be too liberal to remove some valid seeds.
% we will check if there is a seed totally masked by 3d gap map from 3d
% principal curvature
if nargin == 3
    min_seed_sz = 0;
end
seedsMap = fg;
seedsMap(gap3d) = 0; % remove 3d principal curvature
if min_seed_sz > 0
    seedsMap = bwareaopen(seedsMap,min_seed_sz,26);
end
if ~isempty(find(seedsMap, 1))
    % if 3d map has removed all the seeds, just keep it. This region either
    % is not a cell or too samll to be valid
    seedsMap_foundBy2d = fg;
    seedsMap_foundBy2d(gap2d) = 0;
    [seedsMap_foundBy2d, n] = bwlabeln(seedsMap_foundBy2d | seedsMap, 26);

    exist_seeds = unique(seedsMap_foundBy2d(seedsMap));
    exist_seeds = exist_seeds(exist_seeds>0);
    if n > length(exist_seeds)                                              % find seed not discover by 3d map but by 2d map
        seedsMap_foundBy2d(ismember(seedsMap_foundBy2d, exist_seeds)) = 0;
        seedsMap = seedsMap | seedsMap_foundBy2d > 0;
        %display_masks(comMaps.vidComp, seedsMap, seedsMap_foundBy2d);
        %fprintf('we found possible lost seeds.\n');
    end
    if min_seed_sz > 0
        seedsMap = bwareaopen(seedsMap,min_seed_sz,26);
    end
end



