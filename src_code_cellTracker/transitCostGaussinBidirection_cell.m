function movieInfo = transitCostGaussinBidirection_cell(movieInfo, g)
% update transition cost using existing trajectories using both precursors
% and ancestors
%
phatDist = getStdFromTracks_cell(movieInfo, g);
movieInfo.ovGamma = phatDist;
%% we update cost of all edges rather than part of them
nei = movieInfo.nei;
newCij = cell(length(movieInfo.xCoord),1);
for i=1:length(movieInfo.xCoord)
    if mod(i, 100)==0
        disp(i);
    end
    curNode = i;
    neiUp = nei{i};
    newCijSingle = nan(length(neiUp),1);
    for nn=1:length(neiUp)
        % overlapping_cost
%         devXYZ = [0 0 0];
        ov_dist = movieInfo.CDist{i}(nn);
%         ov_dist = ovDistanceRegion(movieInfo.vox{curNode}+round(devXYZ),...
%             movieInfo.vox{neiUp(nn)});
%         p_oc = 1-gamcdf(ov_dist, phatDist(1), phatDist(2));
%         oc = norminv(1-p_oc/2);
% 
%         edgeCost = oc.^2;%
        edgeCost = overlap2cost(ov_dist, phatDist);
        newCijSingle(nn) = edgeCost;
        if isnan(edgeCost)
            fprintf('We found NaN edge cost!!!\n');
        end
    end
    newCij{i} = newCijSingle;
end
movieInfo.Cij = newCij;

Cji = cell(numel(nei),1);
for i=1:numel(nei)
    for j=1:length(nei{i})
        Cji{nei{i}(j)} = cat(1, Cji{nei{i}(j)}, newCij{i}(j));
    end
end

movieInfo.Cji = Cji;

fprintf('finish cost update!\n');

end