function [res_G, pSet, costs] = successive_shortest_paths(orgG, s, t)
% successive shortest path algorithm
vNum = orgG.numnodes;
res_G = orgG;
pSet = cell(vNum,1);
costs = zeros(vNum, 1);
numTracks = 0;
tic;
while true
%     if i==940
%         keyboard;
%     end
    [P,d] = shortestpath(res_G,s,t); % 1 is source, n_nodes is sink saved in g.excess_node
    if d>=-0
        break;
    end
%     if sum(P==14082*2+1)>0
%         %keyboard;
%         pp = P(mod(P,2)==0);
%         pp = pp/2;
%     end
    %fprintf('The %d paths found with cost: %.2f\n', numTracks, d);
    % change weight as -1*weight
%     for j=1:length(P)-1
%         tmpIdx = find(res_G.Edges.EndNodes(:,1) == P(j) & res_G.Edges.EndNodes(:,2) == P(j+1)); % should only have one
%         if length(tmpIdx) ~= 1
%             error('There should be only one edge be changed.');
%         end
%         res_G.Edges.Weight(tmpIdx) = -1*res_G.Edges.Weight(tmpIdx);
%     end
    edgeIdx = findedge(res_G, P(1:end-1),P(2:end));
    res_G.Edges.Weight(edgeIdx) = -1*res_G.Edges.Weight(edgeIdx);
    % inverse edge
    res_G = flipedge(res_G,P(1:end-1),P(2:end));
    numTracks = numTracks + 1;
    pSet{numTracks} = P;
    costs(numTracks) = d;
end
timeLapse = toc;
fprintf('finish successive shortest path alg with %d paths, %.3fs!\n',numTracks, timeLapse);
pSet = pSet(1:numTracks);
costs = costs(1:numTracks);