function [newLabel, regTest] = gapTest(L, n, vidComp, idComp, eigComp,...
    eigCompMap, connect, cost_design, edgeConnect, p_thres, OrSt)
%

%% 
newLabel = L;
for j=1:n
    sMap = L==j;
    if isempty(find(sMap,1))
        continue;
    end
    tMap = idComp<1 | (L~=j & L~=0);
    [dat_in, src, sink] = graphCut_build(eigComp,eigCompMap, sMap, tMap, connect, cost_design);
    G = digraph(dat_in(:,1),dat_in(:,2),dat_in(:,3));
    [~,~,cs,~] = maxflow(G, src, sink); % cs: fg, ct: bg
    cur_l = newLabel(cs(cs<numel(newLabel)));
    if length(unique(cur_l))>2
        keyboard;
    end
    newLabel(cs(cs<numel(newLabel))) = j;
end
% merge over-segmented regions
inLabel = newLabel;

[newLabel, regTest] = edgeTest(vidComp, inLabel, eigCompMap, edgeConnect, p_thres, OrSt);

end