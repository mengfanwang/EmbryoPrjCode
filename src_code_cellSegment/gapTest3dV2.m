function [newLabel, regTest] = gapTest3dV2(L, comMaps, reg_id, q, OrSt)
% test if the gap is significant enough in a 3D manner
other_id = comMaps.idComp~=reg_id & comMaps.idComp>0 & comMaps.fmapComp;
append_id = nan;
if ~isempty(find(other_id, 1))
    append_id = max(L(:)) + 1;
    L(other_id) = append_id;
end                                                                         % other regions don't change
% key1: we should link boundary to background, otherwise the boundary has
% not penalty to be cut and the max-flow will prefer to grow to boundary if
% foreground cover them.
fmapIn = comMaps.fmapComp;
% [h, w, z] = size(fmapIn);
% tmpMap = false(h,w);
% for i=1:1 % the one cycles of boundary are background
%     tmpMap(i,:) = true;
%     tmpMap(end-i+1,:) = true;
%     tmpMap(:,i) = true;
%     tmpMap(:, end-i+1) = true;
% end
% bndIdx = find(tmpMap);
% for i=1:z
%     cur_bidx = bndIdx + (i-1)*w*h;
%     cur_bidx = cur_bidx(L(cur_bidx)==0 & comMaps.idComp(cur_bidx)==0); % do not affect seeds
%     fmapIn(cur_bidx) = 0;
% end
% key2: here we must consider background as sink, otherwise all fg pixels will
% be assigned to one label; Now we only allow 2d grow, so this will cause
% some problem if e.g. a slice only contains one seed while indeed has two
% cells
bg2sink = false;
newLabel = regionGrow(L, comMaps.score2dMap + comMaps.score3dMap, ...
    fmapIn, q.growConnectInRefine, q.cost_design, bg2sink);
% should we first 2d, then 3d??? Practical not good, but why???
% bg2sink = true;
% newLabel = regionGrow(L, ...
%     comMaps.score2dMap,...
%     fmapIn, q.growConnectInTest, ...
%     q.cost_design, bg2sink);
% bg2sink = false;
% newLabel = regionGrow(newLabel, ...
%     comMaps.score2dMap+comMaps.score3dMap,...
%     fmapIn, q.growConnectInRefine, ...
%     q.cost_design, bg2sink);

if ~isnan(append_id)
    newLabel(newLabel == append_id) = 0;
end

[newLabel, regTest] = edgeTest3dV2(comMaps.vidStbComp, newLabel, comMaps.fmapComp, q.edgeConnect, OrSt.p_thres, OrSt);



end