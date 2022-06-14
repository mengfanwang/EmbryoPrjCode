function orgIm3d = drawTrack_klb_loc(orgIm3d, tracks, p,f)
% using lcoations id to draw tracks
% p: particleSize, lineWidth, cmap (color map);
particleSize = p.particleSize;

lineWidth = p.lineWidth;

clMap = p.cmap;
for i=1:numel(tracks)
    curTrack = tracks{i};
    loc = find(curTrack(:,1) <= f);
    if isempty(loc)
        continue;
    end
    curTrack = curTrack(1:loc(end),:);
    fprintf('%d track of %d tracks.\n', i, numel(tracks));
    particleCl = clMap(i,:);
    RegLineCl = clMap(i + 20,:);
    if size(curTrack,1)<1
        continue;
    end
    pt = curTrack(1,2:4);
    orgIm3d = draw3Dparticle(orgIm3d,  pt, particleSize, particleCl);
    if size(curTrack,1)>1
        for j = 1:size(curTrack,1)-1
            tPt = curTrack(j,2:4);
            hPt = curTrack(j+1,2:4);
            orgIm3d = draw3Dline(orgIm3d, hPt, tPt, lineWidth, RegLineCl);
        end
    end

end