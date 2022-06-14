function orgIm3d = drawTrack_klb(orgIm3d, tracks, detections, p,f)
% using detection id to draw tracks
% p: particleSize, lineWidth, cmap (color map);
particleSize = p.particleSize;

lineWidth = p.lineWidth;

clMap = p.cmap;
for i=1:numel(tracks)
    curTrack = tracks{i};
%     if isempty(find(detections(curTrack,4) == f,1))
%         continue;
%     end
    loc = find(detections(curTrack,end) <= f);
    curTrack = curTrack(1:loc(end));
    fprintf('%d track of %d tracks.\n', i, numel(tracks));
    particleCl = clMap(i,:);
    RegLineCl = clMap(i + 20,:);
    if length(curTrack)<1
        continue;
    end
    pt = detections(curTrack(1),1:3);
    orgIm3d = draw3Dparticle(orgIm3d,  pt, particleSize, particleCl);
    if length(curTrack)>1
        %pt = detections(curTrack(end),1:3);
        %orgIm3d = draw3Dparticle(orgIm3d,  pt, particleSize, particleCl);
        for j = 1:length(curTrack)-1
            tail = curTrack(j);
            head = curTrack(j+1);
            tPt = detections(tail,1:3);
            hPt = detections(head,1:3);
            orgIm3d = draw3Dline(orgIm3d, hPt, tPt, lineWidth, RegLineCl);
        end
    end

end