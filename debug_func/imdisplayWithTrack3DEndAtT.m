function orgIm3d = imdisplayWithTrack3DEndAtT(orgIm3d, in_tracks, timePt, p)
% generate a 3D colorful data with ROI overlay
% gt_tracks: t, y, x, z
particleSize = p.particleSize;
%particleCl = p.particleCl;
stoppedParticleCl = p.stoppedParticleCl;
lineWidth = p.lineWidth;
clMap = p.cmap;
[h,w,~,slice] = size(orgIm3d);
if ~iscell(in_tracks)
    gt_tracks = cell(1,1);
    gt_tracks{1} = in_tracks;
else
    gt_tracks = in_tracks;
end
for i=1:numel(gt_tracks)
    disp(i);
    particleCl = [0 150 0];
    RegLineCl = clMap(i,:);
    fms = gt_tracks{i}(:,1);
    % continuous path
    node = find(fms==timePt);% at most one
    if ~isempty(node)
        if node==1
            pt = round(gt_tracks{i}(1, 2:4));
            pt = max(pt, 1);
            pt = min(pt, [h,w,slice]);
            orgIm3d = draw3Dparticle(orgIm3d,  pt, particleSize, stoppedParticleCl);
        else
            for j = node-1:-1:1
                tPt = gt_tracks{i}(j, 2:4);
                tPt = max(tPt, 1);
                tPt = min(tPt, [h,w,slice]);
                hPt = gt_tracks{i}(j+1, 2:4);
                hPt = max(hPt, 1);
                hPt = min(hPt, [h,w,slice]);
                %orgIm3d = draw3Dparticle(orgIm3d,  tPt, particleSize, particleCl);
                if j+1==node
                    orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, stoppedParticleCl);
                    %orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, particleCl);
                end
                orgIm3d = draw3Dline(orgIm3d, hPt, tPt, lineWidth, RegLineCl);
            end
        end
    end
end


end