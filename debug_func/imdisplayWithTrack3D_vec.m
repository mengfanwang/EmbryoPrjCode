function orgIm3d = imdisplayWithTrack3D_vec(orgIm3d, gt_tracks, timePt, p)
% generate a 3D colorful data with ROI overlay
% gt_tracks: t, y, x, z
particleSize = p.particleSize;
%particleCl = p.particleCl;
stoppedParticleCl = p.stoppedParticleCl;
lineWidth = p.lineWidth;
RegLineCl = p.RegLineCl;
stoppeLineCl = p.stoppeLineCl;
clMap = p.cmap;
[h,w,~,slice] = size(orgIm3d);
cnt = 0;
for i=randperm(numel(gt_tracks))
    disp(i);
    particleCl = [0 150 0];
    RegLineCl = clMap(i,:);
    fms = gt_tracks{i}(:,1);
    % continuous path
    node = find(fms==timePt);% at most one
    if ~isempty(node)
        cnt = cnt + 1;
        if cnt > 30
            break;
        end
        if size(gt_tracks{i}, 1)==1
            pt = round(gt_tracks{i}(2:4));
            pt = max(pt, 1);
            pt = min(pt, [h,w,slice]);
            orgIm3d = draw3Dparticle(orgIm3d,  pt, particleSize, particleCl);
        else
            for j = node:size(gt_tracks{i},1)-1
                tPt = gt_tracks{i}(j, 2:4);
                tPt = max(tPt, 1);
                tPt = min(tPt, [h,w,slice]);
                hPt = gt_tracks{i}(j+1, 2:4);
                hPt = max(hPt, 1);
                hPt = min(hPt, [h,w,slice]);
                %orgIm3d = draw3Dparticle(orgIm3d,  tPt, particleSize, particleCl);
                if j==node
                    orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, stoppedParticleCl);
                    %orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, particleCl);
                end
                orgIm3d = draw3Dline(orgIm3d, hPt, tPt, lineWidth, RegLineCl);
            end
        end
    else
        % if the path is jumping this frame
        if sum(fms>timePt) && sum(fms<timePt)
            node = find(fms<timePt);
            node = node(end);
            if node==1
                hPt = round(gt_tracks{i}(1, 2:4));
                hPt = max(hPt, 1);
                hPt = min(hPt, [h,w,slice]);
                orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, stoppedParticleCl);
            else
                for j = node:-1:2
                    tPt = gt_tracks{i}(j-1, 2:4);
                    tPt = max(tPt, 1);
                    tPt = min(tPt, [h,w,slice]);
                    hPt = gt_tracks{i}(j, 2:4);
                    hPt = max(hPt, 1);
                    hPt = min(hPt, [h,w,slice]);   
                    orgIm3d = draw3Dparticle(orgIm3d,  tPt, particleSize, stoppedParticleCl);
                    if j==node
                        orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, stoppedParticleCl);
                    end
                    orgIm3d = draw3Dline(orgIm3d, hPt, tPt, lineWidth, stoppeLineCl);
                end
            end
        end
    end
end


end