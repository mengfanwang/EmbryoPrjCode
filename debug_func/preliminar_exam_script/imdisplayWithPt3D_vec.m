function orgIm3d = imdisplayWithPt3D_vec(orgIm3d, gt_tracks, timePt, p)
% generate a 3D colorful data with ROI overlay
% gt_tracks: t, y, x, z
particleSize = p.particleSize;
%particleCl = p.particleCl;
clMap = p.cmap;
[h,w,~,slice] = size(orgIm3d);
for i=1:numel(gt_tracks)
    disp(i);
    particleCl = clMap(i,:);
    fms = gt_tracks{i}(:,1);
    % continuous path
    node = find(fms==timePt);% at most one
    if ~isempty(node)
        pt = round(gt_tracks{i}(node, 2:4));
        pt = max(pt, 1);
        pt = min(pt, [h,w,slice]);
        orgIm3d = draw3Dparticle(orgIm3d,  pt, particleSize, particleCl);
    end
end


end