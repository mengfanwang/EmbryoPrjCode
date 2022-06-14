function display_gap_res(movieInfo, refine_res, embryo_vid, eigMaps, ...
    cell_id, gap3d)
% display the cell segmentation results
if nargin ==5
    gap3d = true;
end
voxIdx = movieInfo.voxIdx{cell_id};
frame = movieInfo.frames(cell_id);

vid = crop3D(embryo_vid{frame}, voxIdx, [10 10 5]);
idComp = crop3D(refine_res{frame}, voxIdx, [10 10 5]);
if gap3d
    eigComp = crop3D(eigMaps{frame}{2}, voxIdx, [10 10 5]);
else % 2d principla curvature
    eigComp = crop3D(eigMaps{frame}{1}, voxIdx, [10 10 5]);
end
[h,w,zslice] = size(vid);
sMap = idComp == refine_res{frame}(voxIdx(1));

zz = zeros(h,w, 3,zslice);
vid = scale_image(vid,0,1);
for i=1:zslice
    zz(:,:,1,i) = (eigComp(:,:,i)>0 & sMap(:,:,i)>0)*0.5;
    zz(:,:,2,i) = vid(:,:,i);
    zz(:,:,3,i) = (sMap(:,:,i)>0 & eigComp(:,:,i)<=0)*0.5;
end
zzshow(zz); hold on; title(['frame:', num2str(frame), ...
    ', id:', num2str(cell_id)]);
end