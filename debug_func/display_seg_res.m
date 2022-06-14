function zz = display_seg_res(movieInfo, refine_res, embryo_vid, cell_id, threshold_res)
% display the cell segmentation results

voxIdx = movieInfo.voxIdx{cell_id};
frame = movieInfo.frames(cell_id);

vid = crop3D(embryo_vid{frame}, voxIdx, [20 20 5]);
idComp = crop3D(refine_res{frame}, voxIdx, [20 20 5]);

[h,w,zslice] = size(vid);
sMap = idComp == refine_res{frame}(voxIdx(1));

zz = zeros(h,w, 3,zslice);
vid = scale_image(vid,0,1);
for i=1:zslice
    zz(:,:,1,i) = sMap(:,:,i)*0.5;
    zz(:,:,2,i) = vid(:,:,i);
    zz(:,:,3,i) = (idComp(:,:,i)>0)*0.5;
end
if nargin == 4
    zzshow(zz); hold on; title(['frame:', num2str(frame), ...
        ', id:', num2str(cell_id)]);
else
    tt = unique(threshold_res{frame}(voxIdx));
    zzshow(zz); hold on; title(['frame:', num2str(frame), ...
        ', id:', num2str(cell_id), ', thres: ', num2str(tt')]);
end