function zz = display_seg_res4PE(vid, mask)
% display the cell segmentation results

[h,w,zslice] = size(vid);
zz = zeros(h,w, 3,zslice);
vid = scale_image(vid,0,1)*1;
for i=1:zslice
    %zz(:,:,1,i) = mask(:,:,i)*0.5;
    zz(:,:,2,i) = vid(:,:,i);
    zz(:,:,3,i) = (mask(:,:,i)>0)*0.5;
end
zzshow(zz); 