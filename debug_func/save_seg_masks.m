function save_seg_masks(vid, mask, save_path_and_name)
%DISPLAY_MASKS Summary of this function goes here
%   Detailed explanation goes here

g = scale_image(vid,0,1);
r = g; r(mask<=0) = 0;
b = g;
[h,w,slice] = size(r);
zz = zeros(h,w,3,slice);
for i=1:slice
    zz(:,:,1,i) = r(:,:,i);
    zz(:,:,2,i) = g(:,:,i);
    zz(:,:,3,i) = b(:,:,i);
end

tifwrite(label2rgb3d(mask), [save_path_and_name, '_mask_in_color'])
tifwrite(zz, save_path_and_name);
