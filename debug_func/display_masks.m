function zz = display_masks(vid, mask1, mask2)
%DISPLAY_MASKS Summary of this function goes here
%   Detailed explanation goes here
if nargin == 2
    mask2 = ones(size(mask1));
end
g = scale_image(vid,0,1);
r = g; r(mask1<=0) = 0;
b = g; b(mask2<=0) = 0;
[h,w,slice] = size(r);
zz = zeros(h,w,3,slice);
for i=1:slice
    zz(:,:,1,i) = r(:,:,i);
    zz(:,:,2,i) = g(:,:,i);
    zz(:,:,3,i) = b(:,:,i);
end

zzshow(zz);
