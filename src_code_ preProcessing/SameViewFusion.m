clc;clear;close all;
% draft script
% fuse the two image from the same angle

im1 = tifread('/work/Mengfan/Embryo/22-01-11/data/00248/10.tif');
im2 = tifread('/work/Mengfan/Embryo/22-01-11/data/00248/14.tif');
[x, y, z] = size(im1);
save_folder = '/work/Mengfan/Embryo/22-01-11/data_merge/00248';


% 
im_fuse = zeros(size(im1));
for zz = 1:z
    im_fuse(:,:,zz) = wfusimg(im1(:,:,zz),im2(:,:,zz), 'db4', 5, 'mean', 'max');
end
tifwrite(uint16(im_fuse), fullfile(save_folder, 'wavelet_paper'));