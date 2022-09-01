clc;clear;close all;
% draft script
% fuse the two image from the same angle

im1 = tifread('/work/Mengfan/Embryo/22-01-11/data/00248/10.tif');
% im2 = tifread('/work/Mengfan/Embryo/22-01-11/data/00248/14.tif');
[x, y, z] = size(im1);

% im_fuse = zeros(size(im1));
% for zz = 1:z
%     zz
%     im_fuse(:,:,zz) = wfusimg(im1(:,:,zz),im2(:,:,zz), 'db4', 5, 'mean', 'max');
% end