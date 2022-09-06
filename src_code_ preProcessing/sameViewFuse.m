function im_fuse = sameViewFuse(im1, im2)
% Using wavelet transform to fuse two images in the same view
% We'd suggest parallel computation for 10x faster acceleration
im_fuse = zeros(size(im1));
parfor zz = 1:size(im_fuse,3)
    im_fuse(:,:,zz) = wfusimg(im1(:,:,zz),im2(:,:,zz), 'db4', 5, 'mean', 'max');
end

