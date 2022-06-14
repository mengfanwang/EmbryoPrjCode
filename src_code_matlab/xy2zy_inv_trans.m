function vid_xyz = xy2zy_inv_trans(vid_z, zyx_flag)
% re-arrange the data from z,y,x or z,x,y back to x,y,z
% INPUT:
% vid_z: 3d data with z the first dimension
% zyx_flag: true means y is first and then x dimension, false otherwise
% OUTPUT:
% vid_xyz: 3d data with xyz dimension order

% contact: ccwang@vt.edu, 02/06/2020


if zyx_flag
    [z,h,w] = size(vid_z);
    vid_xyz = zeros(h,w,z);
    for i=1:z
        s = squeeze(vid_z(i,:,:));
        vid_xyz(:,:,i) = s;
    end
else
    [z,w,h] = size(vid_z);
    vid_xyz = zeros(h,w,z);
    for i=1:z
        s = squeeze(vid_z(i,:,:));
        vid_xyz(:,:,i) = s';
    end
end
