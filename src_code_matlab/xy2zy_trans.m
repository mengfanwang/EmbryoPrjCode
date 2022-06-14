function vid_z = xy2zy_trans(vid_xyz, zyx_flag)
% re-arrange the data from y,x,z to z,y,x or z,x,y
% INPUT:
% vid_xyz: input 3d data
% zyx_flag: true means y is first and then x dimension, false otherwise
% OUTPUT:
% vid_z: output 3d data with z the first dimension

% contact: ccwang@vt.edu, 02/06/2020

[h,w,z] = size(vid_xyz);
if zyx_flag
    vid_z = zeros(z,h,w);
    for i=1:w
        s = squeeze(vid_xyz(:,i,:));
        vid_z(:,:,i) = s';
    end
else
    vid_z = zeros(z,w,h);
    for i=1:h
        s = squeeze(vid_xyz(i,:,:));
        vid_z(:,:,i) = s';
    end
end
