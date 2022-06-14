function newRegMap = sliceRefine(vid, prinCv, regMap)
% refine the regions using information across z-slices
% INPUT:
% vid: input gray image data
% prinCv: the pincipla curvature map
% regMap: binary region map indicating the region location
% OUTPUT:
% newRegMap

% contact: ccwang@vt.edu, 02/04/2020


%% handle over-merge on the z-direction
[h,w,z] = size(vid);
