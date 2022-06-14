function neighbor_cells = neighbor_pixels(fg_map, im_sz, gap)
% find the neighbor locations for a given area
%INPUT:
% fg_map: the map indicating the area or a vector indicating the locations
% (index)
% im_sz: size of the image [h,w] or [h,w,z]
% gap: how many circles we want, 1: the most close neighbors, 2: two
% circles neighbors... default is eight neighbors
%OUTPUT:
% neighbor_cells: cells contain the location of each circle

% contact: ccwang@vt.edu, 12/04/2019

if size(fg_map,1)==1 || size(fg_map,2)==1                                   % if fg_map is in list
    idx = fg_map;
    fg_map = zeros(im_sz);
    fg_map(idx) = 1;
end
neighbor_cells = cell(gap, 1);
if length(im_sz)==2                                                         % consider 2D 8-neighbor only
    nhood = [1 1 1; 1 0 1; 1 1 1];
else % 3d image
%      nhood = strel('cube',3);
    nhood = zeros(3,3,3);
    nhood(:,:,2) = [1 1 1; 1 0 1; 1 1 1];
end
for i=1:gap
    d_map = imdilate(fg_map,nhood);
    
    neighbor_cells{i} = find(d_map-fg_map>0);
    fg_map = d_map;
end
end