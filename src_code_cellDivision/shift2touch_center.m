function out_xyz = shift2touch_center(center, in_xyz, bound_xyz)
% Given some points, move the points to touch a given center.
% NOTE: we do not want to directly move the center of given point clouds to
% the given center. We move the point cloud along the direction from them
% to the given center and stop once they touch the center.

if ~isempty(find(isnan(center),1)) || isempty(in_xyz)
    out_xyz = in_xyz;
    return;
end

dif_xyz = in_xyz - center;
dist = sum(dif_xyz.^2, 2);

[~, od] = min(dist);

shift = dif_xyz(od,:);


out_xyz = in_xyz - shift;

if nargin == 3 % remove those out of bound and round them to integers
    out_xyz = round(out_xyz);

    valid_ones = out_xyz(:,1) <= bound_xyz(1) & ...
        out_xyz(:,2) <= bound_xyz(2) & out_xyz(:,3) <= bound_xyz(3) & ...
        out_xyz(:,1) >=1 & out_xyz(:,2) >=1 & out_xyz(:,3) >=1;
    
    out_xyz = out_xyz(valid_ones, :);
    out_xyz = unique(out_xyz, 'rows');
end