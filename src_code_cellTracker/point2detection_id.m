function [id, ids] = point2detection_id(point, label_maps)
% point: t, y, x, z

idx = get_index_from_center(point(2:4), label_maps{1});

ids = label_maps{point(1)}(idx);

ids = ids(ids>0);
if ~isempty(ids)
    id  = mode(ids);
else
    id = nan;
end


end