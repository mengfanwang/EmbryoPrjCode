function display_rgb_time_course(refine_res)
% timelapse 3d data saved into cells, each of which is a frame
max_proj = cellfun(@(x) max(x, [], 3), refine_res,'UniformOutput', false);

max_proj = cat(3, max_proj{:});
zzshow(label2rgb3d(max_proj));
end