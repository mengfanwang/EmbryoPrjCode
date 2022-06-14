function flag = adjframesConsistency(root_node, movieInfo, ...
        refine_res, g, o2m_flag)
% if there are segmentation in consistency in continuous frames
% we should merge them first
% NOTE: we only consider two->two regions
if o2m_flag
    p1 = movieInfo.kids{root_node};
    test1 = movieInfo.kids{p1(1)};
    test2 = movieInfo.kids{p1(2)};
    test_frame = movieInfo.frames(p1(1)) + 1;
else
    p1 = movieInfo.parents{root_node};
    test1 = movieInfo.parents{p1(1)};
    test2 = movieInfo.parents{p1(2)};
    test_frame = movieInfo.frames(p1(1)) - 1;
end
flag = true;

if (~isempty(test1) && ~isempty(test2)) || test_frame<1 || test_frame>numel(refine_res)
    return;
end
test_nodes = unique(cat(1, test1, test2));
if isempty(test1) && isempty(test2)
    % check if there are two nodes in test_frame also have no parents/kids
    % if so, they may be paired with the two nodes in p1
    test1 = bestOvNei(p1(1), movieInfo, test_frame);
    test2 = bestOvNei(p1(2), movieInfo, test_frame);
    if ~isnan(test1) && ~isnan(test2)
        if o2m_flag
            test1_n = movieInfo.parents{test1};
            test2_n = movieInfo.parents{test2};
        else
            test1_n = movieInfo.kids{test1};
            test2_n = movieInfo.kids{test2};
        end
        if isempty(test1_n) && isempty(test2_n)
            test_nodes = unique(cat(1, test1, test2));
        end
    end
end

if length(test_nodes) ~= 2
    return
end

voxIdx1 = cat(1, movieInfo.voxIdx{p1});
voxIdx2 = cat(1, movieInfo.voxIdx{test_nodes});

[mc, ~] = voxIdx2cost(voxIdx1, voxIdx2, ...
         [movieInfo.frames(p1(1)), test_frame], movieInfo, size(refine_res{1}));

if mc > abs(g.observationCost / 2)
    return;
end

[mc1, ~] = voxIdx2cost(voxIdx1, movieInfo.voxIdx{test_nodes(1)}, ...
 [movieInfo.frames(p1(1)), test_frame], movieInfo, size(refine_res{1}));
[mc2, ~] = voxIdx2cost(voxIdx1, movieInfo.voxIdx{test_nodes(2)}, ...
 [movieInfo.frames(p1(1)), test_frame], movieInfo, size(refine_res{1}));
if mc < min(mc1, mc2)
    flag = false;
end