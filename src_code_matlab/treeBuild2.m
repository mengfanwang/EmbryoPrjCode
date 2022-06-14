function [one_cell_flag, oneCellpValue, mergeVSseg] = treeBuild2(root_node, movieInfo, g, dir_idx)
% second way to build a tree given a root node
% dir_idx:  0: means only build using the decendents of root_node
%           1: means only build using the ancestors of root_node
%           2: build both trees, one forward tree and one backward tree
% NOTE: we need verify we are still working on the same cell(s), if only
% one region exist in a frame; TODO: Do we need to consider continutity?

right_cnt = [0 0]; % one, >one
if dir_idx==0 || dir_idx == 2
    % right-direction
%     stack_nodes = root_node;
%     rootVox = movieInfo.vox{root_node};
%     tree_nodes = [];
%     multi_pre_flag = false; % if there is a frame lose one region, we need to test
%     while ~isempty(stack_nodes)
%         tree_nodes = cat(1, tree_nodes, stack_nodes);
%         if length(stack_nodes)==1
%             if stack_nodes ~= root_node && multi_pre_flag
%                 curVox = movieInfo.vox{stack_nodes};
%                 prePostVox = cell(2,1);
%                 prePostVox{1} = cat(1,movieInfo.vox{movieInfo.parents{stack_nodes}});
%                 prePostVox{2} = cat(1,movieInfo.vox{movieInfo.kids{stack_nodes}});
%                 same_flag = verifySameCell(rootVox, curVox, prePostVox, ...
%                     movieInfo.ovGamma, min(g.c_ex, g.c_en));
%                 if ~same_flag
%                     break;
%                 end
%             end
%         end
%         further_kids = [];
%         for j=1:length(stack_nodes)
%             further_kids = cat(1, further_kids, movieInfo.kids{stack_nodes(j)});
%         end
%         if length(stack_nodes) > 1
%             multi_pre_flag = true;
%         else
%             multi_pre_flag = false;
%         end
%         stack_nodes = unique(further_kids);
%     end
    tree_nodes = [];
    stack_nodes = root_node;
    while ~isempty(stack_nodes)
        tree_nodes = cat(1, tree_nodes, stack_nodes);
        further_kids = [];
        for j=1:length(stack_nodes)
            further_kids = cat(1, further_kids, movieInfo.kids{stack_nodes(j)});
        end
        stack_nodes = further_kids;
    end
    %tree_nodes = sort(tree_nodes,'ascend');
    tree_nodes = unique(tree_nodes);% already sorted with ascend order
    tree_frames = movieInfo.frames(tree_nodes);
    freq_cnts = frequency_cnt(tree_frames);
    right_cells = mat2cell(tree_nodes, freq_cnts(:,2),1);
    % test if should remove some levels of the tree
    multi_pre_flag = false;
    rootVox = movieInfo.vox{root_node};
    for i=2:numel(right_cells)% if there is a frame lose one region, we need to test
        if multi_pre_flag && length(right_cells{i}) == 1
            curVox = movieInfo.vox{right_cells{i}};
            prePostVox = cell(2,1);
            prePostVox{1} = cat(1,movieInfo.vox{movieInfo.parents{right_cells{i}}});
            prePostVox{2} = cat(1,movieInfo.vox{movieInfo.kids{right_cells{i}}});
            same_flag = verifySameCell(rootVox, curVox, prePostVox, ...
                movieInfo.ovGamma, min(g.c_ex, g.c_en));
            if ~same_flag
                right_cells = right_cells(1:i-1);
                break;
            end
        end
        if length(right_cells{i}) > 1
            multi_pre_flag = true;
        else
            multi_pre_flag = false;
        end
    end
    freq_cnts = cellfun(@length, right_cells);
    right_cnt(1) = sum(freq_cnts==1);
    right_cnt(2) = sum(freq_cnts>1);
end
left_cnt = [0 0]; % one, >one
if dir_idx==1 || dir_idx == 2
    % left-direction
%     tree_nodes = [];
%     stack_nodes = root_node;
%     multi_post_flag = false;
%     while ~isempty(stack_nodes)
%         tree_nodes = cat(1, tree_nodes, stack_nodes);
%         if length(stack_nodes)==1
%             if stack_nodes ~= root_node && multi_post_flag
%                 curVox = movieInfo.vox{stack_nodes};
%                 prePostVox = cell(2,1);
%                 prePostVox{1} = cat(1,movieInfo.vox{movieInfo.parents{stack_nodes}});
%                 prePostVox{2} = cat(1,movieInfo.vox{movieInfo.kids{stack_nodes}});
%                 same_flag = verifySameCell(rootVox, curVox, prePostVox, ...
%                     movieInfo.ovGamma, min(g.c_ex, g.c_en));
%                 if ~same_flag
%                     break;
%                 end
%             end
%         end
%         further_parents = [];
%         for j=1:length(stack_nodes)
%             further_parents = cat(1, further_parents, movieInfo.parents{stack_nodes(j)});
%         end
%         if length(stack_nodes) > 1
%             multi_post_flag = true;
%         else
%             multi_post_flag = false;
%         end
%         stack_nodes = unique(further_parents);
%     end
    tree_nodes = [];
    stack_nodes = root_node;
    while ~isempty(stack_nodes)
        tree_nodes = cat(1, tree_nodes, stack_nodes);
        further_parents = [];
        for j=1:length(stack_nodes)
            further_parents = cat(1, further_parents, movieInfo.parents{stack_nodes(j)});
        end
        stack_nodes = further_parents;
    end
    tree_nodes = unique(tree_nodes); % ascending order
    tree_frames = movieInfo.frames(tree_nodes);
    freq_cnts = frequency_cnt(tree_frames);
    left_cells = mat2cell(tree_nodes, freq_cnts(:,2),1);
    left_cells = left_cells(end:-1:1); % reverse the cells
    % test if should remove some levels of the tree
    multi_post_flag = false;
    rootVox = movieInfo.vox{root_node};
    for i=2:numel(left_cells)% if there is a frame lose one region, we need to test
        if multi_post_flag && length(left_cells{i}) == 1
            curVox = movieInfo.vox{left_cells{i}};
            prePostVox = cell(2,1);
            prePostVox{1} = cat(1,movieInfo.vox{movieInfo.parents{left_cells{i}}});
            prePostVox{2} = cat(1,movieInfo.vox{movieInfo.kids{left_cells{i}}});
            same_flag = verifySameCell(rootVox, curVox, prePostVox, ...
                movieInfo.ovGamma, min(g.c_ex, g.c_en));
            if ~same_flag
                left_cells = left_cells(1:i-1);
                break;
            end
        end
        if length(left_cells{i}) > 1
            multi_post_flag = true;
        else
            multi_post_flag = false;
        end
    end
    freq_cnts = cellfun(@length, left_cells);
    left_cnt(1) = sum(freq_cnts==1);
    left_cnt(2) = sum(freq_cnts>1);
end

mergeVSseg = left_cnt + right_cnt; % element 1 means one cell, 2 means >1 cells
mergeVSseg(1) = mergeVSseg(1)-1;% root node was counted twice

oneCellpValue = binocdf(mergeVSseg(1), sum(mergeVSseg), 0.5);

if oneCellpValue < 0.05 % probability of one cell in this region is low
    one_cell_flag = 0;
elseif 1-oneCellpValue < 0.05 % probability of multi-cell in this region is low
    one_cell_flag = 1;
else % no promising conclusion can be made
    one_cell_flag = nan;
end
end