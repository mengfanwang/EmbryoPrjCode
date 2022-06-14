function [dat_in, excess_node, c_en, c_ex] = build_embryo_graph(detections)
[lens,~] = cellfun(@size, detections);
incre = cumsum(lens);
det_incre = [0;incre(1:end-1)];


nextDetections = detections(2:end);
next_incre = det_incre(2:end);


T = numel(detections);

neiCells = cell(T,1);
parfor i=1:T-1
    disp(i);
    curDet = detections{i};
    nextDet = nextDetections{i};
    
    incre_i = det_incre(i);
    incre_j = next_incre(i);
    nei = cell(size(curDet,1), 1);
    for j=1:size(curDet,1)
        dist0 = vecnorm((curDet(j,1:3)-nextDet(:,1:3))');
        [zz,od] = sort(dist0);
        nei{j} = [od(1:3)' + incre_j, zz(1:3)'];
    end
    neiCells{i} = nei;
end


% calculate en and ex cost
en_num = lens(1);
ex_num = lens(end);
% en_num = 0;
% ex_num = 0;
for i=2:T
    if lens(i)-lens(i-1) > 0
        en_num = en_num + abs(lens(i)-lens(i-1));
    else
        ex_num = ex_num + abs(lens(i)-lens(i-1));
    end
    
end
en_cost = -log(en_num/sum(lens));
ex_cost = en_cost;
c_i = -(en_cost + ex_cost);



%features.det = zeros(sum(lens),1) - (en_cost+ex_cost);
features.nei = cat(1, neiCells{:});
num_det = sum(lens);
nn_dist = zeros(num_det,1);
for i=1:numel(features.nei)
    if mod(i,1e6)==0
        disp(i);
    end
    if ~isempty(features.nei{i}) 
        nn_dist(i) = features.nei{i}(1,2);
    end
end


nn_dist(isnan(nn_dist)) = [];
tmp_nd = round(nn_dist * 100);
p = zeros(max(tmp_nd),1);
len_dist = length(tmp_nd);

parfor i=1:max(tmp_nd)
    p(i) = -log(1-sum(tmp_nd<=i)./len_dist);
end


n_nodes = 2*num_det+2; 


dat_in = zeros(1e8,3); %% each row represents an edge from node in column 1 to node in column 2 with cost in column 3.
k_dat = 0;
for i = 1:num_det
    k_dat = k_dat+3;
    dat_in(k_dat-2,:) = [1      2*i     en_cost      ];
    dat_in(k_dat-1,:) = [2*i    2*i+1   c_i       ];
    dat_in(k_dat,:)   = [2*i+1  n_nodes ex_cost      ];
end

for i=1:numel(features.nei)
    f2 = features.nei{i};
    for j = 1:size(f2,1)
        p_loc = round(f2(j,2)*100);
        if p_loc<1
            p_loc = 1;
        end
        if p_loc > length(p)
            p_loc = length(p);
        end
        c_ij = p(p_loc);
        k_dat = k_dat+1;
        dat_in(k_dat,:) = [2*i+1 2*f2(j,1) c_ij];
    end
end


dat_in = dat_in(1:k_dat,:);


dat_in = [dat_in repmat([0 1],size(dat_in,1),1)];  %% add two columns: 0 for min capacity in column 4 and 1 for max capacity in column 5 for all edges.



excess_node = [1 n_nodes];  %% push flow in the first node and collect it in the last node.
c_en = en_cost;
c_ex = ex_cost;
%save('embryo_10tb_heuristic_graph.mat', 'dat_in', 'excess_node','c_en', 'c_ex', '-v7.3');


end