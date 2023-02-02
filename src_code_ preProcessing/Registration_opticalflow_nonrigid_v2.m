%% system and path
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
    data_folder = '/work/Mengfan/Embryo/22-01-11/sameViewFusion_050-149_11';
    fore_folder = '/work/public/sameViewFusion/sameViewDetection_050-149_11/synQuant_refine_res';
end

reg_ind = 51;
data1 = tifread(fullfile(data_folder, ['000' num2str(reg_ind) '.tif']));
data2 = tifread(fullfile(data_folder, ['000' num2str(reg_ind+1) '.tif']));

fore2 = load(fullfile(fore_folder, ['000' num2str(reg_ind+1) '.mat']));
fore2 = fore2.refine_res>0;

data1 = data1(:,:,1:160);
data1_backup = imresize3(data1, [960 960 160]);
data2 = data2(:,:,1:160); 
data2_backup = imresize3(data2, [960 960 160]);
fore2 = imresize3(fore2, [960 960 160]);
fore2_backup = fore2(:,:,1:160);
% fore2_backup = ones(size(data2_backup));
% fore2_backup(240:270,830:860,1:11) = 0;
%%             size          block
% layer 0: 960 960 160       512
% layer 1: 480 480 80        64
% layer 2: 240 240 40        8
% layer 3: 120 120 20        1
% layer 4: 60  60  10      
%%
clc;
clearvars -except data1_backup data2_backup fore2_backup reg_ind

sigma_gaussian = 1;
layer_num = 4;   
batch_size = round(size(data1_backup)/2^layer_num);
pad_size = [15 15 5];
step = 1;
smoothness_para = 0.1;
% resize_method = 'nearest';
resize_method = 'linear';
% resize_method = 'cubic';
gradient_method = 1;

% 10-connectivity
xyz_direction = zeros(10,4);
dir_ind = 0;
for xx = -1:1
    for yy = -1:1
        if xx~=0 || yy~=0
            dir_ind = dir_ind + 1;
            xyz_direction(dir_ind,:) = [xx yy 0 1/sqrt(xx^2+yy^2)];
        end
    end
end
for zz = -1:2:1
    dir_ind = dir_ind + 1;
    xyz_direction(dir_ind,:) = [0 0 zz 1/3];
end

tic;
for layer = layer_num:-1:0

    %     data1 = imresize3(data1_backup, round([960 960 160]/2^layer));
    %     data2 = imresize3(data2_backup, round([960 960 160]/2^layer));
    
    data1 = imgaussfilt(data1_backup,sigma_gaussian*(layer+1));
    data1 = imresize3(data1, round([960 960 160]/2^layer));
    data2 = imgaussfilt(data2_backup,sigma_gaussian*(layer+1));
    data2 = imresize3(data2, round([960 960 160]/2^layer));    

    fore2 = imresize3(fore2_backup, round([960 960 160]/2^layer));  
    fore2 = fore2 >= 0.5;

    [x,y,z,t] = size(data1);

    data1_pad = padarray(data1,pad_size,'replicate'); 
    data2_pad = padarray(data2,pad_size,'replicate');
    
    gt2 = data2;

    if layer == layer_num
        batch_num = 1;
        batch_bin = logical([1 1 1]');
        batch_loc = [1 x 1 y 1 z];
        L = gpuArray(zeros(3,3));
        smoothness = 0;
    else
        % build graph relationship
        batch_ind = 0;
        batch_table = zeros(2^(3*(layer_num-layer)), 3);
        batch_loc = zeros(2^(3*(layer_num-layer)), 6);   %[x_s x_e y_s ...]
        batch_bin_ind = 0;
        batch_bin = zeros(2^(3*(layer_num-layer)),1);    % non-zero batch index
        % acquire batch index
        for zz = 1:2^(layer_num - layer)
            for yy = 1:2^(layer_num - layer)
                for xx = 1:2^(layer_num - layer)  % zyx order is important
                    x_start = (xx-1)*batch_size(1) + 1;
                    x_end = xx*batch_size(1);
                    y_start = (yy-1)*batch_size(2) + 1;
                    y_end = yy*batch_size(2);
                    z_start = (zz-1)*batch_size(3) + 1;
                    z_end = zz*batch_size(3);
                    batch_fore = fore2(x_start:x_end, y_start:y_end, z_start:z_end);
                    batch_bin_ind = batch_bin_ind + 1;
                    if sum(batch_fore(:)) > 200
                        batch_ind = batch_ind + 1;
                        batch_table(batch_ind,:) = [xx yy zz];
                        batch_loc(batch_ind,:) = [x_start x_end y_start y_end z_start z_end];
                        batch_bin(batch_bin_ind) = 1;
                    end                  
                end
            end
        end
        batch_bin = repmat(batch_bin,1,3)';
        batch_bin = logical(batch_bin(:));
        batch_num = batch_ind;
        batch_table = batch_table(1:batch_num,:);
        batch_loc = batch_loc(1:batch_num,:);
        % acquire batch relationship
        edge_ind = 0;
        batch_relation = zeros(10*2^(3*(layer_num-layer)),3);
        for ii = 1:size(batch_table,1)
            for dir_ind = 1:size(xyz_direction,1)
                ii_nei = find(ismember(batch_table,batch_table(ii,:)+xyz_direction(dir_ind,1:3),'rows'));
                if ~isempty(ii_nei) && ~ismember([ii*3 ii_nei*3],batch_relation(:,1:2), 'rows') ...
                        && ~ismember([ii_nei*3 ii*3], batch_relation(:,1:2), 'rows')
                    edge_ind = edge_ind + 1;
                    batch_relation((edge_ind-1)*3+1:edge_ind*3,:) = [(ii-1)*3+1 (ii_nei-1)*3+1 xyz_direction(dir_ind,4);...
                        (ii-1)*3+2 (ii_nei-1)*3+2 xyz_direction(dir_ind,4); ii*3 ii_nei*3 xyz_direction(dir_ind,4)];
                end
            end
        end
        smoothness = smoothness_para*numel(data1)/sum(batch_relation(:,3));
        batch_relation = batch_relation(1:edge_ind*3,:);
%         batch_graph = graph(batch_relation(:,1), batch_relation(:,2), batch_relation(:,3));
%         if batch_graph.numnodes < batch_num*3
%             batch_graph = batch_graph.addnode(batch_num*3 - batch_graph.numnodes);
%         end
%         L = gpuArray(full(laplacian(batch_graph)));
        L = zeros(3*batch_num,3*batch_num);
        for ee = 1:edge_ind*3
            L(batch_relation(ee,1),batch_relation(ee,2)) = -batch_relation(ee,3);
            L(batch_relation(ee,2),batch_relation(ee,1)) = -batch_relation(ee,3);
        end
        for ii = 1:3*batch_num
            L(ii,ii) = -sum(L(ii,:));
        end
    end

if layer == layer_num
    phi_current = gpuArray(zeros(x,y,z,3));
    phi_current_vec = zeros(3,1);
else 
%     phi_current = zeros(x,y,z,3);
%     phi_current_vec = zeros(batch_bin_ind ,3);
    
    phi_current_vec = reshape(phi_current_vec,3,[]);
    phi_current = imresize4d(phi_current_vec, [x y z], resize_method);
    phi_current_vec_temp = imresize4d(phi_current_vec, 2, resize_method);
    phi_current_vec = zeros(batch_bin_ind ,3);
    phi_current_vec(:,1) = reshape(phi_current_vec_temp(:,:,:,1),[],1);
    phi_current_vec(:,2) = reshape(phi_current_vec_temp(:,:,:,2),[],1);
    phi_current_vec(:,3) = reshape(phi_current_vec_temp(:,:,:,3),[],1);
    phi_current_vec = phi_current_vec';
    phi_current_vec = phi_current_vec(:);
    
    phi_current = gpuArray(phi_current);
end

% phi_current = gpuArray(phi_initial(:,:,:,:,tt));
% [x_ind,y_ind,z_ind] = ind2sub(size(data1),find(fore2));
[x_ind,y_ind,z_ind] = ind2sub(size(data1),find(ones(size(fore2))));

loss = gpuArray(zeros(100000,1));
time = gpuArray(zeros(100000,1));

for iter = 1:50
    phi_previous = phi_current;
    x_bias = phi_previous(:,:,:,1);
    y_bias = phi_previous(:,:,:,2);
    z_bias = phi_previous(:,:,:,3);


    x_new = x_ind + x_bias(:);
    y_new = y_ind + y_bias(:);
    z_new = z_ind + z_bias(:);
    data1_tran = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));

    if gradient_method == 1
    % method1 to calculate gradient
    data1_x_incre = interp3(data1_pad,y_new+pad_size(2),x_new+step+pad_size(1),z_new+pad_size(3));
    data1_x_decre = interp3(data1_pad,y_new+pad_size(2),x_new-step+pad_size(1),z_new+pad_size(3));
    Ix = (data1_x_incre - data1_x_decre)/(2*step);
    Ix = reshape(Ix,[x y z]);
    clear data1_x_incre data1_x_decre

    data1_y_incre = interp3(data1_pad,y_new+step+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
    data1_y_decre = interp3(data1_pad,y_new-step+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
    Iy = (data1_y_incre - data1_y_decre)/(2*step);
    Iy = reshape(Iy,[x y z]);
    clear data1_y_incre data1_y_decre  

    data1_z_incre = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+step+pad_size(3));
    data1_z_decre = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new-step+pad_size(3));
    Iz = (data1_z_incre - data1_z_decre)/(2*step);
    Iz = reshape(Iz,[x y z]);
    clear data1_z_incre data1_z_decre      
    
    elseif gradient_method == 2
    % method2 to calculate gradient
    Ix = reshape(data1_tran,[x y z]);
    Ix = gradient3(Ix, 'x');
    
    
    Iy = reshape(data1_tran,[x y z]);
    Iy = gradient3(Iy, 'y');
    
  
    Iz = reshape(data1_tran,[x y z]);
    Iz = gradient3(Iz, 'z');
    end

    It = data1_tran-gt2(:);
    It = reshape(It,[x y z]);

    AtA = zeros(batch_num*3,batch_num*3);
    AtIt = zeros(batch_num*3,1);
    for ii = 1:batch_num
        batch_loc_temp = num2cell(batch_loc(ii,:));
        [x_start, x_end, y_start, y_end, z_start, z_end] = deal(batch_loc_temp{:});
        batch_Ix = Ix(x_start:x_end, y_start:y_end, z_start:z_end);
        batch_Iy = Iy(x_start:x_end, y_start:y_end, z_start:z_end);
        batch_Iz = Iz(x_start:x_end, y_start:y_end, z_start:z_end);
        batch_It = It(x_start:x_end, y_start:y_end, z_start:z_end);
        
        batch_fore = fore2(x_start:x_end, y_start:y_end, z_start:z_end);
        batch_Ix = batch_Ix(batch_fore);
        batch_Iy = batch_Iy(batch_fore);
        batch_Iz = batch_Iz(batch_fore);
        batch_It = batch_It(batch_fore);

        A_ii = [batch_Ix(:) batch_Iy(:) batch_Iz(:)];
        AtA((ii-1)*3+1:ii*3, (ii-1)*3+1:ii*3) = A_ii'*A_ii;
        AtIt((ii-1)*3+1:ii*3, 1) = A_ii'*batch_It(:);
    end
    AtA = gpuArray(AtA);
    AtIt = gpuArray(AtIt);
%     phi_gradient = -(AtA + smoothness*L)\AtIt;
    phi_gradient = -(AtA + smoothness*L)\(AtIt + smoothness*L*phi_current_vec(batch_bin));
    if max(abs(phi_gradient)) < 1e-4 || norm(phi_gradient) < 1e-6
        break;
    end
    
%     fprintf('Gradient: %f %f %f\n',phi_gradient(1),phi_gradient(2),phi_gradient(3));
% method 1: local patch rigid transform
%     for ii = 1:batch_num
%         batch_loc_temp = num2cell(batch_loc(ii,:));
%         [x_start, x_end, y_start, y_end, z_start, z_end] = deal(batch_loc_temp{:});
%         phi_current(x_start:x_end,y_start:y_end,z_start:z_end,1) = ...
%             phi_current(x_start:x_end,y_start:y_end,z_start:z_end,1) + phi_gradient((ii-1)*3+1);
%         phi_current(x_start:x_end,y_start:y_end,z_start:z_end,2) = ...
%             phi_current(x_start:x_end,y_start:y_end,z_start:z_end,2) + phi_gradient((ii-1)*3+2);
%         phi_current(x_start:x_end,y_start:y_end,z_start:z_end,3) = ...
%             phi_current(x_start:x_end,y_start:y_end,z_start:z_end,3) + phi_gradient(ii*3);
%     end
% method 2:interpolation
    phi_current_vec(batch_bin) = phi_current_vec(batch_bin) + phi_gradient;
    phi_current = imresize4d(reshape(phi_current_vec, 3, []), [x y z], resize_method);
    

    mse = mean((data1_tran(fore2) - gt2(fore2)).^2);


    if batch_num == 1
        smooth_error = 0;
%         smooth_error2 = 0;
    else
%         smooth_error = 0;
%         phi_x = unique(phi_current(:,:,:,1),'stable');
%         phi_y = unique(phi_current(:,:,:,2),'stable');
%         phi_z = unique(phi_current(:,:,:,3),'stable');
%         phi_x = unique(gather(phi_current(:,:,:,1)),'stable');
%         phi_y = unique(gather(phi_current(:,:,:,2)),'stable');
%         phi_z = unique(gather(phi_current(:,:,:,3)),'stable');
%         phi_all = [phi_x'; phi_y'; phi_z';];
%         smooth_error = smoothness*phi_all(:)'*L*phi_all(:)/numel(data1_tran);

% %         % second method calculate smooth error
%         smooth_error2 = 0;
%         if ~exist('phi_all2', 'var')
%             phi_all2 = zeros(size(L,1),1);
%         end
%         phi_all2 = phi_all2 + phi_gradient;
%         for ee = 1:size(batch_relation,1)
%             ii = batch_relation(ee,1);
%             jj = batch_relation(ee,2);
%             smooth_error2 = smooth_error2 + (phi_all2(ii)-phi_all2(jj)).^2;
%         end
%         smooth_error2 = smooth_error2*smoothness/length(data1_tran(fore2));
        
        smooth_error = smoothness*phi_current_vec(batch_bin)'*L*phi_current_vec(batch_bin)/sum(fore2(:));


    end
    fprintf('Gradient: %f\n',norm(phi_gradient)/length(phi_gradient));
    fprintf('Iteration %d\n Current error:%f MSE: %f Smooth: %f Time:%f\n',iter, mse+smooth_error, mse, smooth_error, toc);
    loss(iter) = mse+smooth_error;
    time(iter) = toc;
end
    mse_ori = mean((data1(fore2) - gt2(fore2)).^2);
    fprintf('Original error: %f\n', mse_ori);
    loss = gather(loss(1:iter));
    time = gather(time(1:iter));
end
% phi_current = gather(phi_current);
% save(['/work/Mengfan/Embryo/Registration/nonrigid_result/' num2str(reg_ind) '_l' num2str(layer_num) '_s'  num2str(smoothness_para)  '_' resize_method '_g' num2str(gradient_method) '.mat'], 'phi_current_vec', 'loss', 'mse_ori', 'time' ,'-v7.3');

function phi_current = imresize4d(phi_current, scale, method)
    % transform to 4d data if table
    % resize 4d data
    % gpu not supported
    flag = 0;
    if isgpuarray(phi_current)
        isgpu = 1;
        phi_current = gather(phi_current);
    else
        isgpu = 0;
    end
    if ismatrix(phi_current)       % if table
        % change to [3 []]
        if size(phi_current,1) == 3
            flag = 1;
        elseif size(phi_current, 2) == 3
            flag = 1;
            phi_current = phi_current';
        end
        if flag
            batch_dim = round(size(phi_current,2)^(1/3));
            if batch_dim^3 ~= size(phi_current,2)
                error('Incorrect number of elements.');
            end
            if length(scale) == 1
                scale = [batch_dim batch_dim batch_dim] * scale;
            end
            phi_current_temp = phi_current';
            phi_current = zeros([scale 3]);
            if size(phi_current_temp,1) == 1
                phi_current(:,:,:,1) = phi_current_temp(:,1);
                phi_current(:,:,:,2) = phi_current_temp(:,2);
                phi_current(:,:,:,3) = phi_current_temp(:,3);
            else
                phi_current(:,:,:,1) = imresize3(reshape(phi_current_temp(:,1), batch_dim, batch_dim, batch_dim), scale, method);
                phi_current(:,:,:,2) = imresize3(reshape(phi_current_temp(:,2), batch_dim, batch_dim, batch_dim), scale, method);
                phi_current(:,:,:,3) = imresize3(reshape(phi_current_temp(:,3), batch_dim, batch_dim, batch_dim), scale, method);
            end
        end
    elseif  ndims(phi_current) == 4
        if size(phi_current,4) == 3
            flag = 1;
            if length(scale) == 1
                scale = size(phi_current(:,:,:,1)) * scale;
            end
            phi_current_temp = phi_current;
            phi_current = zeros([scale 3]);
            phi_current(:,:,:,1) = imresize3(phi_current_temp(:,:,:,1), scale, method);
            phi_current(:,:,:,2) = imresize3(phi_current_temp(:,:,:,2), scale, method);
            phi_current(:,:,:,3) = imresize3(phi_current_temp(:,:,:,3), scale, method);
        end
    end
    if isgpu
        phi_current = gpuArray(phi_current);
    end
    if flag == 0
        warning('Incorrect input. No operations.');
    end
    
end