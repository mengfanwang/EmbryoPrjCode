%% system and path
dbstop if error
for reg_ind = 7:80
clearvars -except reg_ind
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');

    data_folder = '/work/public/Embryo/Amat2014/data';
    fore_folder = '/work/public/Embryo/Amat2014/Detection/Wei_refine_res';
    save_folder = '/work/public/Embryo/Amat2014/Registration';
end

data1_name = num2str(100000 + reg_ind);
data1_name = data1_name(2:end);
data1 = tifread(fullfile(data_folder, [data1_name '.tif']));
data2_name = num2str(100000 + reg_ind + 1);
data2_name = data2_name(2:end);
data2 = tifread(fullfile(data_folder, [data2_name '.tif']));

fore2 = load(fullfile(fore_folder, [data2_name '.mat']));
fore2 = fore2.refine_res>0;

zslices = size(data1,3);
data1_backup = imresize3(data1, [896 909 zslices]);
data2_backup = imresize3(data2, [896 909 zslices]);
fore2_backup = imresize3(fore2, [896 909 zslices]);
layer_num = 5;  
if mod(zslices, 2^layer_num) ~= 0 
    data1_backup = padarray(data1_backup, [960-896 960-909 2^layer_num-mod(zslices, 2^layer_num)], 'replicate','post');
    data2_backup = padarray(data2_backup, [960-896 960-909 2^layer_num-mod(zslices, 2^layer_num)], 'replicate','post');
    fore2_backup = padarray(fore2_backup, [960-896 960-909 2^layer_num-mod(zslices, 2^layer_num)], 0, 'post');
end
if reg_ind == 7
    fore2_backup(541:570, 721:750, 193:200) = 0;
    fore2_backup(811:840, 421:450, 177:184) = 0;
end
%%             size          block
% layer 0: 960 960 160       512
% layer 1: 480 480 80        64
% layer 2: 240 240 40        8
% layer 3: 120 120 20        1
% layer 4: 60  60  10      
%%
clc;
clearvars -except data1_backup data2_backup fore2_backup reg_ind data1_name save_folder

% accelerated version of Registration_opticalflow_nonrigid_v2.m
% results should guranteee the same  

sigma_gaussian = 1;
layer_num = 5;   
batch_size = round(size(data1_backup)/2^layer_num);
pad_size = [15 15 5];
step = 1;
smoothness_para = 1;
% resize_method = 'nearest';
resize_method = 'linear';
% resize_method = 'cubic';

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

time_local = 0;
time_local2 = 0;
tic;
for layer = layer_num:-1:0
    
    data1 = imgaussfilt(data1_backup,sigma_gaussian*(layer+1));
    data1 = imresize3(data1, round(size(data1_backup)/2^layer));
    data2 = imgaussfilt(data2_backup,sigma_gaussian*(layer+1));
    data2 = imresize3(data2, round(size(data1_backup)/2^layer));    
    [x,y,z] = size(data1);

    data1_pad = gpuArray(padarray(data1,pad_size,'replicate')); 
    data2_pad = padarray(data2,pad_size,'replicate');
    gt2 = data2;

    fore2 = imresize3(fore2_backup, round(size(data1_backup)/2^layer));  
    fore2 = fore2 >= 0.5;

    
        time_local = time_local - toc;
    if layer == layer_num
        batch_num = 1;
        batch_bin = logical([1 1 1]');
        batch_loc = [1 x 1 y 1 z];
        L = gpuArray(zeros(3,3));
        smoothness = 0;
        fore_ind = ones(sum(fore2(:)), 1);
        fore_ind_list{1} = find(fore_ind==1);
    else
        % build graph relationship
        batch_ind = 0;
        batch_table = zeros(2^(3*(layer_num-layer)), 3);
        batch_loc = zeros(2^(3*(layer_num-layer)), 6);   %[x_s x_e y_s ...]
        batch_bin_ind = 0;
        batch_bin = zeros(2^(3*(layer_num-layer)),1);    % non-zero batch index
        fore_ind = zeros(size(fore2));
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
                    fore_temp = logical(fore2(x_start:x_end, y_start:y_end, z_start:z_end));
                    batch_bin_ind = batch_bin_ind + 1;
                    if any(fore_temp, 'all')
                        if sum(fore_temp(:)) > 30    % enough samples
                            batch_ind = batch_ind + 1;
                            batch_table(batch_ind,:) = [xx yy zz];
                            batch_loc(batch_ind,:) = [x_start x_end y_start y_end z_start z_end];
                            fore_ind(x_start:x_end, y_start:y_end, z_start:z_end) = fore_temp*batch_ind;
                            batch_bin(batch_bin_ind) = 1;
                        else
                            fore2(x_start:x_end, y_start:y_end, z_start:z_end) = false;
                        end
                    end                  
                end
            end
        end
        batch_bin = repmat(batch_bin,1,3)';
        batch_bin = logical(batch_bin(:));
        batch_num = batch_ind;
        batch_table = batch_table(1:batch_num,:);
        batch_loc = batch_loc(1:batch_num,:);
        [~, ~, fore_ind] = find(fore_ind);
        fore_ind_list = cell(batch_num, 1);
        for ii = 1:batch_num
            fore_ind_list{ii} = find(fore_ind == ii);
        end

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

        L = zeros(3*batch_num,3*batch_num);
        for ee = 1:edge_ind*3
            L(batch_relation(ee,1),batch_relation(ee,2)) = -batch_relation(ee,3);
            L(batch_relation(ee,2),batch_relation(ee,1)) = -batch_relation(ee,3);
        end
        for ii = 1:3*batch_num
            L(ii,ii) = -sum(L(ii,:));
        end
    end
    time_local = time_local + toc;

    if layer == layer_num
        phi_current = zeros(x,y,z,3);
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
        
%         phi_current = gpuArray(phi_current);
    end

    fore2_vec = find(fore2);
    [x_fore, y_fore, z_fore] = ind2sub(size(fore2), fore2_vec);
    [y_batch, x_batch, z_batch] = meshgrid(0:2^(layer_num - layer)+1);
    y_batch = y_batch*batch_size(2) + 0.5 - batch_size(2)/2;
    x_batch = x_batch*batch_size(1) + 0.5 - batch_size(1)/2;
    z_batch = z_batch*batch_size(3) + 0.5 - batch_size(3)/2;

    [x_ind,y_ind,z_ind] = ind2sub(size(data1),fore2_vec);
    x_ind = gpuArray(x_ind);
    y_ind = gpuArray(y_ind);
    z_ind = gpuArray(z_ind);

loss = gpuArray(zeros(1000,1));
time = gpuArray(zeros(1000,1));

    phi_previous = phi_current;
    x_bias = phi_previous(fore2_vec);
    y_bias = phi_previous(fore2_vec + x*y*z);
    z_bias = phi_previous(fore2_vec + 2*x*y*z);
for iter = 1:25
%     phi_previous = phi_current;
%     x_bias = phi_previous(fore2_vec);
%     y_bias = phi_previous(fore2_vec + x*y*z);
%     z_bias = phi_previous(fore2_vec + 2*x*y*z);


    x_new = x_ind + x_bias;
    y_new = y_ind + y_bias;
    z_new = z_ind + z_bias;
    data1_tran = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));  
    
    % method1 to calculate gradient
    data1_x_incre = interp3(data1_pad,y_new+pad_size(2),x_new+step+pad_size(1),z_new+pad_size(3));
    data1_x_decre = interp3(data1_pad,y_new+pad_size(2),x_new-step+pad_size(1),z_new+pad_size(3));
    Ix = (data1_x_incre - data1_x_decre)/(2*step);

    data1_y_incre = interp3(data1_pad,y_new+step+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
    data1_y_decre = interp3(data1_pad,y_new-step+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
    Iy = (data1_y_incre - data1_y_decre)/(2*step);

    data1_z_incre = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+step+pad_size(3));
    data1_z_decre = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new-step+pad_size(3));
    Iz = (data1_z_incre - data1_z_decre)/(2*step);

    It = data1_tran-gt2(fore2);

    
    AtA = zeros(batch_num*3,batch_num*3);
    AtIt = zeros(batch_num*3,1);
    time_local2 = time_local2 - toc;
    for ii = 1:batch_num
        batch_Ix = Ix(fore_ind_list{ii});
        batch_Iy = Iy(fore_ind_list{ii});
        batch_Iz = Iz(fore_ind_list{ii});
        batch_It = It(fore_ind_list{ii});

        A_ii = [batch_Ix(:) batch_Iy(:) batch_Iz(:)];
        AtA((ii-1)*3+1:ii*3, (ii-1)*3+1:ii*3) = A_ii'*A_ii;
        AtIt((ii-1)*3+1:ii*3, 1) = A_ii'*batch_It(:);
    end
    time_local2 = time_local2 + toc;
    AtA = gpuArray(AtA);
    AtIt = gpuArray(AtIt);
    phi_gradient = -(AtA + smoothness*L)\(AtIt + smoothness*L*phi_current_vec(batch_bin));
    if gather(any(isnan(phi_gradient)))
        error('Wrong result!');
        error_ind = find(isnan(gather(AtA)));
        [error_x, error_y] = ind2sub(size(AtA), error_ind);
        batch_loc(ceil(error_x/3),:)

        % if error_ind is empty, try:
        a = abs(phi_current_vec(batch_bin));
        error_x = find(a(1:3:end-2) > pad_size(1));
        error_y = find(a(2:3:end-1) > pad_size(2));
        error_z = find(a(3:3:end) > pad_size(3));
        batch_loc(error_y,:)

        % third way
        error_x = [];
        for ii = 1:batch_num
            a = AtA((ii-1)*3+1:ii*3, (ii-1)*3+1:ii*3); 
            if rank(a) < 3
                error_x = [error_x; ii];
            end
        end
        batch_loc(error_x,:)
    end
    if max(abs(phi_gradient)) < 1e-4 || norm(phi_gradient) < 1e-6
        break;
    end
    

    phi_current_vec(batch_bin) = phi_current_vec(batch_bin) + phi_gradient;
%     phi_current = imresize4d(reshape(phi_current_vec, 3, []), [x y z], resize_method);
    if batch_num > 1
    x_bias_temp = padarray(reshape(phi_current_vec(1:3:end-2), [2^(layer_num - layer) 2^(layer_num - layer) 2^(layer_num - layer)]), [1 1 1], 'replicate'); 
    x_bias = interp3(y_batch, x_batch, z_batch, x_bias_temp, y_fore, x_fore, z_fore, resize_method);
    y_bias_temp = padarray(reshape(phi_current_vec(2:3:end-1), [2^(layer_num - layer) 2^(layer_num - layer) 2^(layer_num - layer)]), [1 1 1], 'replicate'); 
    y_bias = interp3(y_batch, x_batch, z_batch, y_bias_temp, y_fore, x_fore, z_fore, resize_method);
    z_bias_temp = padarray(reshape(phi_current_vec(3:3:end), [2^(layer_num - layer) 2^(layer_num - layer) 2^(layer_num - layer)]), [1 1 1], 'replicate'); 
    z_bias = interp3(y_batch, x_batch, z_batch, z_bias_temp, y_fore, x_fore, z_fore, resize_method);
    else
        x_bias(:) = phi_current_vec(1);
        y_bias(:) = phi_current_vec(2);
        z_bias(:) = phi_current_vec(3);
    end

    mse = mean((data1_tran - gt2(fore2)).^2);
    if batch_num == 1
        smooth_error = 0;
    else
        smooth_error = smoothness*phi_current_vec(batch_bin)'*L*phi_current_vec(batch_bin)/sum(fore2(:));
    end
    fprintf('Gradient: %f Max Translation:%f %f %f\n',norm(phi_gradient)/length(phi_gradient), max(abs(x_bias)), max(abs(y_bias)), max(abs(z_bias)));
    fprintf('Iteration %d\n Current error:%f MSE: %f Smooth: %f Time:%f\n',iter, mse+smooth_error, mse, smooth_error, toc);
    loss(iter) = mse+smooth_error;
    time(iter) = toc;
end
    mse_ori = mean((data1(fore2) - gt2(fore2)).^2);
    fprintf('Original error: %f\n', mse_ori);
    loss = gather(loss(1:iter));
    time = gather(time(1:iter));
end
phi_current_vec = gather(phi_current_vec);
% save(['/work/Mengfan/Embryo/Registration/nonrigid_result/' num2str(reg_ind) '_l' num2str(layer_num) '_s'  num2str(smoothness_para)  '_' resize_method '.mat'],'phi_current_vec', 'loss', 'mse_ori' ,'-v7.3');
if any(isnan(phi_current_vec))
    error('Wrong result!');
end
save(fullfile(save_folder, [data1_name '.mat']),'phi_current_vec', 'loss', 'mse_ori' ,'-v7.3');
end

% function k = find_in_boundary(X, loc)
%     % loc: [x_start x_end y_start y_end z_start z_end]
%     if size(X) ~= 3 || length(loc) ~= 6
%         error('Only 3D matrix.');
%     end
%     batch_loc_temp = num2cell(loc);
%     [x_start, x_end, y_start, y_end, z_start, z_end] = deal(batch_loc_temp{:});
%     X_batch = X(x_start:x_end, y_start:y_end, z_start:z_end);
%     k = find(X_batch);
%     [x_ind, y_ind, z_ind] = ind2sub(size(X_batch), k);
%     k = sub2ind(size(X), x_ind+x_start-1, y_ind+y_start-1, z_ind+z_start-1);
% end

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