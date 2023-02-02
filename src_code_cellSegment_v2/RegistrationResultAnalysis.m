clc;clear;close all

% ground truth
% 1874: [1 100]
load('/work/Mengfan/Embryo/Registration/ground_truth/1823.mat');
gt_x = xCoord(trajectory);   
gt_y = yCoord(trajectory);   
gt_z = zCoord(trajectory);   

%% 1000 tracks loading
tic;load('/work/Mengfan/Embryo/Registration/ground_truth/track_100.mat');toc

%% 1000 tracks prediction
clearvars -except target_track tracks vox voxIdx xCoord yCoord zCoord

ts = 1; te = 100;
% gt_loc = zeros(100,3,length(target_track));
figure(1);
for ii = 1:length(target_track)
    trajectory = tracks{target_track(ii)};
    gt_x = xCoord(trajectory);   
    gt_y = yCoord(trajectory);   
    gt_z = zCoord(trajectory);   
    plot3(gt_x(ts:te), gt_y(ts:te), gt_z(ts:te), [0 0.4470 0.7410], 'LineWidth',1); hold on;
end

%% plot prediction (backforward and should correct)
ts = 1;   
te = 30;
pred_x = zeros(100,1);
pred_y = zeros(100,1);
pred_z = zeros(100,1);
pred_x(te) = gt_y(te);
pred_y(te) = gt_x(te);
pred_z(te) = gt_z(te);
for ii = te:-1:ts+1
    ii
    load(['/work/Mengfan/Embryo/Registration/nonrigid_result_l4_s1_linear/', num2str(ii+48), '.mat']);
    [x_bias, y_bias, z_bias] = moving_predict(phi_current_vec, pred_x(ii), pred_y(ii), pred_z(ii),4);
    pred_x(ii-1) = pred_x(ii) + x_bias;
    pred_y(ii-1) = pred_y(ii) + y_bias;
    pred_z(ii-1) = pred_z(ii) + z_bias;
end


figure(2);
plot3(gt_x(ts:te), gt_y(ts:te), gt_z(ts:te),'LineWidth',2); hold on;
plot3(pred_y(ts:te), pred_x(ts:te), pred_z(ts:te), 'LineWidth',2);
% axis([480 590 370 470 130 150]);

%% plot error of each prediction (backforward and should correct)
load('/work/public/sameViewFusion/sameViewDetection_050-149_11/tform_050-149_11_translation.mat');
xx = zeros(99,1);for ii = 1:99;xx(ii) = tform{ii}.T(4,2);end
yy = zeros(99,1);for ii = 1:99;yy(ii) = tform{ii}.T(4,1);end
zz = zeros(99,1);for ii = 1:99;zz(ii) = tform{ii}.T(4,3);end

ts = 1;         %time start
te = 100;
distance_new = zeros(te-1,1);
distance_old = zeros(te-1,1);
distance_rigid = zeros(te-1,1);
x_bias = zeros(te,1);
y_bias = zeros(te,1);
z_bias = zeros(te,1);
for ii = te-1:-1:ts
    ii
    load(['/work/Mengfan/Embryo/Registration/nonrigid_result_l4_s1_linear/', num2str(ii+49), '.mat']);
    [x_bias(ii), y_bias(ii), z_bias(ii)] = moving_predict(phi_current_vec, gt_y(ii+1), gt_x(ii+1), gt_z(ii+1),4);
    bias = [gt_y(ii+1)-gt_y(ii)+x_bias(ii) gt_x(ii+1)-gt_x(ii)+y_bias(ii) gt_z(ii+1)-gt_z(ii)+z_bias(ii)];
    distance_new(ii) = norm(bias);
    distance_old(ii) = norm([gt_y(ii+1)-gt_y(ii) gt_x(ii+1)-gt_x(ii) gt_z(ii+1)-gt_z(ii)]);
    
    bias = [gt_y(ii+1)-gt_y(ii)-xx(ii) gt_x(ii+1)-gt_x(ii)-yy(ii) gt_z(ii+1)-gt_z(ii)-zz(ii)];
    distance_rigid(ii) = norm(bias);
end
% plot(distance_old);hold on;plot(distance_new);
figure(1);plot(distance_old - distance_new);
figure(2);plot(distance_old - distance_rigid);


%% align patch to a image
ii = 88;
load(['/work/Mengfan/Embryo/Registration/nonrigid_result_l4_s1_linear/', num2str(ii+49), '.mat']);
target_loc = [406 634 122];
patch_loc = ceil(target_loc./[60 60 10]);
patch_loc(2) = patch_loc(2)-1;
vec_ind = sub2ind([16 16 16],patch_loc(1), patch_loc(2), patch_loc(3));

x_start = (patch_loc(1)-1)*batch_size(1) + 1;
x_end = patch_loc(1)*batch_size(1);
y_start = (patch_loc(2)-1)*batch_size(2) + 1;
y_end = patch_loc(2)*batch_size(2);
z_start = (patch_loc(3)-1)*batch_size(3) + 1;
z_end = patch_loc(3)*batch_size(3);

data_folder = '/work/Mengfan/Embryo/22-01-11/sameViewFusion_050-149_11';
data1_name = num2str(100049 + ii);
data1_name = data1_name(2:end);
data1 = tifread(fullfile(data_folder, [data1_name '.tif']));
data1 = data1(:,:,1:160);
data1 = imresize3(data1, [960 960 160]);

phi_current_vec = reshape(phi_current_vec, 3, []);
phi_current = imresize4d(phi_current_vec, [960 960 160], 'nearest');
[x_ind,y_ind,z_ind] = ind2sub(size(data1),find(ones(size(data1))));
phi_previous = phi_current;
x_bias = phi_previous(:,:,:,1);
y_bias = phi_previous(:,:,:,2);
z_bias = phi_previous(:,:,:,3);

data2_name = num2str(100049 + ii + 1);
data2_name = data2_name(2:end);
data2 = tifread(fullfile(data_folder, [data2_name '.tif']));
data2 = data2(:,:,1:160);
data2 = imresize3(data2, [960 960 160]);

im = zeros(960,960,3);
im(:,:,1) = data2(:,:,target_loc(3))/2;
im(:,:,2) = data2(:,:,target_loc(3))/2;
im(:,:,3) = data2(:,:,target_loc(3))/2;
im(target_loc(1), target_loc(2), 2) = 255;
im3 = data2(:,:,target_loc(3));
pad_size = [15 15 5];
data1_pad = padarray(data1,pad_size,'replicate'); 
mse = mean((data1(x_start:x_end,y_start:y_end, z_start:z_end) - ...
                data2(x_start:x_end,y_start:y_end, z_start:z_end)).^2 , 'all')
% local patch transform
x_new = x_ind + x_bias(:);
y_new = y_ind + y_bias(:);
z_new = z_ind + z_bias(:);
data1_tran = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
data1_tran = reshape(data1_tran, [960 960 160]);
im2 = im;
im2(x_start:x_end,y_start:y_end,1) = im2(x_start:x_end,y_start:y_end,1) + data1_tran(x_start:x_end,y_start:y_end,target_loc(3));
im4 = im3;
im4(x_start:x_end,y_start:y_end,1) = abs(im3(x_start:x_end,y_start:y_end,1) - data1_tran(x_start:x_end,y_start:y_end,target_loc(3)));
figure(1);imshow(im2/255);
figure(3);imshow(im4/255);
mse = mean((data1_tran(x_start:x_end,y_start:y_end, z_start:z_end) - ...
                data2(x_start:x_end,y_start:y_end, z_start:z_end)).^2 , 'all')

% rigid transform
load('/work/public/sameViewFusion/sameViewDetection_050-149_11/tform_050-149_11_translation.mat');
x_bias = -tform{ii}.T(4,2);
y_bias = -tform{ii}.T(4,1);
z_bias = -tform{ii}.T(4,3);
x_new = x_ind + x_bias(:);
y_new = y_ind + y_bias(:);
z_new = z_ind + z_bias(:);
data1_tran = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
data1_tran = reshape(data1_tran, [960 960 160]);
im2 = im;
im2(x_start:x_end,y_start:y_end,1) = im2(x_start:x_end,y_start:y_end,1) + data1_tran(x_start:x_end,y_start:y_end,target_loc(3));
im4 = im3;
im4(x_start:x_end,y_start:y_end,1) = abs(im3(x_start:x_end,y_start:y_end,1) - data1_tran(x_start:x_end,y_start:y_end,target_loc(3)));
figure(2);imshow(im2/255);
figure(4);imshow(im4/255);
mse = mean((data1_tran(x_start:x_end,y_start:y_end, z_start:z_end) - ...
                data2(x_start:x_end,y_start:y_end, z_start:z_end)).^2 , 'all')



%% plot accumulate translation of foreground of norigid registration
fore_folder = '/work/public/sameViewFusion/sameViewDetection_050-149_11/synQuant_refine_res';
x_bias = zeros(99,1);
y_bias = zeros(99,1);
z_bias = zeros(99,1);
for ii = 1:99
    ii
    load(['/work/Mengfan/Embryo/Registration/nonrigid_result_l4_s0.1_linear/', num2str(ii+49), '.mat']);
    fore_name = num2str(100049+ii);
    fore_name = fore_name(2:end);
    fore = load(fullfile(fore_folder, [fore_name '.mat']),'refine_res');
    fore = fore.refine_res>0;
    fore = fore(:,:,1:160);
    fore = imresize3(fore, [960 960 160]);
    fore = fore > 0.5;
    
    phi_current_vec = reshape(phi_current_vec, 3, []);
    phi_current = imresize4d(phi_current_vec, [960 960 160], 'linear');
    phi_x = phi_current(:,:,:,1);
    x_bias(ii) = mean(phi_x(fore));
    phi_y = phi_current(:,:,:,2);
    y_bias(ii) = mean(phi_y(fore));
    phi_z = phi_current(:,:,:,3);
    z_bias(ii) = mean(phi_z(fore));
end

%% plot accumulate translation of rigid registration
close all
load('/work/public/sameViewFusion/sameViewDetection_050-149_11/tform_050-149_11_translation.mat');
xx = zeros(99,1);for ii = 1:99;xx(ii) = tform{ii}.T(4,2);end
yy = zeros(99,1);for ii = 1:99;yy(ii) = tform{ii}.T(4,1);end
zz = zeros(99,1);for ii = 1:99;zz(ii) = tform{ii}.T(4,3);end

close all
subplot(3,1,1); plot(cumsum(xx)); hold on;
subplot(3,1,2); plot(cumsum(yy)); hold on;
subplot(3,1,3); plot(cumsum(zz)); hold on;

for ii = 1:99
    load(['/work/Mengfan/Embryo/Registration/nonrigid_result_l4_s0.1_linear/', num2str(ii+49), '.mat']);
    xx(ii) = -mean(phi_current_vec(1:3:end-2));
    yy(ii) = -mean(phi_current_vec(2:3:end-1));
    zz(ii) = -mean(phi_current_vec(3:3:end));
end

subplot(3,1,1); plot(cumsum(xx)); hold on;
subplot(3,1,2); plot(cumsum(yy)); hold on;
subplot(3,1,3); plot(cumsum(zz)); hold on;



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

function [x_bias, y_bias, z_bias] = moving_predict(phi_current_vec, x_start, y_start, z_start, layer)
    [y_batch, x_batch, z_batch] = meshgrid(0:2^layer+1);
    grid_size = [2^layer 2^layer 2^layer];
    batch_size = round([960 960 160]./grid_size);
    resize_method = 'linear';

    y_batch = y_batch*batch_size(2) + 0.5 - batch_size(2)/2;
    x_batch = x_batch*batch_size(1) + 0.5 - batch_size(1)/2;
    z_batch = z_batch*batch_size(3) + 0.5 - batch_size(3)/2;
    x_pad = padarray(reshape(phi_current_vec(1:3:end-2), grid_size), [1 1 1], 'replicate'); 
    x_bias = interp3(y_batch, x_batch, z_batch, x_pad, y_start, x_start, z_start, resize_method);
    y_pad = padarray(reshape(phi_current_vec(2:3:end-1), grid_size), [1 1 1], 'replicate'); 
    y_bias = interp3(y_batch, x_batch, z_batch, y_pad, y_start, x_start, z_start, resize_method);
    z_pad = padarray(reshape(phi_current_vec(3:3:end), grid_size), [1 1 1], 'replicate'); 
    z_bias = interp3(y_batch, x_batch, z_batch, z_pad, y_start, x_start, z_start, resize_method);
end