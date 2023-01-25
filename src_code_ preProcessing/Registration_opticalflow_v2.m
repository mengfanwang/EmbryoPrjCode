%% system and path
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
    data_folder = '/work/Mengfan/Embryo/22-01-11/sameViewFusion_050-149_11';
    fore_folder = '/work/public/sameViewFusion/sameViewDetection_050-149_11/synQuant_refine_res';
end

data1 = tifread(fullfile(data_folder, '00082.tif'));
data2 = tifread(fullfile(data_folder, '00083.tif'));
% fore1 = load(fullfile(fore_folder, '00082.mat'));
fore2 = load(fullfile(fore_folder, '00083.mat'));
fore2 = fore2.refine_res>0;

data1 = data1(:,:,1:160);
data1_backup = imresize3(data1, [960 960 160]);
data2 = data2(:,:,1:160);
data2_backup = imresize3(data2, [960 960 160]);

% fore2_backup = ones(size(data2_backup));
%%                           block
% layer 1: 960 960 160       512
% layer 2: 480 480 80        64
% layer 3: 240 240 40        8
% layer 4: 120 120 20        1
% layer 5: 60  60  10
%%
clc;close all;
clearvars -except data1_backup data2_backup fore2_backup

sigma_gaussian = 1;
layer_num = 2;
tform_rigid = zeros(3,3);
loss_rigid = zeros(1,3);
tic;
for layer = layer_num:-1:0
    
    data1 = imresize3(data1_backup, round([960 960 160]/2^layer));
    data2 = imresize3(data2_backup, round([960 960 160]/2^layer));

%     data1 = imgaussfilt(data1_backup,sigma_gaussian*(layer_num+1));
%     data1 = imresize3(data1, round([960 960 160]/2^layer));
%     data2 = imgaussfilt(data2_backup,sigma_gaussian*(layer_num+1));
%     data2 = imresize3(data2, round([960 960 160]/2^layer));  
    fore2 = imresize3(fore2_backup, round([960 960 160]/2^layer));  
    fore2 = fore2 >= 0.5;

[x,y,z,t] = size(data1);
% lambda = 0.02;
pad_size = [10 10 6];
step = 1;
lr = 0.0001;
decay = 0.9;

% block_size = [120 120 20]; % for non-rigid case

data1_pad = padarray(data1,pad_size,'replicate'); 
data2_pad = padarray(data2,pad_size,'replicate');

% gt2 = imresize3(data2_backup, round([960 960 160]/2^layer));
gt2 = data2;
if layer == layer_num
    phi_current = gpuArray(zeros(x,y,z,3));
else 
    phi_current_temp = phi_current*2;
    phi_current = gpuArray(zeros(x,y,z,3));
    phi_current(:,:,:,1) = unique(phi_current_temp(:,:,:,1));
    phi_current(:,:,:,2) = unique(phi_current_temp(:,:,:,2));
    phi_current(:,:,:,3) = unique(phi_current_temp(:,:,:,3));
end

[x_ind,y_ind,z_ind] = ind2sub(size(data1),find(fore2));

loss = gpuArray(zeros(100000,1));
time = gpuArray(zeros(100000,1));

for iter = 1:30
    phi_previous = phi_current;
%     phi_gradient = gpuArray(zeros(x,y,z,3));
%     x_bias = reshape(phi_previous(:,:,:,1),[1 x*y*z]);
%     y_bias = reshape(phi_previous(:,:,:,2),[1 x*y*z]);
%     z_bias = reshape(phi_previous(:,:,:,3),[1 x*y*z]);

    x_bias = phi_previous(find(fore2));
    y_bias = phi_previous(find(fore2) + numel(fore2));
    z_bias = phi_previous(find(fore2) + numel(fore2)*2);

    x_new = x_ind + x_bias;
    y_new = y_ind + y_bias;
    z_new = z_ind + z_bias;
    data1_tran = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));

    % method 2
    data1_x_incre = interp3(data1_pad,y_new+pad_size(2),x_new+1+pad_size(1),z_new+pad_size(3));
    data1_x_decre = interp3(data1_pad,y_new+pad_size(2),x_new-1+pad_size(1),z_new+pad_size(3));
    Ix = (data1_x_incre - data1_x_decre)/2;
    clear data1_x_incre data1_x_decre

    data1_y_incre = interp3(data1_pad,y_new+1+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
    data1_y_decre = interp3(data1_pad,y_new-1+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
    Iy = (data1_y_incre - data1_y_decre)/2;
    clear data1_y_incre data1_y_decre  

    data1_z_incre = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+1+pad_size(3));
    data1_z_decre = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new-1+pad_size(3));
    Iz = (data1_z_incre - data1_z_decre)/2;
    clear data1_z_incre data1_z_decre      

%     % method 3
%     addpath('../src_code_cellSegment');
%     grad_x = gradient3(data1_pad, 'x');
%     Ix = interp3(grad_x,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
%     grad_y = gradient3(data1_pad, 'y');
%     Iy = interp3(grad_y,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
%     grad_z = gradient3(data1_pad, 'z');
%     Iz = interp3(grad_z,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
    

%     Ix = (data1_x_incre - data1_tran)/step;
%     Iy = (data1_y_incre - data1_tran)/step;
%     Iz = (data1_z_incre - data1_tran)/step;
%     It_norm = sum(data1_tran)/sum(gt2(:))*1.02;
    It = data1_tran-gt2(fore2);

    H = [sum(Ix.^2) sum(Ix.*Iy) sum(Ix.*Iz); 
         sum(Ix.*Iy) sum(Iy.^2) sum(Iy.*Iz);
         sum(Ix.*Iz) sum(Iy.*Iz) sum(Iz.^2)];
    b = -[sum(Ix.*It); sum(Iy.*It); sum(Iz.*It)];
    phi_gradient = H\b;
    if norm(phi_gradient) < 1e-6
        break;
    end
    
    fprintf('Gradient: %f %f %f\n',phi_gradient(1),phi_gradient(2),phi_gradient(3));

    phi_current(:,:,:,1) = phi_current(:,:,:,1) + phi_gradient(1);
    phi_current(:,:,:,2) = phi_current(:,:,:,2) + phi_gradient(2);
    phi_current(:,:,:,3) = phi_current(:,:,:,3) + phi_gradient(3);

%     phi_current(phi_current>5) = 5;
%     phi_current(phi_current<-5) = -5;
    mse = mean((data1_tran(:) - gt2(fore2)).^2);
    fprintf('Iteration %d\n Current error:%f Time:%f\n',iter, mse, toc);
    loss(iter) = mse;
    time(iter) = toc;
end
%     data1 = imresize3(data1_backup, round([960 960 160]/2^layer));
%     data1_pad = padarray(data1,pad_size,'replicate'); 
%     data1_tran = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
%     mse = mean((data1_tran(:) - gt2(:)).^2);
%     tform_rigid(3-layer,:) = [unique(phi_current(:,:,:,1)) unique(phi_current(:,:,:,2)) unique(phi_current(:,:,:,3))];
%     loss_rigid(3-layer) = mse;
end
loss = gather(loss(1:iter));
time = gather(time(1:iter));
phi_current = gather(phi_current);
ux = phi_current(:,:,:,1);
uy = phi_current(:,:,:,2);
uz = phi_current(:,:,:,3);
unique(ux)
unique(uy)
unique(uz)
% save(['/work/Mengfan/Embryo/Registration/mse/sigma_', num2str(sigma_gaussian) ,'.mat'], 'tform_rigid', 'loss_rigid');

%% check Taylor expansion
% data1: orignal data; data1_tran: new data
% a = data1(:) + Ix*phi_gradient(1) + Iy*phi_gradient(2) + Iz*phi_gradient(3);
% close all;
% figure(1);
a = (-It.*Ix)/(sum(Ix.^2));
a = reshape(a, [x y z]);
% a(a>200) = 200;
% a(a<-200) = -200;
imagesc(a(:,:,100));

figure(2);
mse1 = zeros(101,1);
mse2 = zeros(101,1);
for xx = 1:101
    xx
x_new = x_ind + (xx-51)/20;
y_new = y_ind;
z_new = z_ind;
data1_tran = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
mse1(xx) = mean((data1_tran(:) - gt2(fore2)).^2);
mse2(xx) = mean((data1(:) + Ix*((xx-51)/20) - gt2(fore2)).^2);
end
plot(-2.5:0.05:2.5,mse1,'LineWidth',2);hold on;
plot(-2.5:0.05:2.5,mse2,'LineWidth',2);
legend('MSE','1st-order expansion');
%% print grount truth
load('/work/public/sameViewFusion/sameViewDetection_050-149_11/tform_050-149_11_translation.mat');
data1 = data1_backup;
data1_pad = padarray(data1,pad_size,'replicate'); 
data1_tran = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
mse = mean((data1_tran(:) - gt2(:)).^2);
fprintf('Our results: %f', mse);

x_bias = tform{33}.T(4,2);
y_bias = tform{33}.T(4,1);
z_bias = tform{33}.T(4,3);
x_new = x_ind - x_bias;
y_new = y_ind - y_bias;
z_new = z_ind - z_bias;
data1_tran = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
mse = mean((data1_tran(:) - gt2(:)).^2);
fprintf('Ground truth: %f', mse);
mse = mean((data1(:) - gt2(:)).^2);
fprintf('Original error: %f\n', mse);