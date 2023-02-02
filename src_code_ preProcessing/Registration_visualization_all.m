%% system and path
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
    data_folder = '/work/Mengfan/Embryo/22-01-11/sameViewFusion_050-149_11';
    fore_folder = '/work/public/sameViewFusion/sameViewDetection_050-149_11/synQuant_refine_res';
end
addpath /home/mengfan/ForExecute/Tools/MatlabTools/export_fig
%% plot method 1: one singe frame

clc;
reg_ind = 90;
para = 'l4_s0.1_cubic';
1
load(['/work/Mengfan/Embryo/Registration/tmp_result/phi_' num2str(reg_ind) '_' num2str(reg_ind+1) '_rigid.mat']);
load(['/work/Mengfan/Embryo/Registration/nonrigid_result/' num2str(reg_ind) '_' para '.mat']);

data1 = tifread(fullfile(data_folder, ['000' num2str(reg_ind) '.tif']));
% data2 = tifread(fullfile(data_folder, '00083.tif'));
fore1 = load(fullfile(fore_folder, ['000' num2str(reg_ind) '.mat']));
fore1 = fore1.refine_res>0;

data1 = data1(:,:,1:160);
data1_backup = imresize3(data1, [960 960 160]);
fore1 = fore1(:,:,1:160);
fore1 = imresize3(fore1, [960 960 160]);
fore2_backup = fore1 > 0.5;

phi_current_vec = reshape(phi_current_vec,3,[]);
phi_current = imresize4d(phi_current_vec, [960 960 160], 'linear');
im = zeros(1920,1920,3,160);

for zz = 1:4:160
    zz
data1 = data1_backup(:,:,zz);
fore1 = fore2_backup(:,:,zz);
if sum(fore1(:))>0
    ux = phi_current(:,:,zz,1) - phi_rigid(1);
    uy = phi_current(:,:,zz,2) - phi_rigid(2);
    uz = phi_current(:,:,zz,3) - phi_rigid(3);
    
    ux(~fore1) = 0;
    uy(~fore1) = 0;
    uz(~fore1) = 0;
    
    im_temp = plot_vec(data, ux, uy);
    im(:,:,:,zz) = im_temp(1:1920,1:1920,:);
else
    im(:,:,1,zz) = imresize(data1, [1920 1920]);
    im(:,:,2,zz) = imresize(data1, [1920 1920]);
    im(:,:,3,zz) = imresize(data1, [1920 1920]);
end
end
im = im(:,:,:,1:4:160);
tifwrite(uint8(im), ['/work/Mengfan/Embryo/Registration/vectorField_img/' num2str(reg_ind) '_' para]);

%% method 2: plot video stream on a z-slice
clc;close all;
clearvars -except data_folder fore_folder
zz = 100;
t_start = 50; t_end = 148;
im = zeros(1920,1920,3,t_end-t_start+1);
load('/work/public/sameViewFusion/sameViewDetection_050-149_11/tform_050-149_11_translation.mat');

for tt = t_start:t_end
    tt
% load(['/work/Mengfan/Embryo/Registration/tmp_result/phi_' num2str(reg_ind) '_' num2str(reg_ind+1) '_rigid.mat']);
phi_rigid = zeros(3,1);
phi_rigid(1) = -tform{t_start - 49}.T(4,2);
phi_rigid(2) = -tform{t_start - 49}.T(4,1);
phi_rigid(3) = -tform{t_start - 49}.T(4,3);
load(['/work/Mengfan/Embryo/Registration/nonrigid_result_l4_s0.1_linear/' num2str(tt) '.mat']);



data1_name = num2str(100000 + tt);
data1_name = data1_name(2:end);
data1 = tifread(fullfile(data_folder, [data1_name '.tif']));
data1 = data1(:,:,zz);
data1 = imresize(data1, [960 960]);

fore1 = load(fullfile(fore_folder, [data1_name '.mat']));
fore1 = fore1.refine_res>0;
fore1 = fore1(:,:,zz);
fore1 = imresize(fore1, [960 960]);

phi_current_vec = reshape(phi_current_vec,3,[]);
phi_current = imresize4d(phi_current_vec, [960 960 160], 'linear');

ux = phi_current(:,:,zz,1);% - phi_rigid(1);
uy = phi_current(:,:,zz,2);% - phi_rigid(2);
uz = phi_current(:,:,zz,3);% - phi_rigid(3);

ux(~fore1) = 0;
uy(~fore1) = 0;
uz(~fore1) = 0;

im_temp = plot_vec(data1, ux, uy);
im(:,:,:,tt-t_start+1) = imresize(im_temp, [1920 1920]);
end
tifwrite(uint8(im), ['/work/Mengfan/Embryo/Registration/vectorField_img/t' num2str(t_start) '-' num2str(t_end) '_z' num2str(zz) '_s0.1_linear']);

function im = plot_vec(data, ux, uy)
    f=8;
    [x_size, y_size, z_size] = size(ux);
    % Enhance the quiver plot visually by downsizing vectors  
    %   -f : downsizing factor
    v = sqrt(ux.^2 + uy.^2);
    color = colormap;
    v = v(1:f:x_size,1:f:y_size);
    color_quiver = ceil(sqrt(max(v,[],3))*255/sqrt(max(v(:))+1e-6)) + 1; % sqrt: change vector color
    x = ux(1:f:x_size,1:f:y_size);
    y = uy(1:f:x_size,1:f:y_size);
    [X,Y]=meshgrid(1:size(x,2),1:size(x,1));
    scale = f;
    h = imshow(data/255);hold on;
    % h = imshow(imresize(max((data1-data2).^2*2/3,[],3),4));hold on;
    % h = imshow(zeros(size(imresize(max(data1,[],3)/2,scale)))); hold on;
    for ii = 1:256
        ind = color_quiver == ii;
        quiver(X(ind)*scale,Y(ind)*scale,-y(ind)*4,-x(ind)*4,0,'Color',color(ii,:));
    end
    im = export_fig('-transparent','-m3');
end

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