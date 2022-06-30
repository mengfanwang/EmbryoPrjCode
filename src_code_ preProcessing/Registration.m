clc;clear;close all;
%% system and path
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
    path_name = '/work/Mengfan/Embryo/TM0-49';
    input_data_name = fullfile(path_name, 'fusion_0.25/fusion_H2BGFP_TM0-49.h5');
    input_res_folder = fullfile(path_name, 'detection_0.25');
    output_data_folder = fullfile(path_name, 'data_reg_0.25');
    output_res_folder = fullfile(path_name, 'res_reg_0.25');
else
    addpath D:\MatlabTools;
    path_name = 'H:\Embryo\TM0-49\';
end
if ~isfolder(output_data_folder)
    mkdir(output_data_folder);
end
if ~isfolder(output_res_folder)
    mkdir(output_res_folder);
end

%% settings
device = 'GPU'; % GPU 25s CPU 45s parallel 250s

optimizer = registration.optimizer.RegularStepGradientDescent;
metric = registration.metric.MeanSquares;
optimizer.MaximumIterations = 200;
optimizer.MaximumStepLength = 0.04;

%% load data
% different frame may have different size, take the minimum
[h5_struct, ~] = readh5info(input_data_name);
num_time = length(h5_struct);
im_size = [inf inf inf];
for tt = 0:num_time - 1
    dim_temp = h5_struct(tt+1).Groups.Datasets.Dims;
    im_size = min(im_size, dim_temp);
end

%% get registration matrix
trans_mat = cell(num_time,1);
tform = cell(num_time,1);
ref = cell(num_time,1);
x_min = inf; y_min = inf; z_min = inf;
x_max = -inf;y_max = -inf;z_max = -inf;

tic;
%1. get corrdiante range
fprintf('Get transformation matrix...');
for tt = 1:num_time
    tt_ind = num2str(99999+tt);
    tt_ind = tt_ind(2:6);
    fprintf('processing: %s\n',tt_ind);
    if tt == 1
        trans_mat{1} = eye(4,4);
        tform{1} = eye(4,4);
    else
        data_moving = hdf5read(input_data_name,['/t' tt_ind '/0/cells']);
        data_moving(data_moving<0) = 0;
        data_fixed = hdf5read(input_data_name,['/t' tt_ind_old '/0/cells']);
        data_fixed(data_fixed<0) = 0;
        tform{tt} = imregtform(data_moving, data_fixed, 'rigid', optimizer, metric);
        trans_mat{tt} = tform{tt}.T*trans_mat{tt-1};
    end
    ref{tt} = affineOutputView(im_size, affine3d(trans_mat{tt}),'BoundsStyle','FollowOutput');
    
    x_min = min(x_min,round(ref{tt}.XWorldLimits(1)));
    x_max = max(x_max,round(ref{tt}.XWorldLimits(2)));
   
    y_min = min(y_min,round(ref{tt}.YWorldLimits(1)));
    y_max = max(y_max,round(ref{tt}.YWorldLimits(2)));
    
    z_min = min(z_min,round(ref{tt}.ZWorldLimits(1)));
    z_max = max(z_max,round(ref{tt}.ZWorldLimits(2)));

    tt_ind_old = tt_ind;
    toc
end
bound.x_min = x_min;
bound.x_max = x_max;
bound.y_min = y_min;
bound.y_max = y_max;
bound.z_min = z_min;
bound.z_max = z_max;
save(fullfile(output_res_folder, 'registration'), 'bound', 'ref', 'tform', 'trans_mat','im_size');
fprintf('Get transform matrix running time:'); % 14453s 0.25
toc
%%
fprintf('\nRegistration...');
load(fullfile(input_res_folder, 'varianceMap.mat'));
load(fullfile(input_res_folder, 'synQuant_priCvt_res.mat'));
load(fullfile(input_res_folder, 'synQuant_refine_res_4d_v9.mat'));
load(fullfile(output_res_folder, 'registration'));
x_max = bound.x_max; x_min = bound.x_min;
y_max = bound.y_max; y_min = bound.y_min;
z_max = bound.z_max; z_min = bound.z_min;

%%
eig_res_2d_reg = eig_res_2d;
eig_res_3d_reg = eig_res_3d;
varMap_reg = varMap;
tic;
for tt = 1:num_time
    tt_ind = num2str(99999+tt);
    tt_ind = tt_ind(2:6);
    fprintf('processing: %s\n',tt_ind);

    x_start = round(ref{tt}.XWorldLimits(1) - x_min + 1);
    x_end = round(ref{tt}.XWorldLimits(2) - x_min);
    y_start = round(ref{tt}.YWorldLimits(1) - y_min + 1);
    y_end = round(ref{tt}.YWorldLimits(2) - y_min);
    z_start = round(ref{tt}.ZWorldLimits(1) - z_min + 1);
    z_end = round(ref{tt}.ZWorldLimits(2) - z_min);

    % crop data to the same size
    data = hdf5read(input_data_name,['/t' tt_ind '/0/cells']);
    data = data(1:im_size(1),1:im_size(2),1:im_size(3));
    eig_res_2d{tt} = eig_res_2d{tt}(1:im_size(1),1:im_size(2),1:im_size(3));
    eig_res_3d{tt} = eig_res_3d{tt}(1:im_size(1),1:im_size(2),1:im_size(3));
    varMap{tt}{1,1} = varMap{tt}{1,1}(1:im_size(1),1:im_size(2),1:im_size(3));
    varMap{tt}{1,2} = varMap{tt}{1,2}(1:im_size(1),1:im_size(2),1:im_size(3));

    if strcmp(device, 'GPU')
        data = gpuArray(data);
    end
    data_temp = imwarp(data,affine3d(trans_mat{tt}));
    if strcmp(device, 'GPU')
        data_temp = gather(data_temp);
    end
    data_reg = zeros(y_max - y_min, x_max - x_min, z_max - z_min, 'single');
    data_reg(y_start:y_end, x_start:x_end, z_start:z_end) = data_temp;
    
    eig_temp = eig_res_2d{tt};
    if strcmp(device, 'GPU')
        eig_temp = gpuArray(eig_temp);
    end
    eig_temp = imwarp(eig_temp,affine3d(trans_mat{tt}));
    if strcmp(device, 'GPU')
        eig_temp = gather(eig_temp);
    end
    eig_res_2d_reg{tt} = zeros(y_max - y_min, x_max - x_min, z_max - z_min, 'single');
    eig_res_2d_reg{tt}(y_start:y_end, x_start:x_end, z_start:z_end) = single(eig_temp);
    
    eig_temp = eig_res_3d{tt};
    if strcmp(device, 'GPU')
        eig_temp = gpuArray(eig_temp);
    end
    eig_temp = imwarp(eig_temp,affine3d(trans_mat{tt}));
    if strcmp(device, 'GPU')
        eig_temp = gather(eig_temp);
    end
    eig_res_3d_reg{tt} = zeros(y_max - y_min, x_max - x_min, z_max - z_min, 'single');
    eig_res_3d_reg{tt}(y_start:y_end, x_start:x_end, z_start:z_end) = single(eig_temp);
    
    var_temp = fillmissing(fillmissing(fillmissing(varMap{tt}{1,1},'nearest',1),'nearest',2),'nearest',3);
    if strcmp(device, 'GPU')
        var_temp = gpuArray(var_temp);
    end
    var_temp = imwarp(var_temp ,affine3d(trans_mat{tt}));
    if strcmp(device, 'GPU')
        var_temp = gather(var_temp);
    end
    varMap_reg{tt}{1,1} = nan(y_max - y_min, x_max - x_min, z_max - z_min);
    varMap_reg{tt}{1,1}(y_start:y_end, x_start:x_end, z_start:z_end) = var_temp;
    
    var_temp = fillmissing(fillmissing(fillmissing(varMap{tt}{1,2},'nearest',1),'nearest',2),'nearest',3);
    if strcmp(device, 'GPU')
        var_temp = gpuArray(var_temp);
    end
    var_temp = imwarp(var_temp ,affine3d(trans_mat{tt}));
    if strcmp(device, 'GPU')
        var_temp = gather(var_temp);
    end
    varMap_reg{tt}{1,2} = nan(y_max - y_min, x_max - x_min, z_max - z_min);
    varMap_reg{tt}{1,2}(y_start:y_end, x_start:x_end, z_start:z_end) = var_temp;
    
    tifwrite(uint16(data_reg), fullfile(output_data_folder, tt_ind));
    toc
end
save(fullfile(output_res_folder, 'synQuant_pri_var_reg.mat'), 'eig_res_2d_reg', 'eig_res_3d_reg', 'varMap_reg', '-v7.3');
toc
