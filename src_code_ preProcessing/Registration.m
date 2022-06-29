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
for tt = 0:num_time-1
    tt_ind = num2str(100000+tt);
    tt_ind = tt_ind(2:6);
    fprintf('processing: %d\n',tt);
    if tt == 0
        trans_mat{1} = eye(4,4);
        tform{1} = eye(4,4);
    else
        data_moving = hdf5read(input_data_name,['/t' tt_ind '/0/cells']);
        data_fixed = hdf5read(input_data_name,['/t' tt_ind_old '/0/cells']);
        if strcmp(device, 'GPU')
            data_moving = gpuArray(data_moving);
            data_fixed = gpuArray(data_fixed);
        end
        tform{tt} = imregtform(data_moving, data_fixed, 'rigid', optimizer, metric);
        trans_mat{tt} = tform{tt}.T*trans_mat{tt-1};
    end
    ref{tt} = affineOutputView(im_size, affine3d(trans_mat{tt}'),'BoundsStyle','FollowOutput');
    
    x_min = min(x_min,round(ref{tt}.XWorldLimits(1)));
    x_max = max(x_max,round(ref{tt}.XWorldLimits(2)));
   
    y_min = min(y_min,round(ref{tt}.YWorldLimits(1)));
    y_max = max(y_max,round(ref{tt}.YWorldLimits(2)));
    
    z_min = min(z_min,round(ref{tt}.ZWorldLimits(1)));
    z_max = max(z_max,round(ref{tt}.ZWorldLimits(2)));

    tt_ind_old = tt_ind;
end
bound.x_min = x_min;
bound.x_max = x_max;
bound.y_min = y_min;
bound.y_max = y_max;
bound.z_min = z_min;
bound.z_max = z_max;
save(fullfile(output_res_folder, 'registration'), 'bound', 'ref', 'tform', 'trans_mat');
fprintf('Get transform matrix running time:'); % 
toc
%%
fprintf('\nRegistration...');
data_reg = zeros(y_max-y_min,x_max-x_min,z_max-z_min,num_time);
for tt = 1:num_time
    fprintf(' %d', tt);
    data_temp = imwarp(data(:,:,:,tt),affine3d(trans_mat{tt}));
    x_start = round(ref{tt}.XWorldLimits(1) - x_min + 1);
    x_end = round(ref{tt}.XWorldLimits(2) - x_min);
    y_start = round(ref{tt}.YWorldLimits(1) - y_min + 1);
    y_end = round(ref{tt}.YWorldLimits(2) - y_min);
    z_start = round(ref{tt}.ZWorldLimits(1) - z_min + 1);
    z_end = round(ref{tt}.ZWorldLimits(2) - z_min);
    data_reg(y_start:y_end,x_start:x_end,z_start:z_end,tt) = data_temp;
    
    ind = num2str(99+tt);
    ind = ind(2:end);
    tifwrite(uint16(data_reg(:,:,:,tt)), ['..\registration_temporal_data\ImageReg\' ind]);
end
toc;

% save data
save('I:\Embryo\TM0-49\foreground_crop\yinan_data_mengfan_fuse_res\ImageReg','trans_mat','tform','bound','ref','data_reg','-v7.3');


