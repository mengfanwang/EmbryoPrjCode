%% system and path
clc;clear;close all;
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
    data_folder = '/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/view12';
end
tif_files = dir(fullfile(data_folder, '/*.tif'));
ds_scale = 2;

%% global registration
optimizer = registration.optimizer.RegularStepGradientDescent;
metric = registration.metric.MeanSquares;
optimizer.MaximumIterations = 500;
optimizer.MaximumStepLength = 0.01;

tform = cell(1,numel(tif_files)-1);

%%
tic;
for ii = 1:51 %52:numel(tif_files)-1
    fprintf('processing %d/%d file\n', ii, numel(tif_files));

    im_a = tifread(fullfile(tif_files(ii).folder, tif_files(ii).name));
    [h, w, slices] = size(im_a);
%     im_a = imresize3(im_a,round([h/ds_scale w/ds_scale slices*ds_scale]));

    im_b = tifread(fullfile(tif_files(ii+1).folder, tif_files(ii+1).name));
    [h, w, slices] = size(im_b);
%     im_b = imresize3(im_b,round([h/ds_scale w/ds_scale slices*ds_scale]));

    tform_temp = imregtform(im_a, im_b, 'rigid', optimizer, metric);
    tform{ii} = tform_temp;
    toc
end
toc

%% save downsampled data
if ~isfolder([data_folder '_ds'])
    mkdir([data_folder '_ds']);
end
for ii = 1:numel(tif_files)
    fprintf('processing %d/%d file\n', ii, numel(tif_files));

    data = tifread(fullfile(tif_files(ii).folder, tif_files(ii).name));
    [h, w, slices] = size(data);
    data = imresize3(data,round([h/ds_scale w/ds_scale slices]));
    [~, org_name, ~] = fileparts(tif_files(ii).name);
    tifwrite(uint8(data),fullfile([data_folder '_ds'], org_name));
end

%% get registration matrix
optimizer = registration.optimizer.RegularStepGradientDescent;
metric = registration.metric.MeanSquares;
optimizer.MaximumIterations = 500;
optimizer.MaximumStepLength = 0.01;

tform = cell(1,numel(tif_files)-1);
tic;
for ii = 1:numel(tif_files)-1
    fprintf('processing %d/%d file\n', ii, numel(tif_files));

    im_a = tifread(fullfile(tif_files(ii).folder, tif_files(ii).name));
    [h, w, slices] = size(im_a);
    im_a = imresize3(im_a,round([h/ds_scale w/ds_scale slices]));

    im_b = tifread(fullfile(tif_files(ii+1).folder, tif_files(ii+1).name));
    [h, w, slices] = size(im_b);
    im_b = imresize3(im_b,round([h/ds_scale w/ds_scale slices]));

    tform_temp = imregtform(im_a, im_b, 'rigid', optimizer, metric);
    tform{ii} = tform_temp;
    toc
end
toc

%% register data
device = 'GPU'; 
load('/work/Mengfan/Embryo/22-01-11/tform_050-149_11_translation.mat');

x_min = inf; y_min = inf; z_min = inf;
x_max = -inf;y_max = -inf;z_max = -inf;
trans_mat = cell(1,numel(tif_files));
ref = cell(1,numel(tif_files));
[h, w, slices]  = size(tifread(fullfile(tif_files(1).folder, tif_files(1).name)));
im_size = [h/ds_scale w/ds_scale slices];

if ~isfolder([data_folder '_reg'])
    mkdir([data_folder '_reg']);
end

tic;
%1. get corrdiante range
fprintf('Get transformation matrix...');
for tt = numel(tif_files):-1:1
    if tt == numel(tif_files)
        trans_mat{numel(tif_files)} = eye(4,4);
    else
        trans_mat{tt} = tform{tt}.T*trans_mat{tt+1};
    end
    ref{tt} = affineOutputView(im_size, affine3d(trans_mat{tt}),'BoundsStyle','FollowOutput');
    
    x_min = min(x_min,round(ref{tt}.XWorldLimits(1)));
    x_max = max(x_max,round(ref{tt}.XWorldLimits(2)));
   
    y_min = min(y_min,round(ref{tt}.YWorldLimits(1)));
    y_max = max(y_max,round(ref{tt}.YWorldLimits(2)));
    
    z_min = min(z_min,round(ref{tt}.ZWorldLimits(1)));
    z_max = max(z_max,round(ref{tt}.ZWorldLimits(2)));
end
bound.x_min = x_min; 
bound.x_max = x_max;
bound.y_min = y_min;
bound.y_max = y_max;
bound.z_min = z_min;
bound.z_max = z_max;

for tt = 1:numel(tif_files)
    fprintf('processing %d/%d file\n', tt, numel(tif_files));

    x_start = round(ref{tt}.XWorldLimits(1) - x_min + 1);
    x_end = round(ref{tt}.XWorldLimits(2) - x_min);
    y_start = round(ref{tt}.YWorldLimits(1) - y_min + 1);
    y_end = round(ref{tt}.YWorldLimits(2) - y_min);
    z_start = round(ref{tt}.ZWorldLimits(1) - z_min + 1);
    z_end = round(ref{tt}.ZWorldLimits(2) - z_min);

    data = tifread(fullfile(tif_files(tt).folder, tif_files(tt).name));
    [h, w, slices] = size(data);
    data = imresize3(data,round([h/ds_scale w/ds_scale slices]));
    [~, org_name, ~] = fileparts(tif_files(tt).name);

    if strcmp(device, 'GPU')
        data = gpuArray(data);
    end
    data_temp = imwarp(data,affine3d(trans_mat{tt}));
    if strcmp(device, 'GPU')
        data_temp = gather(data_temp);
    end
    data_reg = zeros(y_max - y_min, x_max - x_min, z_max - z_min, 'single');
    data_reg(y_start:y_end, x_start:x_end, z_start:z_end) = data_temp;
    tifwrite(uint8(data_reg),fullfile([data_folder '_reg'], org_name));
end
toc