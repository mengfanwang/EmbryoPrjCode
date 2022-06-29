clc;clear;close all;

data_folder = 'H:\Embryo\TM0-49\multiview_crop\';

% load data
data = zeros(200,350,50,40);
for tt = 1:40
    ind = num2str(99+tt);
    ind = ind(2:end);
    load([data_folder ind '.mat']);
    data(:,:,:,tt) = im_multi;
end

optimizer = registration.optimizer.RegularStepGradientDescent;
metric = registration.metric.MeanSquares;
optimizer.MaximumIterations = 200;
optimizer.MaximumStepLength = 0.04;

t = 40;
trans_mat = cell(t,1);
tform = cell(t,1);
trans_mat{1} = eye(4,4);
tform{1} = eye(4,4);
% registration
tic;
fprintf('Get transformation matrix...');
for tt = 2:t
    fprintf(' %d', tt);
    % tform = imregtform(data(:,:,:,6), data(:,:,:,5), 'rigid', optimizer, metric,'DisplayOptimization',true,'PyramidLevels',1);
    tform{tt} = imregtform(data(:,:,:,tt), data(:,:,:,tt-1), 'rigid', optimizer, metric);
    trans_mat{tt} = tform{tt}.T*trans_mat{tt-1};
end
fprintf('\nGet coordinate...');
ref = cell(t,1);
x_min = inf; y_min = inf; z_min = inf;
x_max = -inf;y_max = -inf;z_max = -inf;
for tt = 1:t
    fprintf(' %d', tt);
    [~, ref{tt}] = imwarp(data(:,:,:,tt),affine3d(trans_mat{tt}));
    
    x_min = min(x_min,round(ref{tt}.XWorldLimits(1)));
    x_max = max(x_max,round(ref{tt}.XWorldLimits(2)));
   
    y_min = min(y_min,round(ref{tt}.YWorldLimits(1)));
    y_max = max(y_max,round(ref{tt}.YWorldLimits(2)));
    
    z_min = min(z_min,round(ref{tt}.ZWorldLimits(1)));
    z_max = max(z_max,round(ref{tt}.ZWorldLimits(2)));
end
%%
mkdir('..\registration_temporal_data\ImageReg');
fprintf('\nRegistration...');
data_reg = zeros(y_max-y_min,x_max-x_min,z_max-z_min,t);
for tt = 1:t
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
bound.x_min = x_min;
bound.x_max = x_max;
bound.y_min = y_min;
bound.y_max = y_max;
bound.z_min = z_min;
bound.z_max = z_max;
save('I:\Embryo\TM0-49\foreground_crop\yinan_data_mengfan_fuse_res\ImageReg','trans_mat','tform','bound','ref','data_reg','-v7.3');


