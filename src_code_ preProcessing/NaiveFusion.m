clc;clear;close all;
% The scirpt implement naive fusion algorithm. That is, using the orignal 
% transforming matrix without any refinment

%% system and path
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
    path_name = '/work/Mengfan/Embryo/21-04-28';
    source_data = 'H2BGFP_21-04-28.h5';
    data_name = fullfile(path_name, source_data);
    xml_name = fullfile(path_name, 'H2BGFP_21-04-28.xml');
    target_folder = 'naiveFusion';
else
    addpath D:\MatlabTools;
    path_name = 'H:\Embryo\TM0-49\';
end
xml = xml2struct(xml_name);
register_info = xml.SpimData.ViewRegistrations.ViewRegistration;
[h5_struct, num_view] = readh5info(fullfile(path_name, source_data));
num_time = length(h5_struct);
num_total = num_time * num_view;
if num_view ~= 8
    error('Current version is for 8 views.');
end
%% parameter setting
device = 'GPU'; % GPU 25s CPU 45s parallel 250s
save_mode = 'tif'; %save as 'tif', 'h5', or 'both' 

downsample_scale = 4;
pad_size = ceil(40 / downsample_scale); % fusion refinement constriant
maxIter = 1000; % max iteration of step 4
baseIntensity = 200; % may not be the best threhsold, refine later

target_folder = [target_folder '_' num2str(1/downsample_scale)];
%% primary fusion
if ~isfolder(fullfile(path_name, target_folder))
    mkdir(fullfile(path_name, target_folder));
end
tic;
for tt = 100:100
    tt_ind = num2str(100000+tt);
    tt_ind = tt_ind(2:6);
    fprintf('processing: %d\n',tt);
    
    % 1. get transformation matrix
    trans = {};
    for vv = 1:4
        trans{vv} = eye(4,4);
        for jj = 1:length(register_info{vv}.ViewTransform)
            trans_temp = cellfun(@str2num, split(register_info{tt*num_view + vv}.ViewTransform{jj}.affine.Text));
            trans_temp = [reshape(trans_temp,4,3)'; 0 0 0 1];
            trans{vv} = trans{vv}*trans_temp;
        end
        trans{vv}(1:3,:) = trans{vv}(1:3,:) / downsample_scale;  % downsample
    end
    
    % 2. get corrdiante range
    ref = {};
    x_min = inf; y_min = inf; z_min = inf;
    x_max = -inf;y_max = -inf;z_max = -inf;
    for vv = 1:4
        vv_ind = num2str(100+vv-1);
        vv_ind = vv_ind(2:3);
        im_size = h5info(data_name, ['/t' tt_ind '/s' vv_ind '/0/cells']);
        im_size = im_size.Dataspace.Size;
        ref{vv} = affineOutputView(im_size, affine3d(trans{vv}'),'BoundsStyle','FollowOutput');

        x_min = min(x_min,round(ref{vv}.XWorldLimits(1)));
        x_max = max(x_max,round(ref{vv}.XWorldLimits(2)));
       
        y_min = min(y_min,round(ref{vv}.YWorldLimits(1)));
        y_max = max(y_max,round(ref{vv}.YWorldLimits(2)));
        
        z_min = min(z_min,round(ref{vv}.ZWorldLimits(1)));
        z_max = max(z_max,round(ref{vv}.ZWorldLimits(2)));
    end
    
    % 3. read and transform data
    x_start = cell(4,1); x_end = cell(4,1);
    y_start = cell(4,1); y_end = cell(4,1);
    z_start = cell(4,1); z_end = cell(4,1);
    im = cell(4,1); im_fuse = cell(4,1);
    for vv = 1:4
        vv_ind = num2str(100+vv-1);
        vv_ind = vv_ind(2:3);
        im{vv} = hdf5read(data_name,['/t' tt_ind '/s' vv_ind '/0/cells']);
        im{vv} = single(im{vv});
        im_size = size(im{vv});

        vv_ind = num2str(100+vv+3);
        vv_ind = vv_ind(2:3);
        im2 = hdf5read(data_name,['/t' tt_ind '/s' vv_ind '/0/cells']);
        im2 = single(im2);
        im{vv} = (im{vv} + im2)/2;
        im{vv} = im{vv} - baseIntensity;
        for zz = 1:size(im{vv},3)
            im{vv}(:,:,zz) = im{vv}(:,:,zz)';
        end

        x_start{vv} = round(ref{vv}.XWorldLimits(1) - x_min + 1);
        x_end{vv} = round(ref{vv}.XWorldLimits(2) - x_min);
        y_start{vv} = round(ref{vv}.YWorldLimits(1) - y_min + 1);
        y_end{vv} = round(ref{vv}.YWorldLimits(2) - y_min);
        z_start{vv} = round(ref{vv}.ZWorldLimits(1) - z_min + 1);
        z_end{vv} = round(ref{vv}.ZWorldLimits(2) - z_min);

        im_fuse{vv} = zeros(y_max-y_min,x_max-x_min,z_max-z_min,'single');
        im_fuse{vv}(y_start{vv}:y_end{vv},x_start{vv}:x_end{vv},z_start{vv}:z_end{vv}) = ...
            imwarp(im{vv},affine3d(trans{vv}'));
        im_fuse{vv} = padarray(im_fuse{vv}, [pad_size pad_size pad_size]);
    end
    if strcmp(device, 'GPU')
        for vv = 1:4
            im_fuse{vv} = gpuArray(im_fuse{vv});
        end
    end

    % 6. final fusion
    x_min = inf; y_min = inf; z_min = inf;
    x_max = -inf;y_max = -inf;z_max = -inf;
    for vv = 1:4 % update transform matrix
        ref{vv} = affineOutputView(size(im{vv}), affine3d(trans{vv}'),'BoundsStyle','FollowOutput');

        x_min = min(x_min,round(ref{vv}.XWorldLimits(1)));
        x_max = max(x_max,round(ref{vv}.XWorldLimits(2)));
       
        y_min = min(y_min,round(ref{vv}.YWorldLimits(1)));
        y_max = max(y_max,round(ref{vv}.YWorldLimits(2)));
        
        z_min = min(z_min,round(ref{vv}.ZWorldLimits(1)));
        z_max = max(z_max,round(ref{vv}.ZWorldLimits(2)));
    end

    im_fusion = zeros(y_max-y_min,x_max-x_min,z_max-z_min,'single');
    im_num = zeros(y_max-y_min,x_max-x_min,z_max-z_min,'single');
    for vv = 1:4
        im_unit = ones(size(im{vv}),'single');
        if strcmp(device, 'GPU')
            im{vv} = gpuArray(im{vv});
            im_unit = gpuArray(im_unit);
        end
        im{vv} = imwarp(im{vv},affine3d(trans{vv}'));
        im_unit = imwarp(im_unit,affine3d(trans{vv}'));
        im{vv} = im{vv} .* im_unit;

        x_start = round(ref{vv}.XWorldLimits(1) - x_min + 1);
        x_end = round(ref{vv}.XWorldLimits(2) - x_min);
        y_start = round(ref{vv}.YWorldLimits(1) - y_min + 1);
        y_end = round(ref{vv}.YWorldLimits(2) - y_min);
        z_start = round(ref{vv}.ZWorldLimits(1) - z_min + 1);
        z_end = round(ref{vv}.ZWorldLimits(2) - z_min);

        if strcmp(device, 'GPU')
            im_fusion(y_start:y_end,x_start:x_end,z_start:z_end) = ...
            im_fusion(y_start:y_end,x_start:x_end,z_start:z_end) + gather(im{vv});
            im_num(y_start:y_end,x_start:x_end,z_start:z_end) = ...
            im_num(y_start:y_end,x_start:x_end,z_start:z_end) + gather(im_unit);
        else
            im_fusion(y_start:y_end,x_start:x_end,z_start:z_end) = ...
            im_fusion(y_start:y_end,x_start:x_end,z_start:z_end) + im{vv};
            im_num(y_start:y_end,x_start:x_end,z_start:z_end) = ...
            im_num(y_start:y_end,x_start:x_end,z_start:z_end) + im_unit;
        end
    end
    im_fusion = im_fusion./im_num;
    im_fusion(isnan(im_fusion)) = 0;
    if strcmp(save_mode, 'h5') || strcmp(save_mode, 'both')
        h5create(fullfile(path_name, target_folder, ['fusion_' source_data]),...
            ['/t' tt_ind '/0/cells'],size(im_fusion),'Datatype','single');
        h5write(fullfile(path_name, target_folder, ['fusion_' source_data]),...
            ['/t' tt_ind '/0/cells'],im_fusion);
    end
    if strcmp(save_mode, 'tif') || strcmp(save_mode, 'both')
        tifwrite(uint16(im_fusion),fullfile(path_name, target_folder, tt_ind));
    end
    clear im im_fusion im_num im_unit
    toc
end