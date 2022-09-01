clc;clear;close all;
% The fusion algothm contain the following steps:
% 1. get the transformation matrix from the xml file.
% 2. get the corrdiante range of all views
% 3. read and transform data
% 4. refine transform matrix, minimizing mean squared error
% 5. refine result again to make sure result consistent
% 6. remove boundary, register data and save

%% system and path
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
    path_name = '/work/Mengfan/Embryo/22-01-11';
    source_data = 'myf5GFP-H2BmCherry.v1.h5';
    data_name = fullfile(path_name, 'deconvolution/deconvolution_myf5GFP-H2BmCherry.v1.h5');
    xml_name = fullfile(path_name, 'myf5GFP-H2BmCherry.v1.xml');
    target_folder = 'fusion_NoRefine';
else
    addpath D:\MatlabTools;
    path_name = 'H:\Embryo\TM0-49\';
end
xml = xml2struct(xml_name);
register_info = xml.SpimData.ViewRegistrations.ViewRegistration;
[h5_struct, num_view, name_view] = readh5info(fullfile(path_name, source_data));
num_time = length(h5_struct);
num_total = num_time * num_view;
if num_view ~= 8
    error('Current version is for 8 views.');
end
%% parameter setting
device = 'GPU'; % GPU 25s CPU 45s parallel 250s
save_mode = 'both';

downsample_scale = 2;
pad_size = ceil(40 / downsample_scale); % fusion refinement constriant
bd_size = 3; % related to the deconvolution (n-1)/2. remove boundary
maxIter = 1000; % max iteration of step 4
baseIntensity = 200; % may not be the best threhsold, refine later

target_folder = [target_folder '_' num2str(1/downsample_scale)];
%% fusion
if ~isfolder(fullfile(path_name, target_folder))
    mkdir(fullfile(path_name, target_folder));
end
% loss_all = zeros(num_time,1);
% intensity_all = zeros(num_time,1); % for measure the loss and intensity
tic;
for tt = 195:195 %0:num_time-1
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
        vv_ind = name_view{vv};
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
        vv_ind = name_view{vv};
        im{vv} = hdf5read(data_name,['/t' tt_ind '/s' vv_ind '/0/cells']);
        im_size = size(im{vv});

        vv_ind = name_view{vv+4};
        im2 = hdf5read(data_name,['/t' tt_ind '/s' vv_ind '/0/cells']);
        im{vv} = (im{vv} + im2)/2;
        im{vv} = im{vv} - baseIntensity;
        im{vv}(im{vv}<0) = 0;
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

    % 4. registration refinement (gradient descent mean squared error)
    order = [1 2;2 3; 3 4; 4 1];
    conn = [0 0 1; 0 0 -1; 0 1 0; 0 -1 0; 1 0 0; -1 0 0;];
    diff = zeros(4,3);
    for oo = 1:4
        data1 = im_fuse{order(oo,1)};
        data2 = im_fuse{order(oo,2)};
        
        xs = max(x_start{order(oo,1)},x_start{order(oo,2)}) + pad_size;
        xe = min(x_end{order(oo,1)},x_end{order(oo,2)}) + pad_size;
        ys = max(y_start{order(oo,1)},y_start{order(oo,2)}) + pad_size;
        ye = min(y_end{order(oo,1)},y_end{order(oo,2)}) + pad_size;
        zs = max(z_start{order(oo,1)},z_start{order(oo,2)}) + pad_size;
        ze = min(z_end{order(oo,1)},z_end{order(oo,2)}) + pad_size;


        loss = inf;
        for iter = 1:maxIter
            diff_old = diff(oo,:);
            diff_candi = diff(oo,:) + conn;
            diff_candi(max(diff_candi,[],2)>pad_size,:) = [];
            diff_candi(min(diff_candi,[],2)<-pad_size,:) = [];

            loss_temp = zeros(size(diff_candi,1),1);               
            im1 = data1(ys:ye,xs:xe,zs:ze);
            for vv = 1:size(diff_candi,1)
                im2 = data2(ys+diff_candi(vv,2):ye+diff_candi(vv,2),...
                    xs+diff_candi(vv,1):xe+diff_candi(vv,1),zs+diff_candi(vv,3):ze+diff_candi(vv,3));
                loss_temp(vv) = mean((im1-im2).^2,'all');
            end
            [loss_min, ind_min] = min(loss_temp);
            if loss_min > loss
                break;
            else
                diff(oo,:) = diff_candi(ind_min,:);
                loss = loss_min;
            end
%             fprintf('Iter: %d Loss: %f\n', iter, loss);   % for debug
        end

%         fprintf('Pair: %d %d Diff: %d %d %d Iter: %d\n', order(oo,1), order(oo,2), diff(oo,1), diff(oo,2), diff(oo,3), iter);
    end
    clear data1 data2 im1 im2 % release gpu memory
    fprintf('Minimize MSE refinement result:\n');
    disp(diff);

    trans_modification = [0 0 0;diff;];
    trans_modification = cumsum(trans_modification);
    trans_modification(end,:) = [];
    clear im_fuse

    % 6. final fusion
    x_min = inf; y_min = inf; z_min = inf;
    x_max = -inf;y_max = -inf;z_max = -inf;
    for vv = 1:4 % update transform matrix
        trans_temp = eye(4,4);
        trans_temp(1:3,4) = -trans_modification(vv,1:3);
        trans{vv} = trans_temp*trans{vv};
        ref{vv} = affineOutputView(size(im{vv}), affine3d(trans{vv}'),'BoundsStyle','FollowOutput');

        x_min = min(x_min,round(ref{vv}.XWorldLimits(1)));
        x_max = max(x_max,round(ref{vv}.XWorldLimits(2)));
       
        y_min = min(y_min,round(ref{vv}.YWorldLimits(1)));
        y_max = max(y_max,round(ref{vv}.YWorldLimits(2)));
        
        z_min = min(z_min,round(ref{vv}.ZWorldLimits(1)));
        z_max = max(z_max,round(ref{vv}.ZWorldLimits(2)));
    end

    if ~strcmp(save_mode, 'colortif')
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
    else
        color = hsv2rgb([0 1 0.5; 0.25 1 0.5; 0.5 1 0.5; 0.75 1 0.5]);
        im_fusion = zeros(y_max-y_min,x_max-x_min,3,z_max-z_min,'single');
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
            if strcmp(device, 'GPU')
                im{vv} = gather(im{vv});
                im_unit = gather(im_unit);
            end

            x_start = round(ref{vv}.XWorldLimits(1) - x_min + 1);
            x_end = round(ref{vv}.XWorldLimits(2) - x_min);
            y_start = round(ref{vv}.YWorldLimits(1) - y_min + 1);
            y_end = round(ref{vv}.YWorldLimits(2) - y_min);
            z_start = round(ref{vv}.ZWorldLimits(1) - z_min + 1);
            z_end = round(ref{vv}.ZWorldLimits(2) - z_min);
            
            for zz = z_start:z_end
                im_fusion(y_start:y_end,x_start:x_end,1,zz) = im_fusion(y_start:y_end,x_start:x_end,1,zz) + im{vv}(:,:,zz-z_start+1) * color(vv,1);
                im_fusion(y_start:y_end,x_start:x_end,2,zz) = im_fusion(y_start:y_end,x_start:x_end,2,zz) + im{vv}(:,:,zz-z_start+1) * color(vv,2);
                im_fusion(y_start:y_end,x_start:x_end,3,zz) = im_fusion(y_start:y_end,x_start:x_end,3,zz) + im{vv}(:,:,zz-z_start+1) * color(vv,3);
            end
            im_num(y_start:y_end,x_start:x_end,z_start:z_end) = ...
            im_num(y_start:y_end,x_start:x_end,z_start:z_end) + im_unit;      
        end
        for zz = 1:size(im_fusion,3)
            im_fusion(:,:,1,zz) = im_fusion(:,:,1,zz)./im_num(:,:,zz);
            im_fusion(:,:,2,zz) = im_fusion(:,:,2,zz)./im_num(:,:,zz);
            im_fusion(:,:,3,zz) = im_fusion(:,:,3,zz)./im_num(:,:,zz);
        end
        im_fusion(isnan(im_fusion)) = 0;
        tifwrite(uint8(im_fusion),fullfile(path_name, target_folder, ['color_' tt_ind]));
    end
    clear im im_fusion im_num im_unit
    toc
end