clc;clear;close all;

% It's old version. Please see ../src_code_proProcessing/Registration.m
% register label map with given transformation matrix
% input:
% refine_res ------- logical foreground
% threshold_res ---- threshold to get the foreground
% tform ------- transformation matrix
% bound ------- maximum bound
% ref --------- image bound reference
% output:
% refine_reg
% threshold_reg
% contact: mengfanw@vt.edu, 05/24/2022

% if isunix
%     a = 1;
% else
%     load('I:\Embryo\TM0-49\foreground_crop\yinan_data_mengfan_fuse_res\synQuant_refine_res_4d_v9.mat');
%     load('I:\Embryo\TM0-49\foreground_crop\yinan_data_mengfan_fuse_res\ImageReg.mat');
% end


load('E:\Embryo\registration_temporal_data\ImageReg.mat');
folder_name = 'E:\Embryo\TM0-49\track_v1\';
load([folder_name 'synQuant_refine_res_4d_v9.mat']);

x_max = bound.x_max; x_min = bound.x_min;
y_max = bound.y_max; y_min = bound.y_min;
z_max = bound.z_max; z_min = bound.z_min;
refine_reg = cell(40,1);
threshold_reg = cell(40,1);
% load('.\tempData\synQuant_refine_res_reg.mat');
for tt = 1:40
    tt
% tt = 2;
x_start = round(ref{tt}.XWorldLimits(1) - x_min + 1);
x_end = round(ref{tt}.XWorldLimits(2) - x_min);
y_start = round(ref{tt}.YWorldLimits(1) - y_min + 1);
y_end = round(ref{tt}.YWorldLimits(2) - y_min);
z_start = round(ref{tt}.ZWorldLimits(1) - z_min + 1);
z_end = round(ref{tt}.ZWorldLimits(2) - z_min);
id_num = max(refine_res{tt}(:));
refine_temp = zeros(y_end - y_start + 1, x_end - x_start + 1, z_end - z_start + 1,id_num,'single');

tic;
parfor ii = 1:id_num
%     fprintf(' %d', ii);
    refine_temp2 = imwarp(double(refine_res{tt} == ii),affine3d(trans_mat{tt}));
    refine_temp(:,:,:,ii) = single(refine_temp2);
end
% refine_temp(refine_temp<0.3) = 0;
% [~, refine_temp] = max(refine_temp,[],4);
refine_temp_max = max(refine_temp,[],4);
refine_temp_ind = zeros(size(refine_temp_max));
for ii = 1:id_num
    refine_temp2 = refine_temp(:,:,:,ii);
    refine_temp_ind(refine_temp2 > 0.5 & refine_temp2 >= refine_temp_max) = ii;
end
clear refine_temp
toc
refine_reg{tt} = uint32(zeros(y_max - y_min, x_max - x_min, z_max - z_min));
refine_reg{tt}(y_start:y_end, x_start:x_end, z_start:z_end) = uint32(refine_temp_ind);
threshold_reg{tt} = uint8(zeros(y_max - y_min, x_max - x_min, z_max - z_min));
for ii = 1:id_num
    thre_ii = unique(threshold_res{tt}(refine_res{tt}==ii));
    if length(thre_ii) > 1
        error('Warning! Multiple thresholds!');
    end
    threshold_reg{tt}(refine_reg{tt}==ii) = thre_ii;
end
end
save([folder_name 'synQuant_refine_res_reg.mat'],'refine_reg','threshold_reg','-v7.3');

%%
clc;clear;close all;
folder_name = 'E:\Embryo\TM0-49\track_v1\';
load([folder_name 'varianceMap.mat']);
load([folder_name 'synQuant_priCvt_res.mat']);
load('E:\Embryo\registration_temporal_data\ImageReg.mat');

x_max = bound.x_max; x_min = bound.x_min;
y_max = bound.y_max; y_min = bound.y_min;
z_max = bound.z_max; z_min = bound.z_min;

eig_res_2d_reg = eig_res_2d;
eig_res_3d_reg = eig_res_3d;
varMap_reg = varMap;
for tt = 1:40
    tt
    x_start = round(ref{tt}.XWorldLimits(1) - x_min + 1);
    x_end = round(ref{tt}.XWorldLimits(2) - x_min);
    y_start = round(ref{tt}.YWorldLimits(1) - y_min + 1);
    y_end = round(ref{tt}.YWorldLimits(2) - y_min);
    z_start = round(ref{tt}.ZWorldLimits(1) - z_min + 1);
    z_end = round(ref{tt}.ZWorldLimits(2) - z_min);
    
    eig_temp = imwarp(eig_res_2d{tt},affine3d(trans_mat{tt}));
    eig_res_2d_reg{tt} = zeros(y_max - y_min, x_max - x_min, z_max - z_min, 'single');
    eig_res_2d_reg{tt}(y_start:y_end, x_start:x_end, z_start:z_end) = single(eig_temp);
   
    eig_temp = imwarp(eig_res_3d{tt},affine3d(trans_mat{tt}));
    eig_res_3d_reg{tt} = zeros(y_max - y_min, x_max - x_min, z_max - z_min, 'single');
    eig_res_3d_reg{tt}(y_start:y_end, x_start:x_end, z_start:z_end) = single(eig_temp);
    
    var_temp = fillmissing(fillmissing(varMap{tt}{1,1},'nearest',1),'nearest',2);
    var_temp = imwarp(var_temp ,affine3d(trans_mat{tt}));
    varMap_reg{tt}{1,1} = nan(y_max - y_min, x_max - x_min, z_max - z_min);
    varMap_reg{tt}{1,1}(y_start:y_end, x_start:x_end, z_start:z_end) = var_temp;
    
    var_temp = fillmissing(fillmissing(varMap{tt}{1,2},'nearest',1),'nearest',2);
    var_temp = imwarp(var_temp ,affine3d(trans_mat{tt}));
    varMap_reg{tt}{1,2} = nan(y_max - y_min, x_max - x_min, z_max - z_min);
    varMap_reg{tt}{1,2}(y_start:y_end, x_start:x_end, z_start:z_end) = var_temp;
end

save([folder_name 'synQuant_pri_var_reg.mat'], 'eig_res_2d_reg', 'eig_res_3d_reg', 'varMap_reg', '-v7.3');