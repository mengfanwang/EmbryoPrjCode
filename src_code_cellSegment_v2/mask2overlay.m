clc;clear;close all;

if isunix
    data_folder = '/work/public/Embryo/Fluo-N3DL-TRIF_train/02_data';
    res_folder = '/work/public/cellTrackingChallenge/Fluo-N3DL-TRIF/02_RES';
    save_folder = '/work/public/Embryo/Fluo-N3DL-TRIF_train/previousTunedResult/02_RES_winTest/tif';
    save_mat_folder = '/work/public/Embryo/Fluo-N3DL-TRIF_train/previousTunedResult/02_RES_winTest/';
end

%% check image num
tif_files = dir(fullfile(data_folder ,'/*.tif'));
mask_files = dir(fullfile(res_folder ,'/*.tif'));

%% merge syquant_res
if ~isfolder(save_folder)
    mkdir(save_folder);
end
for ii = 1:10:79
    ii
    org_im = tifread(fullfile(data_folder, tif_files(ii).name));
    mask = tifread(fullfile(res_folder, mask_files(ii).name));
    [~, org_name, ~] = fileparts(tif_files(ii).name);
%     org_im(org_im>255) = 255;
    org_im = org_im-200;
    labelwrite(uint8(org_im/10), mask, fullfile(save_folder, org_name));
end
