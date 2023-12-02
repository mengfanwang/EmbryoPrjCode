clc;clear;close all;

if isunix
    data_folder = '/work/Mengfan/Embryo/20220930_Joaquin/view9';
    res_folder = '/work/Mengfan/Embryo/20220930_Joaquin/Detection_view9/Wei_refine_res';
    save_folder = '/work/Mengfan/Embryo/20220930_Joaquin/Detection_view9/Wei_refine_res_tif_new';
end

%% check image num
tif_files = dir(fullfile(data_folder ,'/*.tif'));
mat_files = dir(fullfile(res_folder ,'/*.mat'));

%% merge syquant_res
if ~isfolder(save_folder)
    mkdir(save_folder);
end
for ii = 1:80
    ii
    org_im = tifread(fullfile(data_folder, tif_files(ii).name));
    mat_temp = load(fullfile(res_folder, mat_files(ii).name));
    [~, org_name, ~] = fileparts(tif_files(ii).name);
%     org_im(org_im>255) = 255;
    org_im = org_im-200;
    labelwrite(uint8(org_im/2), mat_temp.refine_res, fullfile(save_folder, org_name));
end
