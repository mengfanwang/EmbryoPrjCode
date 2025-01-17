clc;clear;close all;
% Using wavelet transform to fuse two images in the same view

%% system and path
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
%     path_name = '/work/Mengfan/Embryo/22-01-11_Yinan';
%     source_data = 'myf5GFP-H2BmCherry.h5';
%     data_name = fullfile(path_name, source_data);
%     xml_name = fullfile(path_name, 'myf5GFP-H2BmCherry.xml');
%     target_folder = 'sameViewFusion_2';

    path_name = '/work/Mengfan/Embryo/23-11-01_Yinan';
    source_data = 'HRASGFP-H2BmChery_2.h5';
    data_name = fullfile(path_name, 'HRASGFP-H2BmChery_2.h5');
    xml_name = fullfile(path_name, 'HRASGFP-H2BmChery_2.xml');
    target_folder = 'view1';
else
    addpath D:\MatlabTools;
    path_name = 'H:\Embryo\TM0-49\';
end
xml = xml2struct(xml_name);
register_info = xml.SpimData.ViewRegistrations.ViewRegistration;
[h5_struct, num_view, name_view] = readh5info(fullfile(path_name, source_data));
num_time = length(h5_struct);
num_total = num_time * num_view;

%% fusion
for tt = 160:179
tt_ind = num2str(100000+tt);
tt_ind = tt_ind(2:6);
fprintf('processing: %d\n',tt);
if ~isfolder(fullfile(path_name, target_folder))
    mkdir(fullfile(path_name, target_folder));
end
for vv = 2:2
    vv_ind = name_view{vv+4};
    im2 = hdf5read(data_name,['/t' tt_ind '/s' vv_ind '/0/cells']);

    vv_ind = name_view{vv};
    im1 = hdf5read(data_name,['/t' tt_ind '/s' vv_ind '/0/cells']);

    tic;
    im_fuse = sameViewFuse(im1, im2);
    toc;
%     im_fuse = im_fuse - 200;
%     im_fuse(im_fuse<0) = 0;
%     tifwrite(uint16(im_fuse),fullfile(path_name, target_folder, tt_ind, vv_ind));
    tifwrite(uint16(im_fuse),fullfile(path_name, target_folder, tt_ind));
end
end

