clc;clear;close all;
% Using wavelet transform to fuse two images in the same view

%% system and path
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');

    path_name = '/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight';
    source_data = '20220518 isl2b H2Bmcherry overnight.h5';
    data_name = fullfile(path_name, 'deconvolution/deconvolution_20220518 isl2b H2Bmcherry overnight.h5');
%     data_name = fullfile(path_name, '20220920_isl2bGFP_H2BmCherry_6h_ON.h5');
    xml_name = fullfile(path_name, '20220518 isl2b H2Bmcherry overnight.xml');
    target_folder = 'view';
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
for tt = 0:19 %80:num_time
tt_ind = num2str(100000+tt);
tt_ind = tt_ind (2:6);
fprintf('processing: %d\n',tt);
% if ~isfolder(fullfile(path_name, target_folder))
%     mkdir(fullfile(path_name, target_folder));
% end
for vv = 12:12
    if ~isfolder(fullfile(path_name, [target_folder num2str(vv)]))
        mkdir(fullfile(path_name, [target_folder num2str(vv)]));
    end
    vv_ind = name_view{vv+4};
    im2 = hdf5read(data_name,['/t' tt_ind '/s' vv_ind '/0/cells']);

    vv_ind = name_view{vv};
    im1 = hdf5read(data_name,['/t' tt_ind '/s' vv_ind '/0/cells']);

    tic;
    im_fuse = sameViewFuse(im1, im2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    toc;
%     im_fuse = im_fuse - 200;
%     im_fuse(im_fuse<0) = 0;
    tifwrite(uint16(im_fuse),fullfile(path_name, [target_folder num2str(vv)], tt_ind));
%     tifwrite(uint16(im1),fullfile(path_name, target_folder, tt_ind));
end
end

