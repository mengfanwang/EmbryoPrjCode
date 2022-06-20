clc;clear;close all;

%% system and path
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
    path_name = '/work/Mengfan/Embryo/TM0-49';
    source_data = 'H2BGFP_TM0-49.h5';
    target_folder = 'data';
else
    path_name = 'E:\Embryo\TM0-49';
    source_data = 'H2BGFP_TM0-49.h5';
    target_folder = 'data';
end

%% read info from h5 file
[h5_struct, num_view] = readh5info(fullfile(path_name, source_data));

% read h5 data and save as tif file
if ~isfolder(fullfile(path_name, target_folder))
    mkdir(fullfile(path_name, target_folder));
end
num_time = length(h5_struct);
for tt = 1:num_time
    fprintf('Processing time point %d/%d:', tt, num_time);
    tt_ind = num2str(99999+tt);
    tt_ind = tt_ind(2:6);
    mkdir(fullfile(path_name, target_folder, tt_ind));
    for vv = 1:num_view
        fprintf(' %d', vv);
        vv_ind = num2str(99+vv);
        vv_ind = vv_ind(2:3);
        data = hdf5read(fullfile(path_name, source_data),['/t' tt_ind '/s' vv_ind '/0/cells']);
        tifwrite(uint16(data),fullfile(path_name, target_folder,  tt_ind, vv_ind));
    end
    fprintf('\n');
end



