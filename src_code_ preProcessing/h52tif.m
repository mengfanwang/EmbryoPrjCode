clc;clear;close all;

%% system and path
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
    path_name = '/work/Mengfan/Embryo/23-11-01_Yinan';
    source_data = 'HRASGFP-H2BmChery_2.h5';
    target_folder = 'membrane_view1';
    mode = 'merge';  % if mode == merge, combine two views together
else
    path_name = 'E:\Embryo\TM0-49';
    source_data = 'H2BGFP_TM0-49.h5';
    target_folder = 'data';
end

%% read info from h5 file
[h5_struct, num_view, name_view] = readh5info(fullfile(path_name, source_data));

% read h5 data and save as tif file
if strcmp(mode, 'merge')
    target_folder = [target_folder, '_merge'];
end
num_time = length(h5_struct);
for tt = 160:179%num_time
    fprintf('Processing time point %d/%d:', tt, num_time);
    tt_ind = num2str(99999+tt);
    tt_ind = tt_ind(2:6);
%     mkdir(fullfile(path_name, target_folder, tt_ind));
    if ~strcmp(mode, 'merge')
        for vv = 2:2 %1:8 %1:num_view
            vv_ind = name_view{vv};
            fprintf(vv_ind);
            data = hdf5read(fullfile(path_name, source_data),['/t' tt_ind '/s' vv_ind '/0/cells']);
            if ~isfolder(fullfile(path_name, [target_folder vv_ind]))
                mkdir(fullfile(path_name, [target_folder vv_ind]));
            end
            tifwrite(uint16(data),fullfile(path_name, [target_folder vv_ind],  tt_ind));
        end
        fprintf('\n');
    else
        for vv = 9:12 %1:num_view/2
            fprintf(' %d', vv);
            vv_ind = name_view{vv+4};
            data2 = hdf5read(fullfile(path_name, source_data),['/t' tt_ind '/s' vv_ind '/0/cells']);
            vv_ind = name_view{vv};
            data = hdf5read(fullfile(path_name, source_data),['/t' tt_ind '/s' vv_ind '/0/cells']);
            weight_1 = repmat([1:-1/1919:0]',1,1920);
            weight_2 = repmat([0:1/1919:1]',1,1920);
            data = (double(data).*weight_1 + double(data2).*weight_2)/2;
            tifwrite(uint16(data),fullfile(path_name, target_folder,  tt_ind, vv_ind));
        end
        fprintf('\n');
    end
end



