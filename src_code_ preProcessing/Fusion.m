clc;clear;close all;

%% system and path
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
    path_name = '/work/Mengfan/Embryo/TM0-49';
    source_data = 'H2BGFP_TM0-49.h5';
    data_name = fullfile(path_name, 'deconvolution/deconvolution_H2BGFP_TM0-49.h5');
    xml_name = fullfile(path_name, 'H2BGFP_TM0-49.xml');
    target_folder = 'fusion';
else
    addpath D:\MatlabTools;
    path_name = 'H:\Embryo\TM0-49\';
end
xml = xml2struct(xml_name);
register_info = xml.SpimData.ViewRegistrations.ViewRegistration;
[h5_struct, num_view] = readh5info(fullfile(path_name, source_data));
if ~isfolder(fullfile(path_name, target_folder))
    mkdir(fullfile(path_name, target_folder));
end
num_time = length(h5_struct);
num_total = num_time * num_view;
%% parameter setting
downsample_scale = 4;

%%
% get transformation
for tt = 47:num_time-1
    tt_ind = num2str(100000+tt);
    tt_ind = tt_ind(2:6);
    fprintf('processing: %d\n',tt);

    trans = {};
    for vv = 1:4
        trans{vv} = eye(4,4);
        for jj = 1:length(register_info{vv}.ViewTransform)
            trans_temp = cellfun(@str2num, split(register_info{tt*8 + vv}.ViewTransform{jj}.affine.Text));
            trans_temp = [reshape(trans_temp,4,3)'; 0 0 0 1];
            trans{vv} = trans{vv}*trans_temp;
        end
        trans{vv}(1:3,:) = trans{vv}(1:3,:) / downsample_scale;  % downsample
    end
    
    
    % get transformed info
    ref = {};
    for vv = 1:4
        vv_ind = num2str(100+vv-1);
        vv_ind = vv_ind(2:3);
        im_size = h5info(data_name, ['/t' tt_ind '/s' vv_ind '/0/cells']);
        im_size = im_size.Dataspace.Size;
        ref{vv} = affineOutputView(im_size, affine3d(trans{vv}'),'BoundsStyle','FollowOutput');
    end
end