clc;clear;close all;

fore_path = 'E:\Embryo\TM0-49\fore_v1\result\v2\';
data_path = 'E:\Embryo\TM0-49\fore_v1\input\';
fore_name_list = {'synQuant_refine_res_4d_v9'};
tif_files = dir(fullfile(data_path, '/*.tif'));
t = length(tif_files);

for nn = 1:length(fore_name_list)
    fore_name = fore_name_list{nn};
    fprintf('Processing %s\n', fore_name);
    
    load([fore_path fore_name '.mat']);
    data = zeros(200,350,50,t);
    label = zeros(size(data));
    for tt = 1:t
        data(:,:,:,tt) = tifread([data_path tif_files(tt).name]);
        if exist('id_mat','var')
            label(:,:,:,tt) = id_mat{tt};
        elseif exist('refine_res','var')
            label(:,:,:,tt) = refine_res{tt};
        end
    end
    if exist('id_mat','var')
        clear id_mat z_mat fMaps
    elseif exist('refine_res','var')
        clear refine_res threhsold_res
    end
    data(data>255) = 255;
    outIm = imdisplayWithMapColor4D(data,label);
    write4dTiffRGB(uint8(outIm),[fore_path fore_name '.tiff']);
end