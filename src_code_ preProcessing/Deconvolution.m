clc;clear;close all;

%% system and path
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
    path_name = '/work/Mengfan/Embryo/22-01-11';
    source_data = 'myf5GFP-H2BmCherry.v1.h5';
    target_folder = 'deconvolution';    
else
    addpath D:\MatlabTools;
    path_name = '\\rs0001\Mengfan\Embryo\TM0-49\';
    source_data = 'H2BGFP_TM0-49.h5';
end
[h5_struct, num_view, name_view] = readh5info(fullfile(path_name, source_data));
if ~isfolder(fullfile(path_name, target_folder))
    mkdir(fullfile(path_name, target_folder));
end
num_time = length(h5_struct);
num_total = num_time * num_view;
%% parameter setting
device = 'GPU'; % CPU or GPU
save_mode = 'both'; %save as 'tif', 'h5', or 'both' 

z = 2:6:38;
pdf_z = normpdf(z,20,8);
pdf_z = pdf_z/sum(pdf_z(:));
pdf_z = reshape(pdf_z,[1 1 numel(pdf_z)]);
lambda = 0.9;
num_iter = 30;

%%
fprintf('Start processing...\n');
tic;
for ii = 0:num_total-1
    tt = floor(ii/num_view);
    vv = mod(ii, num_view);
    % read data from h5 file
    tt_ind = num2str(100000+tt);
    tt_ind = tt_ind(2:6);
    fprintf('processing: %d %d\n',tt, vv);

    vv_ind = name_view{vv+1};
    data = hdf5read(fullfile(path_name, source_data),['/t' tt_ind '/s' vv_ind '/0/cells']);
    data = single(data);
    z_size = size(data,3);
    data_median = medfilt3(data,[3 3 1]);
    H = fftshift(psf2otf(pdf_z,[1 1 z_size]));
    
    if strcmp(device, 'GPU')
        H = gpuArray(H);
    end
    A = 1 - lambda*conj(H).*H;
    y = fftshift(fft(data_median,[],3),3);
    g = lambda*conj(H).*y;
    x = y;
    for iter = 1:num_iter
        x = A.*x + g;
        data_deconv = real(ifft(ifftshift(x,3),[],3));
        data_deconv = max(data_deconv,0);
        x = fftshift(fft(data_deconv,[],3),3);
    end
    if strcmp(device, 'GPU')
        data_deconv = gather(data_deconv);
    end
    if strcmp(save_mode, 'h5') || strcmp(save_mode, 'both')
        h5create(fullfile(path_name, target_folder, ['deconvolution_' source_data]),...
            ['/t' tt_ind '/s' vv_ind '/0/cells'],size(data_deconv),'Datatype','single');
        h5write(fullfile(path_name, target_folder, ['deconvolution_' source_data]),...
            ['/t' tt_ind '/s' vv_ind '/0/cells'],data_deconv);
    end
    if strcmp(save_mode, 'tif') || strcmp(save_mode, 'both')
        if vv == 0
            mkdir(fullfile(path_name, target_folder, tt_ind));
        end
        tifwrite(uint16(data_deconv),fullfile(path_name, target_folder, tt_ind, vv_ind));
    end
    toc
end
toc
