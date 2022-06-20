clc;clear;close all;

%% system and path
if isunix
    path_name = '/work/Mengfan/Embryo/TM0-49';
    source_data = 'H2BGFP_TM0-49.h5';
    target_folder = 'deconvolution';
else
    addpath D:\MatlabTools;
    path_name = '\\rs0001\Mengfan\Embryo\TM0-49\';
    source_data = 'H2BGFP_TM0-49.h5';
end
[h5_struct, num_view] = readh5info(fullfile(path_name, source_data));
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
% try to use gpu
fprintf('Start processing...\n');
tic;
for ii = 0:num_total-1
    tt = floor(ii/num_view);
    vv = mod(ii, num_view);
    % read data from h5 file
    tt_ind = num2str(100000+tt);
    tt_ind = tt_ind(2:6);
    fprintf('processing: %d %d\n',tt, vv);

    vv_ind = num2str(100+vv);
    vv_ind = vv_ind(2:3);
    data = hdf5read(fullfile(path_name, source_data),['/t' tt_ind '/s' vv_ind '/0/cells']);
    data = single(data);
    z_size = size(data,3);
    data_median = medfilt3(data,[3 3 1]);
    H = fftshift(psf2otf(pdf_z,[1 1 z_size]));
    
    if strcmp(device, 'GPU')
        gpuDevice(mod(ii,2)+1);
        H = gpuArray(H);
    end
    A = 1 - lambda*conj(H).*H;
    y = fftshift(fft(data_median,[],3),3);
    g = lambda*conj(H).*y;
    x = y;
    for iter = 1:num_iter
%             fprintf(' %d ',iter);
        x = A.*x + g;
        data_deconv = real(ifft(ifftshift(x,3),[],3));
        data_deconv = max(data_deconv,0);
        x = fftshift(fft(data_deconv,[],3),3);
    end
    if strcmp(device, 'GPU')
        data_deconv = gather(data_deconv);
    end
    if strcmp(save_mode, 'h5') || strcmp(save_mode, 'both')
        for factor = 0:3
            [x,y,z] = size(data_deconv);
            if factor == 3
                z = floor(z/2);
            end
            data = imresize3(data_deconv,[floor(x/2^factor)  floor(y/2^factor)  z]);
            h5write([file_path target_data],['/t00000/s0' num2str(view) '/' num2str(factor) '/cells'],int16(data));
        end
    end
    if strcmp(save_mode, 'tif') || strcmp(save_mode, 'both')
        if ~isfolder(fullfile(path_name, target_folder)) && tt == 0 && vv == 0
            mkdir(fullfile(path_name, target_folder));
        end
        if vv == 0
            mkdir(fullfile(path_name, target_folder, tt_ind));
        end
        tifwrite(uint16(data_deconv),fullfile(path_name, target_folder, tt_ind, vv_ind));
    end
end
toc

% %%
% for view = 0:7
%     view
%     load([file_path 'deconv\30\' num2str(view)]);
%     for factor = 0:3
%         [x,y,z] = size(data_deconv);
%         if factor == 3
%             z = floor(z/2);
%         end
%         data = imresize3(data_deconv,[floor(x/2^factor)  floor(y/2^factor)  z]);
%         h5write([file_path target_data],['/t00000/s0' num2str(view) '/' num2str(factor) '/cells'],int16(data));
% %         h5write([file_path target_data],['/t00000/s0' num2str(view) '/' num2str(factor) '/cells'],int16(1000*rand(size(data))));
%     end
% end