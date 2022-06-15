clc;clear;close all;

%% system and path
if isunix
    folder_path = '\\rs0001\Mengfan\Embryo\TM0-49\';
    source_data = 'H2BGFP_TM0-49.h5';
    % target_data = 'H2BGFP_21-04-28_deconv.h5';
else
    addpath D:\MatlabTools;
    folder_path = '\\rs0001\Mengfan\Embryo\TM0-49\';
    source_data = 'H2BGFP_TM0-49.h5';
end

%% parameter setting
num_worker = 0;

z = 2:6:38;
pdf_z = normpdf(z,20,8);
pdf_z = pdf_z/sum(pdf_z(:));
pdf_z = reshape(pdf_z,[1 1 numel(pdf_z)]);
z_size = size(data,3);
lambda = 0.9;
num_iter = 30;

%%
tic;
for tt = 0:49
% read data from h5 file
tt_ind = num2str(100+tt);
tt_ind = tt_ind(2:3);
mkdir([folder_path 'deconv\' tt_ind]);
fprintf('\nprocessing %d',tt);
for view = 0:7
    fprintf(' %d ', view);
    data = h5read([folder_path source_data],['/t000' tt_ind '/s0' num2str(view) '/0/cells']);

    data = double(data);
    data_median = medfilt3(data,[3 3 1]);
    


    H = fftshift(psf2otf(pdf_z,[1 1 z_size]));
    A = 1 - lambda*conj(H).*H;
    y = fftshift(fft(data_median,[],3),3);
    g = lambda*conj(H).*y;
    x = y;
    for iter = 1:num_iter
        fprintf(' %d ',iter);
        x = A.*x + g;
        data_deconv = real(ifft(ifftshift(x,3),[],3));
        data_deconv = max(data_deconv,0);
        x = fftshift(fft(data_deconv,[],3),3);
    end
    toc
    save([folder_path 'deconv\' tt_ind '\' num2str(view)],'data_deconv','-v7.3');
    tifwrite(uint16(data_deconv),[folder_path 'deconv\' tt_ind '\' num2str(view)]);
end
end

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