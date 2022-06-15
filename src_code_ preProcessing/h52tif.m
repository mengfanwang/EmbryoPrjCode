clc;clear;close all;

%% system and path
if isunix
    folder_path = '\\rs0001\Mengfan\Embryo\TM0-49\';
    source_data = 'H2BGFP_TM0-49.h5';
    % target_data = 'H2BGFP_21-04-28_deconv.h5';
else
    folder_path = 'E:\Embryo\TM0-49\';
    source_data = 'H2BGFP_TM0-49.h5';
end

%% read info from h5 file
% warning: don't use h5info. It's 200 times slower.
h5_struct = hdf5info([folder_path source_data]);

% bdv data format: tTTTTT/sSS/L/cells
% tTTTTT: time points
% sSS: id of the setup (view)
% L: mipmap level (downsample exponent)



