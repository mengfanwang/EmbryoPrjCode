clc;clear;close all;
% The first step of registration. Warp all views to the same coordinate.

file_path = 'H:\Embryo\TM0-49\';
xml = xml2struct([file_path 'H2BGFP_TM0-49.xml']);
register_info = xml.SpimData.ViewRegistrations.ViewRegistration;

% get transformation
num_angle = 8; num_time = 50;
for tt = 47:num_time-1
fprintf('\nprocessing %d',tt);
tt_ind = num2str(100+tt);
tt_ind = tt_ind(2:3);
mkdir([file_path 'registration_step1\' tt_ind]);

trans = {};
for ii = 1:4
    trans{ii} = eye(4,4);
    for jj = 1:length(register_info{ii}.ViewTransform)
        trans_temp = cellfun(@str2num, split(register_info{tt*8 + ii}.ViewTransform{jj}.affine.Text));
        trans_temp = [reshape(trans_temp,4,3)'; 0 0 0 1];
        trans{ii} = trans{ii}*trans_temp;
    end
    trans{ii}(1:3,:) = trans{ii}(1:3,:) .* 0.25;  % downsample
end


% get transformed info
ref = {};
for ii = 1:4
    im = load([file_path 'deconv\' tt_ind '\' num2str(ii-1)]);
    im = im.data_deconv;
    [~, ref{ii}] = imwarp(im,affine3d(trans{ii}'));
end
clear im
save([file_path 'registration_step1\' tt_ind '\' 'ref.mat'],'ref');
% load([file_path 'ref_0.25_modi.mat']);

% merge together
x_min = inf; y_min = inf; z_min = inf;
x_max = -inf;y_max = -inf;z_max = -inf;
for ii = 1:4
    x_min = min(x_min,round(ref{ii}.XWorldLimits(1)));
    x_max = max(x_max,round(ref{ii}.XWorldLimits(2)));
   
    y_min = min(y_min,round(ref{ii}.YWorldLimits(1)));
    y_max = max(y_max,round(ref{ii}.YWorldLimits(2)));
    
    z_min = min(z_min,round(ref{ii}.ZWorldLimits(1)));
    z_max = max(z_max,round(ref{ii}.ZWorldLimits(2)));
end

color = hsv2rgb([0 1 0.5; 0.25 1 0.5; 0.5 1 0.5; 0.75 1 0.5]);
im_multi = zeros(y_max-y_min,x_max-x_min,3,z_max-z_min,'single');
im_num = zeros(y_max-y_min,x_max-x_min,z_max-z_min,'single');
for ii = 1:4
    fprintf(' %d ', ii);
    im = load([file_path 'deconv\' tt_ind '\' num2str(ii-1)]);
    im = im.data_deconv;
    im_size = size(im);
    im2 = load([file_path 'deconv\' tt_ind '\' num2str(ii+3)]);
    im2 = im2.data_deconv;
    im = (im + im2)/2;
    im = im - 200;
    clear im2
    for zz = 1:size(im,3)
%         im(:,:,zz) = flip(im(:,:,zz)', 1);
        im(:,:,zz) = im(:,:,zz)';
    end
    im = imwarp(im,affine3d(trans{ii}'));
    
    
    x_start = round(ref{ii}.XWorldLimits(1) - x_min + 1);
    x_end = round(ref{ii}.XWorldLimits(2) - x_min);
    y_start = round(ref{ii}.YWorldLimits(1) - y_min + 1);
    y_end = round(ref{ii}.YWorldLimits(2) - y_min);
    z_start = round(ref{ii}.ZWorldLimits(1) - z_min + 1);
    z_end = round(ref{ii}.ZWorldLimits(2) - z_min);
%     im_multi(y_start:y_end,x_start:x_end,z_start:z_end) = im_multi(y_start:y_end,x_start:x_end,z_start:z_end) + im;
%     im = im + im_multi(y_start:y_end,x_start:x_end,z_start:z_end);
%     im_multi(y_start:y_end,x_start:x_end,z_start:z_end) = im;
    for zz = z_start:z_end
        im_multi(y_start:y_end,x_start:x_end,1,zz) = im_multi(y_start:y_end,x_start:x_end,1,zz) + im(:,:,zz-z_start+1) * color(ii,1);
        im_multi(y_start:y_end,x_start:x_end,2,zz) = im_multi(y_start:y_end,x_start:x_end,2,zz) + im(:,:,zz-z_start+1) * color(ii,2);
        im_multi(y_start:y_end,x_start:x_end,3,zz) = im_multi(y_start:y_end,x_start:x_end,3,zz) + im(:,:,zz-z_start+1) * color(ii,3);
    end
    
    im_fuse = zeros(size(im_num));
    im_fuse(y_start:y_end,x_start:x_end,z_start:z_end) = im;
    save([file_path 'registration_step1\' tt_ind '\' 'im_' num2str(ii)],'im_fuse','x_start','x_end','y_start','y_end','z_start','z_end');
    

    
    im_unit = ones(im_size,'single');
    im_unit = imwarp(im_unit,affine3d(trans{ii}'));
    im_num(y_start:y_end,x_start:x_end,z_start:z_end) = im_num(y_start:y_end,x_start:x_end,z_start:z_end) + im_unit;
%     clear im im_unit
    

    
end
% im_multi = im_multi./im_num;
for zz = 1:size(im_multi,3)
    im_multi(:,:,1,zz) = im_multi(:,:,1,zz)./im_num(:,:,zz);
    im_multi(:,:,2,zz) = im_multi(:,:,2,zz)./im_num(:,:,zz);
    im_multi(:,:,3,zz) = im_multi(:,:,3,zz)./im_num(:,:,zz);
end
    
im_multi(isnan(im_multi)) = 0;
tifwrite(uint8(im_multi),[file_path 'registration_step1\' tt_ind '\' 'multiview']);
toc
end




