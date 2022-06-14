% crop a small portion of the embryo data
addpath(genpath('klb_wrapper'));
%save_folder = '../data/downsample_crop_embryo_data_470x350x250x50';
save_folder = '../data/downsample_embryo_data_700x700x500_all';

if ~isfolder(save_folder)
    mkdir(save_folder);
end
crop_flag = true;
crop_scale4downsampledData = [230, -1; -1, -1; -1, -1];
time_scale = [0, 531];%[481 530];
maxIntensity = [];
test_tps = time_scale(1):time_scale(2);
yxz_scale = [0.333,0.333, 0.5];
for f = 1:length(test_tps)
    ff = test_tps(f);
    disp(ff);
    im_path = sprintf('E:/Embryo_data/Mmu_E1_CAGTAG1.TM000%03d_timeFused_blending/',ff);
    %im_name = 'SPM00_TM000082_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb';
    im_name = sprintf('SPM00_TM000%03d_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb',ff);
    [im, header] = readKLBstack(fullfile(im_path, im_name), 20);
    [h,w,z] = size(im);
    out_im = imresize3(im, [round(h*yxz_scale(1)),...
        round(w*yxz_scale(2)),round(z*yxz_scale(3))]);
    
    if crop_flag
        [h,w,z] = size(out_im);
        y0 = get_real_idx(crop_scale4downsampledData(1,1),1);
        y1 = get_real_idx(crop_scale4downsampledData(1,2),h);
        x0 = get_real_idx(crop_scale4downsampledData(2,1),1);
        x1 = get_real_idx(crop_scale4downsampledData(2,2),w);
        z0 = get_real_idx(crop_scale4downsampledData(3,1),1);
        z1 = get_real_idx(crop_scale4downsampledData(3,2),z);
        out_im = out_im(y0:y1, x0:x1, z0:z1);
    end
    tifwrite(out_im, fullfile(save_folder, sprintf('embryo_TM%03d',ff)));
    %tifwrite(im, fullfile(save_folder, sprintf('org_embryo_TM%03d',ff)));
end



function idx = get_real_idx(idx, backup)
    if idx < 0
        idx = backup
    end
end