% crop a small portion of the embryo data
addpath(genpath('klb_wrapper'));
save_folder = 'crop_embryo_data_500x500x30x40';
% f_cnt = 0;
%scale = [761 1260; 301 1400; 588 687; 481 520];
%scale = [761 1260; 901 1400; 658 687; 481 520];
scale = [761 1260; 901 1400; 658 687; 481 520];
crop_embryo_vid = zeros([scale(:,2)-scale(:,1)]'+1);
maxIntensity = [];
test_tps = scale(4,1):scale(4,2);
for f = 1:length(test_tps)
    ff = test_tps(f);
    disp(ff);
    im_path = sprintf('E:/Embryo_data/Mmu_E1_CAGTAG1.TM000%03d_timeFused_blending/',ff);
    %im_name = 'SPM00_TM000082_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb';
    im_name = sprintf('SPM00_TM000%03d_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb',ff);
    [im, header] = readKLBstack(fullfile(im_path, im_name), 20);
    im = im(scale(1,1):scale(1,2),scale(2,1):scale(2,2),scale(3,1):scale(3,2));
    maxIntensity(f) = (max(im(:)));
    %im_double = double(im)/double(max(im(:)));
    disp(maxIntensity(f));
    tifwrite(im, fullfile(save_folder, sprintf('embryo_TM%03d',ff)));
    
%     f_cnt = f_cnt + 1;
    
    crop_embryo_vid(:,:,:,f) = im;
end

dets = cell(length(test_tps),1); % y, x, z
load('embryo_10tb_detections.mat');
pre_id = [];
for f = 1:length(test_tps)
    disp(f);
    ff = test_tps(f);
    tmp_det = detections{ff};
    vd_id = tmp_det(:,1) > scale(1,1) & tmp_det(:,1) < scale(1,2) & ...
        tmp_det(:,2) > scale(2,1) & tmp_det(:,2) < scale(2,2) &...
        tmp_det(:,3) > scale(3,1) & tmp_det(:,3) < scale(3,2);
    dets{f} = tmp_det(vd_id,:);
    dets{f}(:,1) = dets{f}(:,1) - scale(1,1);
    dets{f}(:,2) = dets{f}(:,2) - scale(2,1);
    dets{f}(:,3) = dets{f}(:,3) - scale(3,1);
    if isempty(pre_id)
        dets{f}(:,4) = nan; % no valid parent
    else
        for i=1:size(dets{f},1)
            locs = find(pre_id(:,1) == dets{f}(i,4),1);
            if isempty(locs)
                dets{f}(i,4) = nan; % no valid parent
            else
                dets{f}(i,4) = pre_id(locs,2);
            end
        end
    end
    tmp_ids = find(vd_id);
    pre_id = [tmp_ids-1, [1:length(tmp_ids)]'];
end
save(fullfile(save_folder,'embryo_image_detection.mat'),'crop_embryo_vid','dets','-v7.3');