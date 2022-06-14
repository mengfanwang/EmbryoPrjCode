function [crop_embryo_vid, dets] = cropEmbryoData(scale, test_tps)
% INPUT:
% scale:     [y_min, y_max; x_min, x_max; z_min, z_max]
% test_tps:  [t_min t_max]
% OUTPUT:
% crop_embryo_vid: h*w*z*t matrix saving the cropped video
% dets: cells with t*1 size saving detection and linking results

test_tps = test_tps(1):test_tps(2);
crop_embryo_vid = zeros([scale(:,2)-scale(:,1)+1; length(test_tps)]');
for f = 1:length(test_tps)
    ff = test_tps(f);
    disp(ff);
    if isunix
        im_path = sprintf('/media/congchao/New Volume/ImageData/Mmu_E1_CAGTAG1.TM000%03d_timeFused_blending/',ff);
    else
        im_path = sprintf('E:\\Embryo_data\\Mmu_E1_CAGTAG1.TM000%03d_timeFused_blending\\',ff);
    end
    %im_name = 'SPM00_TM000082_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb';
    im_name = sprintf('SPM00_TM000%03d_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb',ff);
    [im, header] = readKLBstack(fullfile(im_path, im_name), 20);
    im = im(scale(1,1):scale(1,2),scale(2,1):scale(2,2),scale(3,1):scale(3,2));
%     im_double = double(im)/double(max(im(:)));
%     
%     tifwrite(im_double, fullfile(save_folder, sprintf('embryo_TM%03d',ff)));
    
    %     f_cnt = f_cnt + 1;
    
    crop_embryo_vid(:,:,:,f) = im;
end

dets = cell(length(test_tps),1);
load('embryo_10tb_detections.mat');
pre_id = [];
for f = 1:length(test_tps)
    disp(f);
    ff = test_tps(f);
    tmp_det = detections{ff};
    if isempty(tmp_det)
        dets{f} = [];
        continue;
    end
    vd_id = tmp_det(:,1) > scale(1,1) & tmp_det(:,1) < scale(1,2) & ...
        tmp_det(:,2) > scale(2,1) & tmp_det(:,2) < scale(2,2) &...
        tmp_det(:,3) > scale(3,1) & tmp_det(:,3) < scale(3,2);
    dets{f} = tmp_det(vd_id,:);
    dets{f}(:,1) = dets{f}(:,1) - scale(1,1) + 1;
    dets{f}(:,2) = dets{f}(:,2) - scale(2,1) + 1;
    dets{f}(:,3) = dets{f}(:,3) - scale(3,1) + 1;
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
%save(fullfile(save_folder,'embryo_image_detection.mat'),'crop_embryo_vid','dets','-v7.3');

end