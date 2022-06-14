% This function crop our embryo data into batches. Temporally, it will 
% split the data every T time points; For each time window with size T,
% it will divide the data into 4 batches spatially in x and z direction 
% with user-defined size.

input_folder = '../data/embryo_data'; % original images in tif format with name embryo_TM***.tif (e.g. embryo_TM001.tif)
save_folder = '../data/embryo_data_partitions'; % output folder
if ~isfolder(save_folder)
    mkdir(save_folder);
end

org_ims = dir(fullfile(input_folder, '*.tif'));

crop_size_w = 366; w_ov = 50;
crop_size_z = 259; z_ov = 25;
T = 50;

f_num = 0;
curr_folder_f_num = 0;
while f_num < numel(org_ims)
    if curr_folder_f_num == 0
        if f_num > 5 % temporally 5 frames overlap
            f_num = f_num - 5;
        end
        folder_name = sprintf('TM%03d_%03d', f_num, f_num+49);
        if ~isfolder(fullfile(save_folder, folder_name))
            mkdir(fullfile(save_folder, folder_name));
        end
        fl = fullfile(save_folder, folder_name, 'frontleft');
        if ~isfolder(fl)
            mkdir(fl);
        end
        fr = fullfile(save_folder, folder_name, 'frontright');
        if ~isfolder(fr)
            mkdir(fr);
        end
        bl = fullfile(save_folder, folder_name, 'backleft');
        if ~isfolder(bl)
            mkdir(bl);
        end
        br = fullfile(save_folder, folder_name, 'back_right');
        if ~isfolder(br)
            mkdir(br);
        end
    end
    imname = fullfile(input_folder, org_ims(f_num+1).name);
    im = uint16(tifread(imname));
    
    
    out_im = im(:,1:crop_size_w, 1:crop_size_z);
    tifwrite(out_im, fullfile(fl, sprintf('embryo_TM%03d',f_num)));
    r_start = crop_size_w - w_ov + 1;
    out_im = im(:,r_start:end, 1:crop_size_z);
    tifwrite(out_im, fullfile(fr, sprintf('embryo_TM%03d',f_num)));
    b_start = crop_size_z - z_ov + 1;
    out_im = im(:,1:crop_size_w, b_start:end);
    tifwrite(out_im, fullfile(bl, sprintf('embryo_TM%03d',f_num)));
    out_im = im(:,r_start:end, b_start:end);
    tifwrite(out_im, fullfile(br, sprintf('embryo_TM%03d',f_num)));
    
    f_num = f_num + 1;
    curr_folder_f_num = curr_folder_f_num + 1;
    if curr_folder_f_num == T % one folder contains 50 frames
        curr_folder_f_num = 0;
    end
end