% crop a small portion of the embryo data
addpath(genpath('klb_wrapper'));
input_folder = '../data/downsample_embryo_data_700x700x500_all';
save_folder = '../data/embryo_data_partitions';
if ~isfolder(save_folder)
    mkdir(save_folder);
end

org_ims = dir(fullfile(input_folder, '*.tif'));

f_num = 0;
curr_folder_f_num = 0;
while f_num < numel(org_ims)
    if curr_folder_f_num == 0
        if f_num > 5
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
    
    out_im = im(:,1:366, 1:259);
    tifwrite(out_im, fullfile(fl, sprintf('embryo_TM%03d',f_num)));
    out_im = im(:,317:end, 1:259);
    tifwrite(out_im, fullfile(fr, sprintf('embryo_TM%03d',f_num)));
    out_im = im(:,1:366, 236:end);
    tifwrite(out_im, fullfile(bl, sprintf('embryo_TM%03d',f_num)));
    out_im = im(:,317:end, 236:end);
    tifwrite(out_im, fullfile(br, sprintf('embryo_TM%03d',f_num)));
    
    f_num = f_num + 1;
    curr_folder_f_num = curr_folder_f_num + 1;
    if curr_folder_f_num == 50 % one folder contains 50 frames
        curr_folder_f_num = 0;
    end
end