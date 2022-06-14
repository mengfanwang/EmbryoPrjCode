addpath('klb_wrapper');
if isunix
    addpath('/home/ccw/Dropbox/cc_ImHandle/');
    im_path = '/home/congchao/Desktop/ImageData/Mmu_E1_CAGTAG1.TM000394_timeFused_blending/';
else
    addpath('C:\Users\Congchao\Dropbox\cc_ImHandle\');
    im_path = 'E:\Embryo_data\Mmu_E1_CAGTAG1.TM000200_timeFused_blending\';
end
%im_name = 'SPM00_TM000082_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb';
%im_name = 'SPM00_TM000394_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb';
im_name = 'SPM00_TM000200_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb';
[im, header] = readKLBstack(fullfile(im_path, im_name), 4);
%'E:\Embryo_data\Mmu_E1_CAGTAG1.TM000000_timeFused_blending\SPM00_TM000000_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb'
save_f = 'C:\Users\Congchao\Desktop\cell_detection_samples\crop_embryo_data_500x500x30x40\track_visual';
tifwrite(im,fullfile(save_f, 't200_ch0_org'));

im_name = 'SPM00_TM000481_CM00_CM01_CHN01.fusedStack.corrected.shifted.klb';
[im_ch1, header] = readKLBstack(fullfile(im_path, im_name), 4);

%im_fg = tifread('/home/congchao/Desktop/t250_ch0_fg-1.tif');
im_fg = tifread('t250_ch0_fg.tif');
tifwrite(uint16(im_fg),'t250_ch0_fg_rm');
im_fg = uint8(255*im_fg/max(im_fg(:)));
luts;
[h,w,z] = size(im_fg);
out_cl_im = uint8(zeros(h,w,3,z));
for i=1:z
    disp(i);
    cur_im = im_fg(:,:,i);
    
    tmp = uint8(zeros(h,w));
    tmp(:) = royal(cur_im(:)+1,2);
    out_cl_im(:,:,1,i) = tmp;
    
    tmp(:) = royal(cur_im(:)+1,3);
    out_cl_im(:,:,2,i) = tmp;
    
    tmp = uint8(zeros(h,w));
    tmp(:) = royal(cur_im(:)+1,4);
    out_cl_im(:,:,3,i) = tmp;
end

tifwrite(out_cl_im,'t250_ch0_fg_small_cl');
