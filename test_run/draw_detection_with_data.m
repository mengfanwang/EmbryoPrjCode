if isunix
    addpath('/home/ccw/Dropbox/cc_ImHandle/');
else
    addpath('C:\Users\Congchao\Dropbox\cc_ImHandle');
end
addpath('../../particleTracking_linux');
% draw the detected cells with the org data
%data_path = 'crop_embryo_data_500x1100x100x40'; 
data_path = 'crop_embryo_data_800x1000x200x40_sample'; 
load(fullfile(data_path,'embryo_image_detection.mat'));
% use crop_embryo_data_500x500x30x40 as an example

if strcmp(data_path, 'crop_embryo_data_500x1100x100x40')
    scale = [761 1260; 301 1400; 588 687; 481 520];
    for f=1:numel(dets)
        dets{f}(:,1) = dets{f}(:,1) - scale(1,1) + 1;
        dets{f}(:,2) = dets{f}(:,2) - scale(2,1) + 1;
        dets{f}(:,3) = dets{f}(:,3) - scale(3,1) + 1;
    end
end
p = [];
p.particleSize = 2;
p.lineWidth = 1;
p.cmap = hsv(100);
p.cmap = uint8(p.cmap*255);
max_in = max(crop_embryo_vid(:));
tt = numel(dets);
for ff = 2:17
    im = crop_embryo_vid(:,:,:,ff);
    im_fg = uint8(255*im/max_in);

    [h,w,z] = size(im_fg);
    orgIm3d = uint8(zeros(h,w,3,z));
    for i=1:z
        disp(i);
        tmp = im_fg(:,:,i);

        orgIm3d(:,:,1,i) = tmp;
        orgIm3d(:,:,2,i) = tmp;
        orgIm3d(:,:,3,i) = tmp;
    end

    for j=1:size(dets{ff},1)
        if mod(j,30)==0
            disp([j size(dets{ff},1)]);
        end
        pt = dets{ff}(j,1:3);
        orgIm3d = draw3Dparticle(orgIm3d,  pt, p.particleSize, p.cmap(1,:));
    end
    tifwrite(orgIm3d,fullfile(data_path,sprintf('track_dispaly_%03d',10+(ff-1)*20)));

end
