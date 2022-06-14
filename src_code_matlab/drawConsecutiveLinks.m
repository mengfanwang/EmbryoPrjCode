function vidout = drawConsecutiveLinks(vidin, detections, tif_savefolder)
% INPUT:
% vid: h*w*z*t matrix saving the cropped video
% dets: cells of t by 1; each cell is n*4 with location and parents id
% e.g. parent id = 1 means its parent is the first detection in the
% tif_savefolder: if not empty means save the data into tif files, [] means
% return to vidout

% OUTPUT:
% vidout: if tif_flag is 1, vid is []; otherwise vid is h*w*3*z*t matrix save
% the colored image saving links among adjacent frames

vidin = double(vidin)./double(max(vidin(:)));
% trajectory number
tra_num = 0;
for i=1:numel(detections)
    tra_num = tra_num + sum(isnan(detections{i}(:,end)));
end

p = [];
p.particleSize = 2;
p.lineWidth = 0.7;
p.cmap = hsv(tra_num + 1);
p.cmap = uint8(p.cmap*255);

%pt_cl = p.cmap(end,:);
cnt = 0;% trajectory count
for f = 1:numel(detections)
    disp(f);
    im = vidin(:,:,:,f);
    im_fg = uint8(255*im);

    [h,w,z] = size(im_fg);
    orgIm3d = uint8(zeros(h,w,3,z));
    for i=1:z
        cur_im = im_fg(:,:,i);
        orgIm3d(:,:,1,i) = cur_im;
        orgIm3d(:,:,2,i) = cur_im;
        orgIm3d(:,:,3,i) = cur_im;
    end
    detections{f} = cat(2, detections{f}, nan(size(detections{f},1), 1));
    for i=1:size(detections{f},1)
        cur_pt = detections{f}(i,:);
        if isnan(cur_pt(4))
            cnt = cnt + 1;
            cl = p.cmap(cnt,:);
            detections{f}(i,5) = cnt;
        else
            parent_pt = detections{f-1}(cur_pt(4), :);
            tr_id = parent_pt(5);
            cl = p.cmap(tr_id,:);
            detections{f}(i,5) = tr_id;
            orgIm3d = draw3Dline(orgIm3d, parent_pt(1:3), cur_pt(1:3), p.lineWidth, cl);
        end
        orgIm3d = draw3Dparticle(orgIm3d,  cur_pt(1:3), p.particleSize, cl);
    end
    tifwrite(orgIm3d,fullfile(tif_savefolder, sprintf('adjFrameTrack_%03d',f)));
end