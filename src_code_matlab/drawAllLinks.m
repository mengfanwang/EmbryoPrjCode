function vidout = drawAllLinks(vidin, detections, tif_savefolder)
% We save the trajectories starting at time point 1 as samples to show the
% performance of current testing method

% INPUT:
% vid: h*w*z*t matrix saving the cropped video
% dets: cells of t by 1; each cell is n*4 with location and parents id
% e.g. parent id = 1 means its parent is the first detection in the
% tif_savefolder: if not empty means save the data into tif files, [] means
% return to vidout

% OUTPUT:
% vidout: if tif_flag is 1, vid is []; otherwise vid is h*w*3*z*t matrix save
% the colored image saving links of all frames


vidin = double(vidin)./double(max(vidin(:)));
% trajectory number starting from frame 1
tra_num = sum(isnan(detections{1}(:,end)));

p = [];
p.particleSize = 2;
p.lineWidth = 0.7;
p.cmap = hsv(tra_num + 1); % particle color
p.cmap = uint8(p.cmap*255);

p.cmap = p.cmap(randperm(size(p.cmap,1)),:); % shuffle the color
pt_cl = p.cmap(end,:);
cnt = 0;% trajectory count
[h,w,z,~] = size(vidin);
labelMap = uint8(zeros(h,w,3,z));
track_pix = cell(tra_num, 1);
val_ids = [1:tra_num]';
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
    
    % remove those stopped ones
    if f>1
        parent_ids = detections{f}(:,end);
        %parent_ids = parent_ids(~isnan(parent_ids));
        %inval_ids = setdiff(val_ids, parent_ids);
        org_val_ids = val_ids;
        nxt_val_ids = true(size(detections{f},1),1);
        for i=1:length(nxt_val_ids)
            loc = find(val_ids==parent_ids(i),1);
            if isempty(loc)
                nxt_val_ids(i) = false;
            end
        end
        val_ids = find(nxt_val_ids);
        
        inval_ids = setdiff(org_val_ids, unique(parent_ids(nxt_val_ids)));
        for i=1:length(inval_ids)
            tr_id = detections{f-1}(inval_ids(i),end);
            labelMap(track_pix{tr_id}) = 0;
        end
        if ~isempty(val_ids)
            validIdx = labelMap>0;
            orgIm3d(validIdx) = labelMap(validIdx);
        end
    end
    detections{f} = cat(2, detections{f}, nan(size(detections{f},1), 1));
    for i=1:size(detections{f},1)
        cur_pt = detections{f}(i,:);
        if isnan(cur_pt(4))
            cnt = cnt + 1;
            detections{f}(i,5) = cnt;
            if cnt <= tra_num % starting from the first frame
                cl = p.cmap(cnt,:);
            else
                cl = pt_cl;
            end
            pix_l = {};% starting of the track, no line now
        else
            parent_pt = detections{f-1}(cur_pt(4), :);
            tr_id = parent_pt(5);
            if tr_id <= tra_num % starting from the first frame
                cl = p.cmap(tr_id,:);
                detections{f}(i,5) = tr_id;
                [orgIm3d, pix_l] = draw3Dline(orgIm3d, parent_pt(1:3), cur_pt(1:3), p.lineWidth, cl);
                labelMap(pix_l{1}) = cl(1);
                labelMap(pix_l{2}) = cl(2);
                labelMap(pix_l{3}) = cl(3);
            else
                cl = pt_cl;
            end
        end
        [orgIm3d, pix_p] = draw3Dparticle(orgIm3d,  cur_pt(1:3), p.particleSize, cl);
        if detections{f}(i,5) <= tra_num
            labelMap(pix_p{1}) = cl(1);
            labelMap(pix_p{2}) = cl(2);
            labelMap(pix_p{3}) = cl(3);
            
            tr_id = detections{f}(i,5);
            tmp_pix = cat(1, pix_p{:}, pix_l{:});
            track_pix{tr_id} = cat(1, track_pix{tr_id}, tmp_pix);
        end
    end
    tifwrite(orgIm3d,fullfile(tif_savefolder, sprintf('adjFrameTrack_%03d',f)));
end