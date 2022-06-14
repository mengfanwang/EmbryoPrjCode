% read all the gt files to generate groundtruth tracks

detection_path = 'crop_embryo_data_500x500x30x40';
load(fullfile(detection_path, 'embryo_image_detection.mat'), 'dets');
for i=1:numel(dets) % label each detection their frames
    dets{i} = cat(2, dets{i}, i*ones(size(dets{i},1),1));
end
det_all = cat(1, dets{:});

gt_im_path = '/home/congchao/Desktop/crop_embryo_data_500x500x30x40_k2_iter_2';

gt_id = [1 2 4 6 9 11 12 13 14 16 19 21 22 23 25 26 27 28 30 34 44 45 46 52 ...
    59 63 66 68 73 74 75 76 83 86 92 94 108 109 122 125 126 132 141 142 147 ...
    152 154 162 167 169 176 180 187 190 194 197 200 202 203 210 218 231 232 ...
    236 238 239 248 249 252 256 257 263 401 402 403 404 405 406 407 408 409 ...
    410 411 412 413 414 415 416 417 418 419 420 421 422 423 424 425 426 427 ...
    428];
gt_mat = cell(length(gt_id),1);
for i=1:length(gt_id)
    disp(i);
    
    tr_id = gt_id(i);
    im = tifread(fullfile(gt_im_path, ['track_',num2str(tr_id),'.tif']));
    [h,w,~,t] = size(im);
    if tr_id < 400
        cur_z = det_all(trajectories{tr_id}(1),3);
        fs = det_all(trajectories{tr_id},end);
    else
        cur_z = 0;
        fs = nan; % will not use for >400
    end
    
    track_inf = zeros(1e3,4);
    label_inf = zeros(1,2); % label number, frame,
    cnt = 0;
    for j=1:t
        % image info
        r = im(:,:,1,j);
        g = im(:,:,2,j);
        b = im(:,:,3,j);
        map = r==255 & g==255 & b==255;
        reg_pr = regionprops(map, 'Area', 'Centroid','BoundingBox');
        vv = [reg_pr(:).Area];
        [val, pos] = max(vv);
        
        if val>10 
            bb = reg_pr(pos).BoundingBox;
            bb(3) = bb(3) + bb(1);
            bb(4) = bb(4) + bb(2);
            val_pt = dets{j}(:,1) > bb(2) & dets{j}(:,1) < bb(4) & dets{j}(:,2) > bb(1) & dets{j}(:,2) < bb(3);
            loc_valpt = find(val_pt);
            if ~isempty(loc_valpt)
                if size(loc_valpt,1)>1
                    [~,tmp_pos] = min(abs(dets{j}(loc_valpt,3)-cur_z));
                    loc_valpt = loc_valpt(tmp_pos);
                end
                cnt = cnt + 1;
                track_inf(cnt,:) = [j, dets{j}(loc_valpt,1:3)];
                cur_z = track_inf(cnt,4);
            else
                cnt = cnt + 1;
                track_inf(cnt,:) = [j, reg_pr(pos).Centroid(2),reg_pr(pos).Centroid(1) cur_z];
            end
            label_inf(1) = label_inf(1) + 1;
            label_inf(2) = j;
        else
            % track infor
            loc = find(fs==j);
            if ~isempty(loc)
                cnt = cnt + 1;
                track_inf(cnt,:) = [j, det_all(trajectories{tr_id}(loc),1:3)];
                cur_z = track_inf(cnt,4);
            end
        end
    end
    loc_label = find(fs==label_inf(2));
    if label_inf(1)==1 && ~isempty(loc_label) % the label is end point
        loc_stop = find(track_inf(:,1)==label_inf(2));
        track_inf = track_inf(1:loc_stop-1,:);
    elseif tr_id==25
        loc_stop = find(track_inf(:,1)==label_inf(2));
        track_inf = track_inf(1:loc_stop,:);
        
    else
        track_inf = track_inf(1:cnt,:);
    end
    if tr_id > 400
        locs = find(track_inf(:,4) == 0);
        if ~isempty(locs)
            track_inf(locs, 4)  = track_inf(locs(end)+1, 4);
        end
    end
    gt_mat{i} = track_inf;
    
end

save(fullfile(detection_path,'gt.mat'),'gt_mat');

%% draw ground truth

max_rgb = zeros(h,w,3,t);
for j=1:t
    disp(j);
    f = j;
    
    orgIm3d = drawTrack_klb_loc(rgb_org{f}, gt_mat, p,f);
    tifwrite(orgIm3d, fullfile(save_folder, ['gt_tracks_',num2str(j)]));
    max_rgb(:,:,:,j) = max(orgIm3d,[],4);
end
tifwrite(max_rgb, fullfile(save_folder, 'gt_tracks'));


% write gt to folder
[ll,~] = cellfun(@size, gt_mat);
result_for_kitti_eval = zeros(sum(ll), 7);
cnt = 0;
for i=1:numel(gt_mat)
    for j=1:size(gt_mat{i},1)
        cnt = cnt + 1;
        result_for_kitti_eval(cnt, 1)   = gt_mat{i}(j,1);
        result_for_kitti_eval(cnt, 2)   = i;
        cc = gt_mat{i}(j,2:4); % simulate a bounding box: y, x, z;
        % only use x, y direction
        result_for_kitti_eval(cnt, 3:6) = max([cc(1)-2, cc(2)-2, 4 4], 1);
        result_for_kitti_eval(cnt, 7)   = 0.9;
    end
end
dlmwrite(fullfile(data_folder,'embryo_gt.txt'),result_for_kitti_eval,...
    'delimiter',' ', 'precision',8);



% for i=1:numel(gt_mat)
%     val_pt = gt_mat{i}(:,2) > 245 & gt_mat{i}(:,2) <262 & gt_mat{i}(:,3) > 376 & gt_mat{i}(:,3) < 394;
%     if sum(val_pt)>0
%         disp(gt_id(i));
%         break;
%     end
% end
cl_stack = zeros(h,w,z,t);
for i=1:t
    cl_stack(:,:,:,i) = squeeze(rgb_org{i}(:,:,1,:));

end
stack = reshape(cl_stack, h,w, []);
tifwrite(stack,fullfile(save_folder, 'embryo_3d_vid'));