addpath('/home/ccw/Dropbox/cc_ImHandle/');
addpath('../../Tracking work related/codes/');

% load the detections data
data_folder = 'crop_embryo_data_500x500x30x40';
load(fullfile(data_folder,'embryo_image_detection.mat'),'dets');
for i=1:numel(dets) % label each detection their frames
    dets{i} = cat(2, dets{i}, i*ones(size(dets{i},1),1));
end
% build graph and solve with circulation framework
[dat_in, excess_node, c_en, c_ex] = build_embryo_graph_k2(dets);
g.n_nodes = excess_node(2);
g.c_en = c_en;
g.c_ex = c_ex;
newStr = 'embryo_crop';
out_name = writeCirculationInputFile(dat_in, g, newStr);
curPath = pwd;
cd /home/congchao/Downloads/cs2-master
[status, cmdOut] = system(['./cs2 < ', out_name]);
if status ~= 0
    error('cs2 function failed!\n');
end
cd(curPath);
[cost, dat1] = parsecs2OutputFile(0);
g.particleNum = g.n_nodes/2-1;
[trajectories, particle2track] = backtrack(dat1, g.particleNum, 1);
loopCnt = 0;
while loopCnt < 4
    [dat_in, excess_node, c_en, c_ex] = build_embryo_graph_iter_k2(dets, trajectories);
    g.n_nodes = excess_node(2);
    g.c_en = c_en;
    g.c_ex = c_ex;
    newStr = 'embryo_crop';
    out_name = writeCirculationInputFile(dat_in, g, newStr);
    curPath = pwd;
    cd /home/congchao/Downloads/cs2-master
    [status, cmdOut] = system(['./cs2 < ', out_name]);
    if status ~= 0
        error('cs2 function failed!\n');
    end
    cd(curPath);
    [cost, dat1] = parsecs2OutputFile(0);
    g.particleNum = g.n_nodes/2-1;
    [trajectories, particle2track] = backtrack(dat1, g.particleNum, 1);
    loopCnt = loopCnt + 1;
end
%% write to .txt file
ll = cellfun(@length, trajectories);
result_for_kitti_eval = zeros(sum(ll), 7);
cnt = 0;
track_cnt = 0;
for i=1:numel(trajectories)
    if length(trajectories{i}) < 5
        continue;
    end
    track_cnt = track_cnt + 1;
    for j=1:size(trajectories{i},1)
        cnt = cnt + 1;
        cur_pt = trajectories{i}(j);
        result_for_kitti_eval(cnt, 1)   = det_all(cur_pt, end);
        result_for_kitti_eval(cnt, 2)   = track_cnt;
        cc = det_all(cur_pt,1:3); % simulate a bounding box: y, x, z;
        % only use x, y direction
        result_for_kitti_eval(cnt, 3:6) = max([cc(1)-2, cc(2)-2, cc(1)+2, cc(2)+2], 0);
        result_for_kitti_eval(cnt, 7)   = 0.9;
    end
end
result_for_kitti_eval = result_for_kitti_eval(1:cnt,:);
invalid = isnan(result_for_kitti_eval(:,2));
result_for_kitti_eval(invalid,:) = [];

res_folder = '/home/congchao/Desktop/Tracking work related/MOTBeyondPixels-master/results/train';
dlmwrite(fullfile(res_folder,'iter_res_10.txt'),result_for_kitti_eval,'delimiter',' ', 'precision',8);



%% write results to images
save_folder = fullfile('/home/congchao/Desktop/', [data_folder,'_k2_iter_2']);
load(fullfile(data_folder,'embryo_image_detection.mat'), 'crop_embryo_vid');
if ~isfolder(save_folder)
    mkdir(save_folder);
end

det_all = cat(1, dets{:});
p = []; % particle and line infor
p.particleSize = 2;
p.lineWidth = 1;
p.cmap = hsv(200);
p_blank = p;
p_blank.cmap = p_blank.cmap(100:end,:);

% build the rgb graph cells
crop_embryo_vid = crop_embryo_vid./max(crop_embryo_vid(:));
[h,w,z,t] = size(crop_embryo_vid);
rgb_org = cell(t,1);
particleSize = p.particleSize;
parfor i=1:t
    disp(i);
    tmp_rgb = zeros(h,w,3,z);
    for j=1:z
        tmp_rgb(:,:,1,j) = crop_embryo_vid(:,:,j,i);
        tmp_rgb(:,:,2,j) = crop_embryo_vid(:,:,j,i);
        tmp_rgb(:,:,3,j) = crop_embryo_vid(:,:,j,i);
    end
    for j=1:size(dets{i})
        pt = dets{i}(j,1:3);
        tmp_rgb = draw3Dparticle(tmp_rgb,  pt, particleSize, [1 1 0]);
    end
    rgb_org{i} = tmp_rgb;
end

for i=1:numel(trajectories)
    if length(trajectories{i}) < 5
        continue;
    end
    disp(i);
    cur_folder = fullfile(save_folder, ['track_', num2str(i)]);
    if ~isfolder(cur_folder)
        mkdir(cur_folder);
    end
    if det_all(trajectories{i}(1),end) > 1
        f = det_all(trajectories{i}(1),end) - 1;
        tifwrite(rgb_org{f}, fullfile(cur_folder, ['track_',num2str(i),'_',num2str(f)]));
    end
    pre_f = det_all(trajectories{i}(1),end);
    for j=1:length(trajectories{i})
        cur_det = det_all(trajectories{i}(j),:);
        f = cur_det(end);
        if f-pre_f == 2
            orgIm3d = drawTrack_klb(rgb_org{f-1}, trajectories(i), det_all, p_blank,f-1);
            tifwrite(orgIm3d, fullfile(cur_folder, ['track_',num2str(i),'_',num2str(f-1)]));
        end
        orgIm3d = drawTrack_klb(rgb_org{f}, trajectories(i), det_all, p,f);
        tifwrite(orgIm3d, fullfile(cur_folder, ['track_',num2str(i),'_',num2str(f)]));
        pre_f = f;
    end
    if det_all(trajectories{i}(end),end) < t
        f = det_all(trajectories{i}(end),end) + 1;
        tifwrite(rgb_org{f}, fullfile(cur_folder, ['track_',num2str(i),'_',num2str(f)]));
    end
end

for i=1:numel(trajectories)
    if length(trajectories{i}) < 5
        continue;
    end
    disp(i);
    max_rgb = zeros(h,w,3,t);
    for j=1:t
        loct = find(det_all(trajectories{i},end) == j,1);
        if isempty(loct)
            max_rgb(:,:,:,j) = max(rgb_org{j}, [], 4);
        else
            f = j;
            orgIm3d = drawTrack_klb(rgb_org{j}, trajectories(i), det_all, p,f);
            max_rgb(:,:,:,j) = max(orgIm3d,[],4);
        end
    end
    tifwrite(max_rgb, fullfile(save_folder, ['track_',num2str(i)]));
end



