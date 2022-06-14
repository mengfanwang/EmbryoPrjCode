% for all purpose tmperoary code/bug testing
ss = 6;
f1 = fmapComp(:,:,ss);
v1 = vid_sm(:,:,ss);
zz = v1(f1);
thres5p = zz(round(0.05*length(zz)))


f1(v1<=thres5p) = false;
figure;imshow(label2rgb3d(double(fmapComp(:,:,ss))+double(f1)))

figure;imshow(cat(3, f1*0, vid_sm(:,:,ss)/50, f1*0.6));


fid = fopen('graph_1.txt','r');
dat_in = zeros(100000,3);
cnt = 0;
tline = fgetl(fid);
while ischar(tline)
    tline = fgetl(fid);
    if strcmp(tline(1), 'a')
        curl = sscanf(tline,'a %d %d %f');
        cnt = cnt + 1;
        dat_in(cnt, :) = curl';
    end
end
fclose(fid);
dat_in=dat_in(1:cnt, :);

dat_in(:,3) = floor(dat_in(:,3)*1e13);
fileID = fopen('graph_1_cs2.txt','w');
for i=1:size(dat_in,1)
    fprintf(fileID,'a %d %d 0 1  %d\n', dat_in(i,1), dat_in(i,2), dat_in(i,3));
end
fclose(fileID);

p min 2034 9515
n 1 3
n 2034 -3



im1 = imread('C:\Users\Congchao\Google Drive\Reports\Phd_exams\Preliminary exam\figures\gap_testing_rmGaps.png');
seeds = im1(:,:,1) + im1(:,:,3);
seeds = bwlabel(seeds>0,8);
zzshow(seeds)
org_data = im1(:,:,2);
sigma = 3; % 2d principal curvature smooth scale is smaller
[eig2d, ~] = principalCv2d(double(org_data), ones(size(org_data)), sigma, ones(size(org_data)));

im1 = imread('C:\Users\Congchao\Google Drive\Reports\Phd_exams\Preliminary exam\figures\gap_testing_cellrefine.png');
fg = im1(:,:,3)>0;
eig2d(eig2d<0) = 0;
newLabel = regionGrow(seeds, ...
            scale_image(eig2d, 1e-3, 1),...
            fg, q.growConnectInTest, ...
            q.cost_design, false);

[h,w] = size(org_data);
zz = zeros(h,w, 3);
vid = scale_image(double(org_data),0,1);
zz(:,:,1) = (newLabel==1)*0.5;
zz(:,:,2) = vid;
zz(:,:,3) = (newLabel==2)*0.5;

figure;imshow(zz)
tifwrite(zz, 're_detect_cells');



%% debug
for i=1:numel(gt_tracks)
    if gt_tracks{i}(1,2) > 206 && gt_tracks{i}(1,2) < 240 &&...
            gt_tracks{i}(1,3) > 82 && gt_tracks{i}(1,3) < 100
        disp(i);
    end
end

for i=1:numel(gt_tracks_all)
    gt_tracks_all{i}(:,2:3) = gt_tracks_all{i}(:,2:3)*2;
end
for i=1:numel(gt_tracks_no_split)
    gt_tracks_no_split{i}(:,2:3) = gt_tracks_no_split{i}(:,2:3)*2;
end
for i=1:numel(gt_track_split_parent)
    gt_track_split_parent{i}(:,2:3) = gt_track_split_parent{i}(:,2:3)*2;
    gt_track_split_kids{i,1}(:,2:3) = gt_track_split_kids{i,1}(:,2:3)*2;
    gt_track_split_kids{i,2}(:,2:3) = gt_track_split_kids{i,2}(:,2:3)*2;
end

save('../data/crop_embryo_data_500x500x30x40/gt_tracks_tyxz_with_split.mat',...
    'gt_tracks_all', 'gt_tracks_no_split', 'gt_track_split_parent', ...
    'gt_track_split_kids');


profile on;
for i=1:10000
    if mod(i,1000) == 0
        disp(i);
    end
    zz = randi([1 250*250],1,700);
    [yy0, xx0] = ind2sub([250, 250], zz);
    [yy1, xx1] = ind2sub_direct([250, 250], zz);
    if max(abs(yy0-yy1)) ~= 0
        disp(i);
    end
end
profile viewer;


profile on;
for i=1:10000
    if mod(i,1000) == 0
        disp(i);
    end
    yy = randi([1 250],1,700);
    xx = randi([1 250],1,700);
    zz = randi([1 30],1,700);
    id1 = sub2ind([250, 250, 30], yy, xx, zz);
    id2 = sub2ind_direct([250, 250, 30], yy, xx, zz);
    if max(abs(id1-id2)) ~= 0
        disp(i);
    end
end
profile viewer;
%%
for i=1:numel(refine_resAll)
    disp(i);
    disp(unique(refine_resAll{i}{25}(movieInfo.voxIdx{3337})));
    disp(unique(refine_resAll{i}{25}(movieInfo.voxIdx{3457})));
end