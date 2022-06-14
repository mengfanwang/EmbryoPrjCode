addpath('klb_wrapper\');
%% crop a small portion of the embryo data
save_folder = 'crop_embryo_data_multi_samples';
if ~exist(save_folder,'file')
    mkdir(save_folder);
end

% detail of the cropped region (spatial and temporal locations)
bs_scale = [761 1560; 351 1350; 401 600; 481 520];
detailed_sc = {[225 555; 220 785; 11 60], [220 570; 120 590; 1 50], [100 550; 350 730; 151 200],...
    [140 530; 1 450; 151 200],[120 580; 1 540; 151 200]};
test_tps = {[11 50], [111 150], [211 250], [311 350], [411 450]};
for i = 1:numel(test_tps)
    scale = detailed_sc{i} + bs_scale(1:3,1);
    tps = test_tps{i};
   [crop_embryo_vid, dets] = cropEmbryoData(scale, tps);
    save(fullfile(save_folder,sprintf('embryo_crop_%03d.mat', tps(1))),...
        'crop_embryo_vid','dets','-v7.3');
    clear crop_embryo_vid dets;
end
%% save each frame into a single frame
files = dir(fullfile(save_folder, '*.mat'));
% save several frames for desingning detection algorithm

for i=numel(files):-1:1
    disp(i);
    load(fullfile(files(i).folder, files(i).name));
    video_folder = fullfile(save_folder, files(i).name(1:end-4));
    if ~exist(video_folder,'file')
        mkdir(video_folder);
    end
    for f = 1:size(crop_embryo_vid,4)
        disp(f);
        im = squeeze(crop_embryo_vid(:,:,:,f));
        tifwrite(uint16(im),fullfile(video_folder, sprintf('detection_sample_%03d',f)));
    end
end


%% draw the tracking results and save in the corresponding folders
% contains all the linkings among "adjacent frames"
files = dir(fullfile(save_folder, '*.mat'));

for i=numel(files):-1:1
    disp(i);
    load(fullfile(files(i).folder, files(i).name));
    video_folder = fullfile(save_folder, files(i).name(1:end-4));
    if ~exist(video_folder,'file')
        mkdir(video_folder);
    end
    drawConsecutiveLinks(crop_embryo_vid, dets, video_folder);
end

%% draw the tracking results and save in the corresponding folders
% contains the intact trajectories that starts from the 1st frame

files = dir(fullfile(save_folder, '*.mat'));

for i=numel(files)% just test the data with time points: 411-450
    disp(i);
    load(fullfile(files(i).folder, files(i).name));
    video_folder = fullfile(save_folder, [files(i).name(1:end-4),'_intact_trace']);
    if ~exist(video_folder,'file')
        mkdir(video_folder);
    end
    drawAllLinks(crop_embryo_vid, dets, video_folder);
end

%% draw the tracking results and save in the corresponding folders
% contains the intact trajectories that starts from the 1st frame
% each trajectory saved in a folder

files = dir(fullfile(save_folder, '*.mat'));
saveFolder4singleTrack = 'C:\Users\Congchao\Desktop\trajectory_samples';
for i=numel(files)% just test the data with time points: 411-450
    disp(i);
    load(fullfile(files(i).folder, files(i).name));
    video_folder = fullfile(saveFolder4singleTrack, [files(i).name(1:end-4),'_all']);
    if ~exist(video_folder,'file')
        mkdir(video_folder);
    end
    drawTracksStart1st(crop_embryo_vid, dets,video_folder);
end