% test the tracking results of mouse embryo data
addpath('C:\Users\Congchao\Dropbox\cc_ImHandle\');
addpath('C:\Users\Congchao\Desktop\Probjects_Google_Drive\ParticleTracking\src_code\');
addpath('C:\Users\Congchao\Desktop\Probjects_Google_Drive\ParticleTracking\src_uTrack\');
addpath('C:\Users\Congchao\Desktop\Probjects_Google_Drive\ParticleTracking\CINDA\');

% load detections from TGMM 2.0
detection_path = fullfile('crop_embryo_data_multi_samples','embryo_crop_411.mat');
load(detection_path, 'dets');
% build inputs
movieInfo = tree2tracks(dets, false);% false means we do not use existing tracking results
xCoord = cell(numel(dets),1);
yCoord = cell(numel(dets),1);
zCoord = cell(numel(dets),1);
for i=1:numel(dets)
    yCoord{i} = dets{i}(:,1);
    xCoord{i} = dets{i}(:,2);
    zCoord{i} = dets{i}(:,3);
end
% main function
[movieInfo,movieInfoAll] = mcfTracking(movieInfo, xCoord,yCoord,zCoord);

% save results
test_round = 3;
saveFolder4singleTrack = 'C:\Users\Congchao\Desktop\trajectory_samples';
video_folder = fullfile(saveFolder4singleTrack, ['embryo_crop_411_em_test_'...
    , num2str(test_round)]);
if ~exist(video_folder,'file')
    mkdir(video_folder);
end
load(detection_path, 'crop_embryo_vid');

movieInfo.xCoord = movieInfo.orgCoord(:,1);
movieInfo.yCoord = movieInfo.orgCoord(:,2);
movieInfo.zCoord = movieInfo.orgCoord(:,3);

drawTracksStart1st(crop_embryo_vid, dets,video_folder, movieInfo);