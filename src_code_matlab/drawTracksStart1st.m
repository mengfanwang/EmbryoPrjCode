function drawTracksStart1st(vidin, save_folder,movieInfo, dataWithParticles, write_ids)
% draw the trajectories starting at the first time point


if ~exist('dataWithParticles','var')
    dataWithParticles=[];
end
if ~exist('write_ids','var')
    write_ids = [];
end
% remove the nan id in tracks
for i=1:numel(movieInfo.tracks)
    cur_track = movieInfo.tracks{i};
    cur_track(isnan(cur_track)) = [];
    movieInfo.tracks{i} = cur_track;
end
vidin = double(vidin)./(double(max(vidin(:)))/255);


numTracks = min(numel(movieInfo.tracks),sum(movieInfo.frames==1)); % number of trjectories starting at 1 time piont
p = track_disp_para(numTracks);
% generate the data with all particles
if isempty(dataWithParticles)
    dataWithParticles = generate3DImWithParticle(vidin, movieInfo, p);
end
% if ~exist(fullfile(save_folder,'dataWithParticles.mat'),'file')
%     tic;
%     dataWithParticles = generate3DImWithParticle(vidin, movieInfo, p);
%     toc;
%     %save(fullfile(save_folder,'dataWithParticles.mat'),'dataWithParticles','-v7.3');
% else
%     tic;
%     load(fullfile(save_folder,'dataWithParticles.mat'),'dataWithParticles');
%     toc;
% end
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
if isempty(write_ids)
    write_ids = 1:20;
end
for j = 1:length(write_ids)
    p.filePath = fullfile(save_folder, sprintf('track_%04d', write_ids(j)));

    p.fileName = 'frame';

    write3DVidwithOneTrack(vidin,movieInfo,write_ids(j),p,dataWithParticles);
end

end