%This is the main function for tracking results refinement
if isunix
    addpath('/home/ccw/Dropbox/cc_ImHandle/');
    save_folder = '/home/ccw/Desktop/embryo_res_folder';
else
    addpath('C:\Users\Congchao\Dropbox\cc_ImHandle\');
    save_folder = 'C:\Users\Congchao\Desktop\cell_detection_samples';
end

addpath('dt');
addpath(genpath('CINDA/'));
addpath('src_code_matlab');
addpath('src_code_cellTracker');
addpath('src_code_cellSegment');
addpath('debug_func/');
addpath('../../ParticleTracking/src_code/');
addpath('../../ParticleTracking/src_uTrack/');
% set saving folder
data_folder = '../data';
tif_files = dir(fullfile(data_folder, 'crop_embryo_data_500x500x30x40/embryo_TM*.tif'));
res_folder = fullfile(save_folder,'crop_embryo_data_500x500x30x40');
if ~exist(res_folder,'dir')
    mkdir(res_folder);
end
disp_folder = 'synQuant_refine_res';

% load original data and detection results and ground truth
if ~exist('refine_res','var') % segmentation results
    %load(fullfile(res_folder, 'synQuant_refine_res.mat'),'refine_res');
    load(fullfile(res_folder, 'synQuant_refine_res_4d_v9plus.mat'),...
        'refine_res', 'threshold_res');
    org_refine_res = refine_res;
    org_threshold_res = threshold_res;
end
scale_term = 5000;
embryo_vid_org = cell(numel(tif_files), 1); % original data
for i=1:numel(tif_files)
    embryo_vid_org{i} = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    embryo_vid_org{i} = 255*embryo_vid_org{i}./scale_term;
end

if ~exist('varMap','var') % segmentation results
    %load(fullfile(res_folder, 'synQuant_refine_res.mat'),'refine_res');
    load(fullfile(res_folder, 'varianceMap.mat'),'varMap');
    load(fullfile(res_folder, 'synQuant_priCvt_res.mat'),...
        'eig_res_2d', 'eig_res_3d');
    
    org_varMap = varMap;    
    org_eigMaps = cell(numel(eig_res_2d),1);
    for i=1:numel(org_eigMaps)
        org_eigMaps{i} = cell(2,1);
        org_eigMaps{i}{1} = eig_res_2d{i};
        org_eigMaps{i}{2} = eig_res_3d{i};
    end
end


if ~exist('gt_mat','var')% tracking groundtruth
    load('../data/crop_embryo_data_500x500x30x40/gt_tracks_tyxz_with_split.mat','gt_tracks_all');
    gt_mat_org = gt_tracks_all;
    %load('../data/crop_embryo_data_500x500x30x40/gt_tracks_tyxz_with_split.mat','gt_track_split_kids');
    %gt_mat_org = gt_track_split_kids';
    %gt_mat_org = gt_mat_org(:);
end


% let's first downsample the detection results
[h, w, z] = size(org_refine_res{1});
sc_f = 2; % we resize the data to [h/sc_f, w/sc_f, z, t]

st_loc = [];
sz_crop = [];
% st_loc = [251, 1, 1];
% sz_crop = [200, 250, z];
[refine_res_in, embryo_vid, gt_mat, threshold_res_in, varMaps, ...
    eigMaps] = data_scaling(sc_f, st_loc, ...
    sz_crop, org_refine_res, embryo_vid_org, gt_mat_org, ...
    org_threshold_res, org_varMap, org_eigMaps);

g = graphPara_cell(sum(cellfun(@(x) max(x(:)), refine_res_in)));%1

% reload gt_mat
g.gt_mat = gt_mat;

q = initial_q(sc_f, true);

% save the data for debug C++ project
% save_data_in_bin_file(refine_res_in, threshold_res_in, varMaps, eigMaps);

profile off;
profile on;
[movieInfo, movieInfoAll, out_refine_res, refine_resAll, threshold_resAll] = ...
    mcfTracking_cell(refine_res_in, embryo_vid, threshold_res_in, ...
    varMaps, eigMaps, g, q);
profile viewer;

%% save the iterative results
it = numel(movieInfoAll)-1;
[gt_tracks, gt_tracks_ov] = loc2id(gt_mat, refine_resAll{it});
[~, ratio_mat, track_id_head, ratio_mat2] = validate_res(movieInfoAll{it}, gt_tracks_ov);
bestOvRatio = nanmax(ratio_mat,[],2);
nanmean(bestOvRatio)
% save the tracks with ground truth
movieInfo_tmp = tree2tracks_cell(refine_resAll{it}, q, false);
movieInfo_tmp.tracks = cell(length(track_id_head), 1);
for i = 1:length(track_id_head)
    tr_id = track_id_head(i);
    if isnan(tr_id)
        movieInfo_tmp.tracks{i} = track_id_head(i);
    else
        movieInfo_tmp.tracks{i} = movieInfoAll{it}.tracks{tr_id};
    end
end

if abs(g.observationCost)>10
    video_folder = 'newSegV19_addMiss_24';
else
    video_folder = 'newSeg_10';
end
%video_folder = 'No_track_link';
%write_id(res{cnt,1}(write_id)<0.8) = [];
% generate the data with originally detected cells
p = track_disp_para();
%if ~exist('dataWithParticles', 'var')
dataWithParticles = generate3DImWithParticle(...
    cat(4, embryo_vid{:}), movieInfoAll{it}, p);
%end
write_id = [105,61,247,177,40,237,10,20,194, 7 222 22 181 315 266 307 236];
%[5 31 47 102 258 254 162];%find(res{cnt,1}<0.8); % consider no overlap
drawTracksStart1st(cat(4, embryo_vid{:}), fullfile(res_folder,disp_folder,video_folder),...
    movieInfoAll{it}, dataWithParticles, write_id);

video_folder = 'tmp_disp';
write_id = 1;
movieInfo_tmp.tracks{end+1} = movieInfo.tracks{write_id};
drawTracksStart1st(cat(4, embryo_vid{:}), fullfile(res_folder,disp_folder,video_folder),...
    movieInfo_tmp, dataWithParticles, numel(movieInfo_tmp.tracks));

% video_folder = 'way_2_5';
% if ~exist(fullfile(res_folder,disp_folder,video_folder),'dir')
%     mkdir(fullfile(res_folder,disp_folder,video_folder));
% end
% color_im = display_cl_frame(movieInfo, refine_res);
% color_im = cat(5, color_im{:});
% color_im = color_im*255;
% 
% drawTracksStart1st(cat(4, embryo_vid{:}), fullfile(res_folder,disp_folder,video_folder),...
%     movieInfo, color_im);

% tr_id = 35;
% track = movieInfo_tmp.tracks{35};
% coods = movieInfo_tmp.vox(track);
% centers = zeros(numel(coods), 3);
% for i=1:numel(coods)
%     centers(i,:) = mean(coods{i},1);
% end
% 
% tmp = centers(:,1);
% centers(:,1) = centers(:,2);
% centers(:,2) = tmp;
% centers = cat(2, [1:size(centers,1)]', centers);
% centers(:,2:4) = centers(:,2:4) ./ (ds_sz./[h,w,z]);