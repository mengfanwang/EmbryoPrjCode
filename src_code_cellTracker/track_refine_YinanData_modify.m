% %This is the main function for tracking results refinement
% if isunix
%     addpath('/home/ccw/Dropbox/cc_ImHandle/');
%     save_folder = '/home/ccw/Desktop/embryo_res_folder';
% else
%     addpath('C:\Users\Congchao\Dropbox\cc_ImHandle\');
%     save_folder = 'C:\Users\Congchao\Desktop\cell_detection_samples';
% end

% addpath('dt');
% addpath(genpath('CINDA/'));
% addpath('src_code_matlab');
% addpath('src_code_cellTracker');
% addpath('src_code_cellSegment');
% addpath('debug_func/');
% addpath('../../ParticleTracking/src_code/');
% addpath('../../ParticleTracking/src_uTrack/');

clc;clear;close all;
dbstop if error
addpath('D:\Congchao''s code\cc_ImHandle\');
addpath('..\dt');
addpath(genpath('..\CINDA\'));
addpath('..\debug_func');
addpath('..\src_code_matlab');
addpath('..\src_code_cellSegment');
addpath('..\src_code_cellTracker');
addpath('..\src_code_visualization');
addpath('D:\Congchao''s code\ParticleTracking\src_code');
addpath('D:\Congchao''s code\ParticleTracking\src_uTrack');

% set saving folder
save_folder = 'E:\Embryo\TM0-49\track_v1';
data_folder = 'E:\Embryo\TM0-49\data_crop_reg';
tif_files = dir(fullfile(data_folder, '*.tif'));
res_folder = 'E:\Embryo\TM0-49\track_v1';
if ~exist(res_folder,'dir')
    mkdir(res_folder);
end
% disp_folder = 'synQuant_refine_res';

% load original data and detection results and ground truth
if ~exist('refine_res','var') % segmentation results
    %load(fullfile(res_folder, 'synQuant_refine_res.mat'),'refine_res');
    load(fullfile(res_folder, 'synQuant_refine_res_reg.mat'),...
        'refine_reg', 'threshold_reg');
    refine_res = refine_reg; threshold_res = threshold_reg;
    for ii = 1:length(refine_res) % rearrange id for cropped data
        refine_res{ii} = rearrange_id(refine_res{ii});
    end
    org_refine_res = refine_res;
    org_threshold_res = threshold_res;
end
scale_term = 300;
embryo_vid_org = cell(numel(tif_files), 1); % original data
for i=1:numel(tif_files)
    embryo_vid_org{i} = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    embryo_vid_org{i} = 255*embryo_vid_org{i}./scale_term;
end

if ~exist('varMap','var') % segmentation results
    %load(fullfile(res_folder, 'synQuant_refine_res.mat'),'refine_res');
%     load(fullfile(res_folder, 'varianceMap.mat'),'varMap');
%     load(fullfile(res_folder, 'synQuant_priCvt_res.mat'),...
%         'eig_res_2d', 'eig_res_3d');
    load([res_folder '\synQuant_pri_var_reg.mat']);
    eig_res_2d = eig_res_2d_reg;
    eig_res_3d = eig_res_3d_reg;
    varMap = varMap_reg;
    
    org_varMap = varMap;    
    org_eigMaps = cell(numel(eig_res_2d),1);
    for i=1:numel(org_eigMaps)
        org_eigMaps{i} = cell(2,1);
        org_eigMaps{i}{1} = eig_res_2d{i};
        org_eigMaps{i}{2} = eig_res_3d{i};
    end
end



% let's first downsample the detection results
[h, w, z] = size(org_refine_res{1});
sc_f = 1; % we resize the data to [h/sc_f, w/sc_f, z, t]

st_loc = [];
sz_crop = [];
% st_loc = [251, 1, 1];
% sz_crop = [200, 250, z];
gt_mat_org = {};
[refine_res_in, embryo_vid, gt_mat, threshold_res_in, varMaps, ...
    eigMaps] = data_scaling(sc_f, st_loc, ...
    sz_crop, org_refine_res, embryo_vid_org, gt_mat_org, ...
    org_threshold_res, org_varMap, org_eigMaps);

g = graphPara_cell(sum(cellfun(@(x) max(x(:)), refine_res_in)));%1


q = initial_q(sc_f, true);

profile off;
profile on;
[movieInfo, movieInfoAll, out_refine_res, refine_resAll,...
    threshold_res, threshold_resAll] = ...
    mcfTracking_cell(refine_res_in, embryo_vid, threshold_res_in, ...
    varMaps, eigMaps, g, q);
profile viewer;

% display_rgb_time_course(out_refine_res);

if q.saveInterMediateRes
    save(fullfile(save_folder, 'movieInfo.mat'), 'movieInfo',...
        'movieInfoAll','-v7.3');
    save(fullfile(save_folder, 'refine_res.mat'), 'out_refine_res',...
        'refine_resAll','-v7.3');
    save(fullfile(save_folder, 'threshold_res.mat'), threshold_res,...
        'threshold_resAll','-v7.3');
else
    save(fullfile(save_folder, 'movieInfo.mat'), 'movieInfo', '-v7.3');
    save(fullfile(save_folder, 'track_refine_res.mat'), 'out_refine_res','-v7.3');
    save(fullfile(save_folder, 'track_threshold_res.mat'), 'threshold_res','-v7.3');
end

return;
% save seg data for Yinan
seg_res_dir = fullfile(save_folder, 'seg_before_refine');
for i = 1:numel(tif_files)
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    save_seg_masks(org_im/2000, refine_res{i}, fullfile(seg_res_dir, ...
        tif_files(i).name(1:length(tif_files(i).name)-4)))
end

% save track data for Yinan
track_res_dir = fullfile(save_folder, 'track_res');
l = cellfun(@length, movieInfo.tracks);
[~, idx] = sort(l, 'descend');
save_cnt = 1;
for i = 1:numel(movieInfo.tracks)
    track_id = idx(i);
    if movieInfo.frames(movieInfo.tracks{track_id}(1)) == 1
        saved = save_maxProj_track(movieInfo, i, out_refine_res,...
            embryo_vid, track_res_dir, save_cnt);
        if saved
            save_cnt = save_cnt + 1;
        end
    end
end