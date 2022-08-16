clc;clear;close all;
dbstop if error
addpath('../');
addpath('../src_code_matlab');
addpath('../src_code_cellSegment');
addpath('../src_code_cellTracker');
addpath('../src_code_visualization');

%% system and path
if isunix
    addpath('/home/mengfan/ForExecute/Tools/MatlabTools');
    addpath('/home/mengfan/ForExecute/cc_ImHandle');
    data_folder = '/work/Mengfan/Embryo/22-01-11/fusion_0.5/';
    res_folder = '/work/Mengfan/Embryo/22-01-11/detection_0.5/';
else
    addpath('D:\Congchao''s code\cc_ImHandle\');
    addpath D:\MatlabTools;
    data_folder  = 'E:\Embryo\TM0-49\debug_v2\input\';
    res_folder = fullfile('E:\Embryo\TM0-49\debug_v2\');
end

tif_files = dir(fullfile(data_folder, '/*.tif'));
if ~exist(res_folder,'dir')
    mkdir(res_folder);
end
minIntensity = 50; % The middle of two Gaussian intensity distributions (
                    % should learn from data)

%% synQuant
tic;
% add synQuant java path
Pij = fullfile('../src_synquant/ij-1.52i.jar');
javaaddpath(Pij);
p1 = fullfile('../src_synquant/commons-math3-3.6.1.jar');
javaaddpath(p1);
p0 = fullfile('../src_synquant/SynQuantVid_v1.2.5.1.jar');
javaaddpath(p0);

z_mat = cell(numel(tif_files), 1);
id_mat = cell(numel(tif_files), 1);
fMaps = cell(numel(tif_files), 1);
q.minIntensity = minIntensity;
for i=1:numel(tif_files)
    fprintf('processing %d/%d file\n', i, numel(tif_files));
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    [h, w, slices] = size(org_im);
    %out_ims = SliceImage(in_im);
    %org_im = imresize3(org_im,round([h/2 w/2 slices/2]));
    sigma = 0.8;
    sm_im = imgaussfilt3(org_im,sigma);
    %q.posEigMap = eig_res_3d{i}>0;
    % 3D version
    [zMap, synId, fMap] = Synquant4Embryo_Paramater(sm_im, q);
    
    z_mat{i} = single(zMap);
    id_mat{i} = uint16(synId);
    fMaps{i} = fMap;
    toc
end
save(fullfile(res_folder, 'synQuant_res.mat'), 'z_mat', 'id_mat','fMaps','-v7.3');
% tifwrite(uint8(id_mat{1}*255), [res_folder 'result_1']);
% remove java path
javarmpath(p0);
javarmpath(p1);
javarmpath(Pij);
fprintf('Synquant part running time:') % around 12000s 0.25
toc

%% refine results from synQuant
tic;
load(fullfile(res_folder, 'synQuant_res.mat'));
eig_res_2d = cell(numel(tif_files), 1);
eig_res_3d = cell(numel(tif_files), 1);
eig_overlay = cell(numel(tif_files), 1);
for i=1:numel(tif_files)
    fprintf('cal priCvt %d/%d file\n', i, numel(tif_files));
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    synId = id_mat{i};
    
    fMaps = ones(size(org_im));
    
    [eig2d, ~] = principalCv2d(org_im, synId, 1, fMaps);   % 1: sigma
    sigma = [1 1 1];
    %   grad3d = imgradient3(imgaussfilt3(org_im,sigma));
    [eig3d, overlay_cl] = principalCv3d(org_im, synId, sigma, fMaps);
    
    % use eig_res_2d to store grad3d to keep strcut uniform
    eig_res_2d{i} = single(eig2d);
    eig_res_3d{i} = single(eig3d);
    eig_overlay{i} = [];
end
save(fullfile(res_folder, 'synQuant_priCvt_res.mat'), 'eig_res_2d',...
    'eig_res_3d','eig_overlay','-v7.3');
fprintf('Principal curvature running time:'); % around 5700s 0.25
toc
%% calculate the variance map of all frames
tic;
varMap = cell(numel(tif_files), 1);
scale_term = 300;
for i=1:numel(tif_files)
    disp(i);
    
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    vid = 255*org_im/scale_term;
    varMap{i} = cell(3,2);
    [varMap{i}{1,1}, varMap{i}{2,1},varMap{i}{3,1}] = ...
        calVarianceStablizationBY(vid, 0.8, 3);
    vid_stb = sqrt(vid+3/8);
    [varMap{i}{1,2}, varMap{i}{2,2},varMap{i}{3,2}] = ...
        calVarianceStablizationBY(vid_stb, 0.8, 3);
end
save(fullfile(res_folder, 'varianceMap.mat'), 'varMap','-v7.3');
fprintf('Variance running time:'); % around 7500s 0.25
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % variance map by ConvexVST from Mengfan
% addpath D:\MatlabTools\
% varMap = cell(numel(tif_files), 1);
% scale_term = 300;
% for i=1:numel(tif_files)
%     disp(i);
%     tic;
%     org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
%     org_im(org_im<0) = 0;
%     org_im(org_im>300) = 300;
% %     vid = 255*org_im/scale_term;
%     varMap{i} = cell(3,2);
%     
%     options.histEdges = 0.5:299.5;
%     options.binEdges = -0.5:300.5;
%     options.sampleSize = 200;
%     options.ratio = 0.03;
%     options.display = false;
%     [his,variance,parameters] = histogramCount(org_im,[3 3 3],options);
%     variance = variance*255^2/scale_term^2;
%     varMap{i}{2,1} = median(variance);
%     varMap{i}{3,1} = variance;
%     varMap{i}{1,1} = interp1(parameters.histCenters,variance,org_im);
% 
%     [stabilizeFunction, variance] = convexOptimization(his,parameters);
%     stabilizeFunction(1) = stabilizeFunction(2);
%     stabilizeFunction(end) = stabilizeFunction(end-1);
%     stabilizeFunction = 300*(stabilizeFunction - min(stabilizeFunction))/(max(stabilizeFunction) - min(stabilizeFunction));
%     stab_im = interp1(0:300,stabilizeFunction,org_im);
%     [his,variance,parameters] = histogramCount(stab_im,[3 3 3],options);
%     variance = variance*255^2/scale_term^2;
%     varMap{i}{2,2} = median(variance);
%     varMap{i}{3,2} = variance;
%     varMap{i}{1,2} = interp1(parameters.histCenters,variance,stab_im);
% 
%     toc;
% end
% save(fullfile(res_folder, 'varianceMap.mat'), 'varMap','-v7.3');

%% region refine based on 4d information (infor across >1 frame) around 1 hour
tic;
scale_term = 300;
load(fullfile(res_folder, 'synQuant_priCvt_res.mat'));
load(fullfile(res_folder, 'synQuant_res.mat'));
load(fullfile(res_folder, 'varianceMap.mat'));
refine_res = cell(numel(tif_files), 1);
threshold_res = cell(numel(tif_files), 1);
multi_frames_flag = false; % use multiple frames for segmentation
cell_wise_save_flag = false; % save each cell segment
for i=1:numel(tif_files)
    fprintf('Processing the frame %d ', i);
    
    if cell_wise_save_flag
        seg_res_folder = fullfile(res_folder,'segment_cells');
        seg_res_folder = fullfile(seg_res_folder,['frame_', num2str(i)]);
        if ~exist(seg_res_folder,'dir')
            mkdir(seg_res_folder);
        end
    end
    if multi_frames_flag % use 3 consecutive frames
        org_im = cell(3,1);
        synId = cell(3,1);
        eigAll = cell(3,1);
        varMapAll = cell(3,1);
        for j=i-1:i+1
            if j>0 && j<=numel(tif_files)
                tmp_im = tifread(fullfile(tif_files(j).folder, tif_files(j).name));
                tmp_im = 255*tmp_im/scale_term;
                org_im{j-i+2} = tmp_im;
                synId{j-i+2} = id_mat{j};
                eigAll{j-i+2} = cell(2,1); % save both 2d and 3d pincipal curvature
                eigAll{j-i+2}{1} = eig_res_2d{i};
                eigAll{j-i+2}{2} = eig_res_3d{i};
                varMapAll{j-i+2} = varMap{j};
            end
        end
    else
        tmp_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
%         org_im = 255*tmp_im/scale_term;
        org_im = tmp_im;
        synId = id_mat{i};
        eigAll = cell(2,1); % save both 2d and 3d pincipal curvature
        eigAll{1} = eig_res_2d{i};
        eigAll{2} = eig_res_3d{i};
        varMapAll = varMap{i};
    end
    
    %profile on;
    [newIdMap, thresholdMap] = m_regionWiseAnalysis4d(synId, ...
            eigAll,org_im, varMapAll, []);%, i
    %profile viewer;
    %profile off;
    refine_res{i} = uint32(newIdMap);
    threshold_res{i} = uint8(thresholdMap);
    toc
end
save(fullfile(res_folder, 'synQuant_refine_res_4d_v9.mat'), 'refine_res',...
    'threshold_res','-v7.3');
fprintf('Refinement running time:'); % around 183000s 0.25
toc
% tifwrite(uint8((refine_res{1}>0)*255), [res_folder 'result_2']);
%% write 3d segmentation result with label
mkdir(fullfile(res_folder,'synQuant_refine_res_4d_v9'));
for tt = 1:numel(tif_files)
    tt_ind = num2str(100000+tt);
    tt_ind = tt_ind(2:6);
    org_im_all = tifread(fullfile(tif_files(tt).folder, tif_files(tt).name));
    refine_res_all = refine_res{tt};
    labelwrite(uint8(org_im_all),refine_res_all, fullfile(res_folder,'synQuant_refine_res_4d_v9',tt_ind));
end
%% write 4d segmentatin result with label
% org_im_all = zeros(h, w, slices, numel(tif_files));
% refine_res_all = zeros(size(org_im_all));
% for tt = 1:numel(tif_files)
%         org_im_all(:,:,:,tt) = tifread(fullfile(tif_files(tt).folder, tif_files(tt).name));
%         refine_res_all(:,:,:,tt) = refine_res{tt};
% end
% labelwrite(uint8(org_im_all),refine_res_all, [res_folder '\synQuant_refine_res_4d_v9']);