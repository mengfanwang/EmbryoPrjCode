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
%     data_folder = '/work/Mengfan/Embryo/23-11-01_Yinan/view9';
%     res_folder = '/work/Mengfan/Embryo/23-11-01_Yinan/Detection_view9';

    data_folder = '/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/view9';
    res_folder = '/work/Mengfan/Embryo/20220518 isl2b H2Bmcherry overnight/Detection_view9';
    
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
minIntensity = 50; 

%% synQuant
tic;
% add synQuant java path
Pij = fullfile('../src_synquant/ij-1.52i.jar');
javaaddpath(Pij);
p1 = fullfile('../src_synquant/commons-math3-3.6.1.jar');
javaaddpath(p1);
p0 = fullfile('../src_synquant/SynQuantVid_v1.2.5.1.jar');
javaaddpath(p0);

if ~isfolder(fullfile(res_folder, 'synQuant_res'))
    mkdir(fullfile(res_folder, 'synQuant_res'));
end
if ~isfolder(fullfile(res_folder, 'synQuant_res_tif'))
    mkdir(fullfile(res_folder, 'synQuant_res_tif'));
end
q.minIntensity = minIntensity;
ds_scale = 1; % down sample scale
for i=1:20 %numel(tif_files)
    fprintf('processing %d/%d file\n', i, numel(tif_files));
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    [~, org_name, ~] = fileparts(tif_files(i).name);
    [h, w, slices] = size(org_im);
    org_im = imresize3(org_im,round([h/ds_scale w/ds_scale slices]));
    org_im = org_im - 200;  
    org_im(org_im < 0) = 0;
    
    sigma = [3 3 1];
    sm_im = imgaussfilt3(org_im,sigma);
    [zMap, synId, fMap] = m_Synquant4Embryo_Paramater(sm_im, q);
        
    z_mat = single(zMap);
    id_mat = uint16(synId);
    fMaps = fMap;

    z_mat = imresize3(z_mat,[h w slices],'nearest');
    id_mat = imresize3(id_mat,[h w slices],'nearest');
    fMaps = imresize3(fMaps,[h w slices],'nearest');
    toc
    save(fullfile(res_folder, 'synQuant_res', [org_name '.mat']), 'z_mat', 'id_mat','fMaps','-v7.3');
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    org_im = org_im - 200;  
    org_im(org_im < 0) = 0;
    labelwrite(uint8(org_im/2), id_mat, fullfile(res_folder, 'synQuant_res_tif', org_name));
end

% 
% labelwrite(org_im, id_mat{1}, fullfile(res_folder, 'synQuant_res'));
% remove java path
javarmpath(p0);
javarmpath(p1);
javarmpath(Pij);
fprintf('Synquant part running time:') % around 12000s 0.25
toc

%% refine results from synQuant
tic;
if ~isfolder(fullfile(res_folder, 'synQuant_priCvt_res'))
    mkdir(fullfile(res_folder, 'synQuant_priCvt_res'));
end
for i=154:193 % 1:numel(tif_files)
    fprintf('cal priCvt %d/%d file\n', i, numel(tif_files));
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    [~, org_name, ~] = fileparts(tif_files(i).name);
    load(fullfile(res_folder, 'synQuant_res', [org_name '.mat']));

    synId = id_mat;
    fMaps = ones(size(org_im));
    
    % start of video: 8/[8 8 2]  % end of video: 4/[4 4 1]
    % current setting: sigma = 4*(2-tt/250)
    sigma = 4*(2-str2num(org_name)/250);
    [eig2d, ~] = principalCv2d(org_im, synId, sigma, fMaps);
    sigma = [sigma sigma sigma/4];
    %   grad3d = imgradient3(imgaussfilt3(org_im,sigma));
    [eig3d, overlay_cl] = principalCv3d(org_im, synId, sigma, fMaps);
    
    % use eig_res_2d to store grad3d to keep strcut uniform
    eig_res_2d = single(eig2d);
    eig_res_3d = single(eig3d);
    save(fullfile(res_folder, 'synQuant_priCvt_res', [org_name '.mat']), 'eig_res_2d',...
    'eig_res_3d','-v7.3');
end
fprintf('Principal curvature running time:'); % around 5700s 0.25
toc
%% calculate the variance map of all frames
tic;
if ~isfolder(fullfile(res_folder, 'varianceMap'))
    mkdir(fullfile(res_folder, 'varianceMap'));
end
scale_term = 500;
for i= 1:numel(tif_files)  
    disp(i);
    
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    [~, org_name, ~] = fileparts(tif_files(i).name);
    org_im = org_im - 200;  
    org_im(org_im < 0) = 0;
    vid = 255*org_im/scale_term;
    varMap = cell(3,2);
    [varMap{1,1}, varMap{2,1},varMap{3,1}] = ...
        calVarianceStablizationBY(vid, 0.8, 3);
    vid_stb = sqrt(vid+3/8);
    [varMap{1,2}, varMap{2,2},varMap{3,2}] = ...
        calVarianceStablizationBY(vid_stb, 0.8, 3);
    save(fullfile(res_folder, 'varianceMap', [org_name '.mat']), 'varMap','-v7.3');
end
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
if ~isfolder(fullfile(res_folder, 'synQuant_refine_res'))
    mkdir(fullfile(res_folder, 'synQuant_refine_res'));
end
if ~isfolder(fullfile(res_folder, 'synQuant_refine_res_tif'))
    mkdir(fullfile(res_folder, 'synQuant_refine_res_tif'));
end
scale_term = 300;
multi_frames_flag = false; % use multiple frames for segmentation
cell_wise_save_flag = false; % save each cell segment
for i=1:numel(tif_files)
    fprintf('Processing the frame %d ', i);
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    [~, org_name, ~] = fileparts(tif_files(i).name);
    load(fullfile(res_folder, 'synQuant_res', [org_name '.mat']));
    load(fullfile(res_folder, 'synQuant_priCvt_res', [org_name '.mat']));
    load(fullfile(res_folder, 'varianceMap', [org_name '.mat']));

    synId = id_mat;
    eigAll = cell(2,1); % save both 2d and 3d pincipal curvature
    eigAll{1} = eig_res_2d;
    eigAll{2} = eig_res_3d;
    varMapAll = varMap;

    %profile on;
%     [newIdMap, thresholdMap] = m_regionWiseAnalysis4d(synId, ...
%             eigAll,org_im, varMapAll, []);%, i
    [newIdMap, thresholdMap] = m_regionWiseAnalysis4d_parallel(synId, ...
        eigAll,org_im, varMapAll, []);%, i
    %profile viewer;
    %profile off;
    refine_res = uint32(newIdMap);
    threshold_res = uint8(thresholdMap);
    toc
    save(fullfile(res_folder, 'synQuant_refine_res', [org_name '.mat']), 'refine_res',...
    'threshold_res','-v7.3');
    org_im(org_im>255) = 255;
    labelwrite(uint8(org_im), refine_res, fullfile(res_folder, 'synQuant_refine_res_tif', org_name));
end
fprintf('Refinement running time:'); % around 183000s 0.25
toc
% tifwrite(uint8((refine_res{1}>0)*255), [res_folder 'result_2']);