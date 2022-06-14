if isunix
    addpath('/home/ccw/Dropbox/cc_ImHandle/');
else
    addpath('C:\Users\Congchao\Dropbox\cc_ImHandle\');
end
addpath('src_code_matlab');
addpath('src_code_cellSegment');
addpath('src_code_cellTracker');


mydir  = 'C:\Users\Congchao\Google Drive\Reports\Phd_exams\Preliminary exam\test_data\CTChallenge_TRIC\';
addpath('./');
% data folder
data_folder = mydir;
tif_files = dir(fullfile(data_folder, 't*.tif'));
gt_file = tifread(fullfile(data_folder, 'man_seg_021_004.tif'));
fg = imdilate(gt_file>0, strel('disk', 30));
fg = repmat(fg, 1,1, 13);
% results folder
save_folder = mydir;
res_folder = fullfile(save_folder,'t021');
if ~exist(res_folder,'dir')
    mkdir(res_folder);
end

%% synQuant for the 40 frames
% add synQuant java path
Pij = fullfile('../../ParticleTracking', 'src_synquant/ij-1.52i.jar');
javaaddpath(Pij);
p1 = fullfile('../../ParticleTracking','src_synquant/commons-math3-3.6.1.jar');
javaaddpath(p1);%
p0 = fullfile('../../ParticleTracking', 'src_synquant/SynQuantVid_v1.2.4.jar');
javaaddpath(p0);%

z_mat = cell(numel(tif_files), 1);
id_mat = cell(numel(tif_files), 1);
fMaps = cell(numel(tif_files), 1);
q.minIntensity = 0;
for i=1:numel(tif_files)
    fprintf('processing %d/%d file\n', i, numel(tif_files));
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    org_im(~fg) = 0;
    % 3D version
    [zMap, synId, fMap] = Synquant4Embryo_Paramater(org_im, q, ...
        500, 5000, 0, 0.1, 2);
    
    z_mat{i} = single(zMap);
    id_mat{i} = uint16(synId);
    fMaps{i} = fMap;
end
save(fullfile(res_folder, 'synQuant_res.mat'), 'z_mat', 'id_mat','fMaps','-v7.3');
% remove java path
javarmpath(p0);
javarmpath(p1);
javarmpath(Pij);
%% refine results from synQuant
eig_res_2d = cell(numel(tif_files), 1);
eig_res_3d = cell(numel(tif_files), 1);
eig_overlay = cell(numel(tif_files), 1);
for i=1:numel(tif_files)
    fprintf('cal priCur %d/%d file\n', i, numel(tif_files));
    tic;
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    synId = id_mat{i};
    sigma = 3; % 2d principal curvature smooth scale is smaller
    [eig2d, ~] = principalCv2d(org_im, synId, sigma, fMaps{i});
    sigma = [5,1];
    [eig3d, overlay_cl] = principalCv3d(org_im, synId, sigma, fMaps{i});

    eig_res_2d{i} = single(eig2d);
    eig_res_3d{i} = single(eig3d);
    eig_overlay{i} = overlay_cl;
    toc;
end
save(fullfile(res_folder, 'synQuant_priCvt_res.mat'), 'eig_res_2d',...
    'eig_res_3d','eig_overlay','-v7.3');
%% calculate the variance map of all frames
varMap = cell(numel(tif_files), 1);
scale_term = 2000;
parfor i=1:numel(tif_files)
    disp(i);
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    vid = 255*org_im/scale_term;
    varMap{i} = cell(3,2);
    [varMap{i}{1,1}, varMap{i}{2,1},varMap{i}{3,1}] = ...
        calVarianceStablizationBY(vid, 0.8);
    vid_stb = sqrt(vid+3/8);
    [varMap{i}{1,2}, varMap{i}{2,2},varMap{i}{3,2}] = ...
        calVarianceStablizationBY(vid_stb, 0.8);
end
save(fullfile(res_folder, 'varianceMap.mat'), 'varMap','-v7.3');

%% region refine based on 4d information (infor across >1 frame)
load(fullfile(res_folder, 'synQuant_priCvt_res.mat'));
load(fullfile(res_folder, 'synQuant_res.mat'));
load(fullfile(res_folder, 'varianceMap.mat'));
scale_term = 2000;
refine_res = cell(numel(tif_files), 1);
threshold_res = cell(numel(tif_files), 1);
seg_res_folder = [];
parfor i=1:numel(tif_files)
    seg_res_folder = fullfile(res_folder,'segment_cells');
    seg_res_folder = fullfile(seg_res_folder,['frame_', num2str(i)]);
    if ~exist(seg_res_folder,'dir')
        mkdir(seg_res_folder);
    end
    org_im = cell(3,1);
    synId = cell(3,1);
    eigAll = cell(3,1);
    varMapAll = cell(3,1);
    for j=i-1:i+1
        if j>0 && j<=numel(tif_files)
            tmp_im = tifread(fullfile(tif_files(j).folder, tif_files(j).name));
            tmp_im = 255*tmp_im/scale_term;
            tmp_im(~fg) = 0;
            org_im{j-i+2} = tmp_im;
            
            synId{j-i+2} = id_mat{j};
            eigAll{j-i+2} = cell(2,1); % save both 2d and 3d pincipal curvature
            eigAll{j-i+2}{1} = eig_res_2d{i};
            eigAll{j-i+2}{2} = eig_res_3d{i};
            varMapAll{j-i+2} = varMap{j};
        end
    end
    % pre-frame
    tic;
%     [newIdMap, thresholdMap] = regionWiseAnalysis4d(synId, ...
%         eigAll,org_im, varMapAll, seg_res_folder);
    [newIdMap, thresholdMap] = regionWiseAnalysis4d(synId, ...
        eigAll,org_im, varMapAll, []);
    toc;
    
    refine_res{i} = uint32(newIdMap);
    threshold_res{i} = uint8(thresholdMap);
end
save(fullfile(res_folder, 'synQuant_refine_res_4d_v6.mat'), 'refine_res',...
    'threshold_res','-v7.3');
% v2: grow fg if foreground touching boundary in xy direction
% v3: grow fg if the detected cell touches boundary in xy direction
% v4: grow fg if the detected cell touches boundary in any 3 direction
% v5: grow fg in corresponding direction if the detected cell touches x, y
%     or z direction
% v6: shift = [10 10 2] ==> [20 20 4]

