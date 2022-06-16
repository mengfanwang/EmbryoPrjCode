clc;clear;close all;
dbstop if error
if isunix
    addpath('../src_code_cellSegment');
else
    addpath('..\src_code_cellSegment');
end


%% simulation
len = 1000;
tic;
for ii = 1:51
    smoothness = (ii - 1)/10;
    fprintf('Working on the smothness %f', smoothness);
    bound_size = 2*ceil(2*smoothness)+1;
    vid = normrnd(0, 1, [len+2*bound_size,len+2*bound_size, len+2*bound_size]);
    [eig_all, ~] = principalCv3d(vid, zeros(size(vid)), smoothness, ones(size(vid)));
    eig_all = eig_all(bound_size+1:end-bound_size, bound_size+1:end-bound_size, bound_size+1:end-bound_size);
    
    max_eig = max(eig_all(:));
    min_eig = min(eig_all(:));
    binLength = (max_eig-min_eig)/9999;
    h = hist(eig_all(:), min_eig:binLength:max_eig);
    
    PriCvtSim(ii).smoothness = smoothness;
    PriCvtSim(ii).max_eig = max_eig;
    PriCvtSim(ii).min_eig = min_eig;
    PriCvtSim(ii).binCenters = min_eig:binLength:max_eig;
    PriCvtSim(ii).hist = h;
    clear eig_all vid
    toc
end