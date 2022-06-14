function newIdMap = regionWiseAnalysis3d(idMap, eigMap,vid, fmap, save_folder)
% region grow: shrink the 3d gaps detected by 3d principal curvature
% the idea is based on boyu's graph-cut (max-flow)
% INPUT:
% idMap: the label map of h*w*z
% eigMap: the principal curvature values of foreground in idMap include 2d
% and 3d (cells with two elements)
% eigThres: the threshold of principal curvature for eigMap
% vid: the orgriginal data
% save_folder: the folder to save images for sanity check
% over_seg_flag: true, oversegment the region and then merge, false
% otherwise
% OUTPUT:
% newIdMap: detected regions with each region one unique id
%
% contact: ccwang@vt.edu, 02/07/2020
if nargin == 4
    save_folder = [];
end
cost_design = [1, 2];
connect = 6;

edgeConnect = 48; % three circle of neighbors

p_thres = 0.01; % pvalue threshold

% deal with one component each time
s = regionprops3(idMap, {'VoxelIdxList'});
neiMap = 26;
% neiMap = zeros(3,3,3);
% neiMap(:,:,2) = 1;
% neiMap(2,2,:) = 1;

shift = [10 10 2];
%% using different way to do gap test
% data stablization
vid_stb = sqrt(vid);
if 1
    % use noise variance
    OrSt.noiseVar = calVarianceStablizationBY(vid_stb, 0.8);
    OrSt.NoStbVar = calVarianceStablizationBY(vid, 0.8);
    muSigma = paraP3D(0.05, 0,0,255,4, 500, 0.1, 10);
    OrSt.NoStbMu = muSigma.mu;
    OrSt.NoStbSig = muSigma.sigma;
else
    % use variance of regions
    OrSt.noiseVar = [];
end
if 0
    % way 1: ttest
    OrSt = [];
elseif 0
    % way 2: truncated Gaussian
    muSigma = paraP3D(0.05, 0,0,255,4, 500, 0.1, 10);
    OrSt.mu = muSigma.mu;
    OrSt.sigma = muSigma.sigma;
else
    % way 3: order statistics
    OrSt.mu = [];
    OrSt.sigma = [];
end
%% get the gradient map
% F = vid * 0;
% for i=1:size(F,3)
%     F(:,:,i) = imgaussian(vid(:,:,i).^2, 2);
% end
% dx = gradient3(F,'x');
% dy = gradient3(F,'y');
% %dz = gradient3(F,'z');
% dMap = sqrt(dx.^2 + dy.^2);
eig2d = eigMap{1};
eig3d = eigMap{2};
loc_cells = cell(numel(s.VoxelIdxList),1);
minIn = sort(vid(fmap),'ascend');
minIn = minIn(round(length(minIn)*.001));
for i=79%[64 60 59 53 51]%1:numel(s.VoxelIdxList)%[4 23 62 67 78 79 81 84 87]%
    yxz = s.VoxelIdxList{i};

    fprintf('process %d out of %d\n', i, numel(s.VoxelIdxList));
    [vidStbComp, linerInd] = crop3D(vid_stb, yxz, shift);
    idComp = crop3D(idMap, yxz, shift);
    regComp = idComp == i;
    fmapComp = crop3D(fmap, yxz, shift);
    % build 3d/2d principla curveture score map
    eigComp = crop3D(eig3d, yxz, shift);
    eigCompMap = eigComp>0; % gaps detected by principal curvature
    eigComp(eigComp<0) = 0;
    score3dMap = scale_image(eigComp, 1e-3,1);
    eig2dComp = crop3D(eig2d, yxz, shift);
    eig2dComp(eig2dComp<0) = 0;
    score2dMap = scale_image(eig2dComp, 1e-3,1);
    %% grow the region firstly immediately after synQuant
%     regComp_rg = double(regComp);
%     regComp_rg(idComp ~= i & idComp > 0) = 1.5;
%     regComp_rg = regionGrow(regComp_rg, score3dMap, ...
%         fmapComp, connect, cost_design);

    %% over-segment regions first and try merge them
    newIdComp = regComp;
    newIdComp(eigCompMap) = 0; % if we segment current region further
    [L, n] = bwlabeln(newIdComp, neiMap);
    S = regionprops3(L, 'Volume');
    inval_ids = find([S.Volume] <= 10);
    if ~isempty(inval_ids)
        invalidMap = ismember(L, inval_ids);
        L(invalidMap) = 0;
        L = rearrange_id(L);
    end
    % for small region, do not grow them, merge them to principal
    % curvature maps
    v_l = regionprops3(L,'VoxelIdxList');
    if numel(v_l.VoxelIdxList) > 1
        %% for such region, we analyse them one by one
        loopCnt = 1;
        while true
            [newLabel, regTest] = gapTest3d(L, i, vidStbComp, idComp, ...
                score2dMap,fmapComp,...
                4, cost_design, edgeConnect, p_thres, OrSt);
            loopCnt = loopCnt + 1;
            if loopCnt > 1
                break;
            end
            n = max(newLabel(:));
            L = newLabel;
        end
    else
        L = regComp;
        inLabel = double(regComp);
        regTest = regComp*0;
        newLabel = inLabel;
    end
    %% update the id map
    scoreMap = score3dMap;%(score3dMap + score2dMap)/2;
    vidComp = crop3D(vid, yxz, shift);
    [newLabel] = region_refine(newLabel, L, ...
        vidComp, idComp, i, fmapComp, scoreMap,...
        OrSt, connect, cost_design);
    [newLabel, ~] = region_sanity_check(newLabel, 20);
    valid_newLabel = newLabel(:)>0;
    loc_cells{i} = [linerInd(valid_newLabel), newLabel(valid_newLabel)];
    %% write data into images for error check
    %inLabel = seedRegion; % seed label
    if ~isempty(save_folder) % write image data into disk
        %         score3dMap = scale_image(score3dMap, 0, 1);
        vidComp = scale_image(vidComp, 0, 1);
        [h,w,z] = size(regTest);
        maxId = max(newLabel(:));
        if maxId==0
            newLabel = regComp;
            maxId = 1;
            for ll = 1:maxId
                combined = zeros(h,w,3,z);
                for j=1:z
                    combined(:,:,1,j) = 0.5*(newLabel(:,:,j)==ll); % detected region
                    combined(:,:,2,j) = vidComp(:,:,j);
                end
                tifwrite(combined, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_', num2str(ll),'_comp']));
            end
        else
            for ll = 1:maxId
                combined = zeros(h,w,3,z);
                for j=1:z
                    combined(:,:,3,j) = 0.5*(newLabel(:,:,j)==ll); % detected region
                    combined(:,:,2,j) = vidComp(:,:,j);
                end
                tifwrite(combined, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_', num2str(ll),'_comp']));
            end
        end
%         for ll = 1:maxId
%             %             combined1 = zeros(h,w,3,z);
%             combined = zeros(h,w,3,z);
%             %             combined2 = zeros(h,w,3,z);
%             %             gradCombined = zeros(h,w,3,z);
%             for j=1:z
%                 combined(:,:,1,j) = 0.8*(regTest(:,:,j)>0);
%                 combined(:,:,3,j) = 0.5*(newLabel(:,:,j)==ll); % detected region
%                 combined(:,:,2,j) = vidComp(:,:,j);
%                 %                 combined1(:,:,1,j) = 0.5*regTest(:,:,j)>0; % detected gap
%                 %                 combined1(:,:,2,j) = vidStbComp(:,:,j);
%                 %                 %             combined(:,:,2,j) = L(:,:,j)>0;
%                 %                 combined2(:,:,1,j) = 0.1*eigCompMap(:,:,j)>0; % principal cur
%                 %                 combined2(:,:,2,j) = vidStbComp(:,:,j);
%                 %
%                 %                 gradCombined(:,:,2,j) = vidStbComp(:,:,j);
%                 %                 gradCombined(:,:,3,j) = 0.5*(inLabel(:,:,j)>0);
%             end
%             %zzshow(combined)
%             tifwrite(combined, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_', num2str(ll),'_comp']));
%             %         tifwrite(combined1, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_gap']));
%             %         tifwrite(vidComp, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_org']));
% 
%             %         tifwrite(combined2, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_pc']));
%             %         tmprgb = label2rgb3d(newLabel,'jet', [0 0 0], 'shuffle');
%             %         tifwrite(tmprgb, ...
%             %             fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_newLabel']));
%             %         tmprgb = label2rgb3d(inLabel,'jet', [0 0 0], 'shuffle');
%             %         tifwrite(tmprgb, ...
%             %             fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_inLabel']));
%             %         tifwrite(scale_image(vidComp,0,1), ...
%             %             fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_input']));
%             %         tifwrite(gradCombined, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_grad']));
%         end
    end
end
% label the new ID Map
newIdMap = zeros(size(idMap));
regCnt = 0;
for i=1:numel(loc_cells)
    cur_locs = loc_cells{i}(:,1);
    cur_labels = loc_cells{i}(:,2);
    if ~isempty(cur_locs)
        newIdMap(cur_locs) = cur_labels + regCnt;
        regCnt = regCnt + max(cur_labels);
    end
end

end