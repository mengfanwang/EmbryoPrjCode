function newIdMap = regionGrow_graphCut2d(idMap, eigMap,eigThres, vid, save_folder,over_seg_flag)
% region grow: shrink the 2d gaps detected by 2d principal curvature
% the idea is based on boyu's graph-cut (max-flow)
% INPUT:
% idMap: the label map of h*w*z
% eigMap: the principal curvature values of foreground in idMap
% eigThres: the threshold of principal curvature for eigMap
% vid: the orgriginal data
% save_folder: the folder to save images for sanity check
% over_seg_flag: true, oversegment the region and then merge, false
% otherwise
% OUTPUT:
% newIdMap: detected regions with each region one unique id
%
% contact: ccwang@vt.edu, 02/01/2020

if nargin < 5
    save_folder = [];
    over_seg_flag = false; % 
elseif nargin < 6
    over_seg_flag = false; %
end
newIdMap = zeros(size(idMap));
cost_design = [1 1];
connect = 8;

edgeConnect = 24;

p_thres = 0.01; % pvalue threshold

% deal with one component each time
s = regionprops3(idMap, {'VoxelIdxList'});
neiMap = zeros(3,3,3);
neiMap(:,:,2) = 1;
shift = [10 10 2];
regCnt = 0;

%% using different way to do gap test
% data stablization
vid = sqrt(vid);
if 1
    % use noise variance
    OrSt.noiseVar = calVarianceStablizationBY(vid, 0.8);
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
ok_id = [5,7:12,14,15,17:19,22:24,26,27,30,31,33:39];
z_merge = [9, 17, 20, 24, 26, 31, 36, 37, 39];
for i=4%z_merge%numel(s.VoxelIdxList)
    fprintf('%d out of %d\n', i, numel(s.VoxelIdxList));
    yxz = s.VoxelIdxList{i};
    [eigComp, linerInd, edge_flag] = crop3D(eigMap, yxz, shift);
    if edge_flag && size(vid,1) > 200% the region is on the edge
        continue;
    end
    vidComp = crop3D(vid, yxz, shift);
    idComp = crop3D(idMap, yxz, shift);
    idComp = idComp == i;

    %% way 1: over-segment regions first and try merge them
    if over_seg_flag
        eigCompMap = eigComp>eigThres; % gaps detected by principal curvature
        % !!! demean of principal curvature value
        eigComp = eigComp - eigThres;
%         newIdComp = shapeRefine(idComp, eigCompMap);
        newIdComp = idComp;
        newIdComp(eigCompMap) = 0; % if we segment current region further
        [L, n] = bwlabeln(newIdComp, neiMap);
        % for small region, do not grow them, merge them to principal
        % curvature maps
        v_l = regionprops3(L,'VoxelIdxList');
        com_len = cellfun(@length, v_l.VoxelIdxList);
        invalid_comps = find(com_len<=10);
        for j=1:length(invalid_comps)
            iv_px = v_l.VoxelIdxList{invalid_comps(j)};
            eigCompMap(iv_px) = 1;
            eigComp(iv_px) = 0;
            newIdComp(iv_px) = 0;
            L(iv_px) = 0;
        end
        if sum(com_len>10) == 1 % still one component
            regCnt = regCnt + 1;
            newIdMap(yxz) = regCnt;
            continue;
        end
        inLabel = L;
        
        loopCnt = 1;
        while true
            [newLabel, regTest] = gapTest(L, n,vidComp, idComp, eigComp,eigCompMap,...
                connect, cost_design, edgeConnect, p_thres, OrSt);
            loopCnt = loopCnt + 1;
            if loopCnt >3
                break;
            end
            n = max(newLabel(:));
            L = newLabel;
        end
        newLabel = region_refine(newLabel, n,  eigComp, connect, [1 1]);
    else
    %% way 2: simply try to segment regions using slice information
        newLabel = sliceRefine(vidComp, eigComp, idComp);
    end
    %% update the id map    
    [newLabel,n] = bwlabeln(newLabel>0);
    %n = max(newLabel(:));
    valid_newLabel = newLabel(:)>0;
    newIdMap(linerInd(valid_newLabel)) = newLabel(valid_newLabel) + regCnt;
    
    regCnt = regCnt + n;
    if true %sum(ok_id==i)==0 || (sum(ok_id==i)>0 && n > 1)
        %% write data into images for error check
        if ~isempty(save_folder) % write image data into disk
            vidComp = crop3D(vid, yxz, shift);
            vidComp = scale_image(vidComp, 0, 1);
            [h,w,z] = size(L);
            combined1 = zeros(h,w,3,z);
            combined = zeros(h,w,3,z);
            combined2 = zeros(h,w,3,z);
            for j=1:z
                combined(:,:,1,j) = 0.8*(regTest(:,:,j)>0);
                combined1(:,:,1,j) = 0.5*regTest(:,:,j)>0;
                combined1(:,:,2,j) = vidComp(:,:,j);
                combined(:,:,3,j) = 0.5*(idComp(:,:,j)>0);
                combined(:,:,2,j) = vidComp(:,:,j);
                %             combined(:,:,2,j) = L(:,:,j)>0;
                combined2(:,:,1,j) = 0.1*eigCompMap(:,:,j)>0;
                combined2(:,:,2,j) = vidComp(:,:,j);
            end
            %zzshow(combined)
            tifwrite(combined, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_comp']));
            %         tifwrite(combined1, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_gap']));
            %         tifwrite(vidComp, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_org']));
            tifwrite(combined2, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_pc']));
            tmprgb = label2rgb3d(newLabel,'jet', [0 0 0], 'shuffle');
            tifwrite(tmprgb, ...
                fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_newLabel']));
            tmprgb = label2rgb3d(inLabel,'jet', [0 0 0], 'shuffle');
            tifwrite(tmprgb, ...
                fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_inLabel']));
            tifwrite(scale_image(vidComp,0,1), ...
                fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_input']));
            
        end
    end
end

end