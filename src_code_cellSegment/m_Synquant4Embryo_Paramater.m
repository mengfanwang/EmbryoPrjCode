function [zMap, synId, fMap] = Synquant4Embryo_Paramater(imgIn, q, minSz, maxSz, flag_2d, minfill, maxWHRatio)
% cell detection using synQuant for 3D image
% INPUT:
% vid: is YXZ 3D image
% q: minIntensity
% minSz: the minimum size of valid region (default 100 for 3d, 20 for 2d)
% maxSz: the maximum size of valid region
% flag_2d: true means we process slice one by one, default is false
% OUTPUT:
% zMap: the zscore map of detected regions
% synId: the id map of the detection regions
% fMap: foreground maps
% contact: ccwant@vt.edu 02/04/2020
if nargin < 5 % the input args number is less than 5
    flag_2d = false;
end
if nargin < 3
    if flag_2d % TODO: 2d funciton has not been implemented
        minSz = 20;
        maxSz = 450;
        minfill = 0.2;
        maxWHRatio = 10;
    else
        minSz = 100;
        %minSz = 10; % indeed 5 should be OK, since 0.5^5<0.05
        maxSz = 3000;
        minfill = 0.0001;
        maxWHRatio = 100;
    end
end
fMap =  imgIn > q.minIntensity;
if isfield(q, 'posEigMap')
    posEigMap = q.posEigMap;
end
maxVid = max(imgIn(:));
if maxVid>255
    imgIn = round(255 * imgIn / maxVid);
end
q = paraQ3D(1,0,0.8); % number of channels (always 1); way to combine multiple channel (always 0);ratio used for noise estimation, small value==>small variance
% INPUTS of paraP3D
% 1. fdr threshold; 2. z-score threshold; 3. lowest intensity; 4. highest
% intensity; 5. minimum size; 6. largest size; 7. min fill; 8. max WHRatio
p = paraP3D(0, 0,0,255,minSz, maxSz, minfill, maxWHRatio);
vox_x = 2e-7; % useless indeed
% detection for one channel
[H1,W1,D1] = size(imgIn);
datx = zeros(D1,H1*W1,'uint8');
for ii=1:D1
    tmp = imgIn(:,:,ii)';
    datx(ii,:) = tmp(:);
end
det_res = ppsd3D(datx, W1, H1, vox_x, p,q);
zMap1 = det_res.ppsd_main.zMap;
synId1 = det_res.ppsd_main.kMap;

zMap = zeros(size(imgIn));
synId = zeros(size(imgIn));
for i=1:size(zMap,3)
    zMap(:,:,i) = zMap1(i,:,:);
    synId(:,:,i) = synId1(i,:,:);
end

synId(zMap<=1 | ~fMap) = 0;
synId = rearrange_id(synId);
s = regionprops3(synId, {'VoxelIdxList'});
cnt = 0;
synId = zeros(size(imgIn));
cell_level = zeros(numel(s.VoxelIdxList), 2);
for i=1:numel(s.VoxelIdxList)
    if length(s.VoxelIdxList{i}) >= minSz
        cnt = cnt + 1;
        synId(s.VoxelIdxList{i}) = cnt;
        cell_level(cnt,1) = mean(imgIn(s.VoxelIdxList{i}));
        cell_level(cnt,2) = length(s.VoxelIdxList{i});
    end
end
cell_level = cell_level(1:cnt,:);
%% get the cell information
% od = round(0.05*cnt);
% st_levels = sort(cell_level(:,1), 'ascend');
% min_level = st_levels(od);
% st_szs = sort(cell_level(:,2), 'ascend');
% minSz = st_szs(od);
% imgIn_org = imgIn;
%% mask out all the detected cells and re-do everything
% imgIn = imgIn_org;
% mask_out = posEigMap | imdilate(synId>0, strel('sphere',3));%
% imgIn(mask_out) = 0;
% 
% [H1,W1,D1] = size(imgIn);
% datx = zeros(D1,H1*W1,'uint8');
% useN = zeros(D1,H1*W1,'uint8');
% for ii=1:D1
%     tmp = imgIn(:,:,ii)';
%     datx(ii,:) = tmp(:);
%     tmp = synId(:,:,ii)';
%     useN(ii,:) = tmp(:)*2;
% end
% useN = useN(:);
% ppsd_main = fastppsdcore3D(datx, W1, H1, D1, p,q, useN);
% zMap1 = ppsd_main.zMap;
% synId1 = ppsd_main.kMap;
% 
% zMap_append = zeros(size(imgIn));
% synId_append = zeros(size(imgIn));
% for i=1:size(zMap,3)
%     zMap_append(:,:,i) = zMap1(i,:,:);
%     synId_append(:,:,i) = synId1(i,:,:);
% end
% 
% synId_append(zMap_append<=1 | ~fMap) = 0;
% synId_append = rearrange_id(synId_append);
% s = regionprops3(synId_append, {'VoxelIdxList','Centroid'});
% %cnt = 0;
% synId2 = zeros(size(imgIn));
% cell_level = cat(1, cell_level, zeros(numel(s.VoxelIdxList), 2));
% for i=1:numel(s.VoxelIdxList)
%     tmp_level = mean(imgIn(s.VoxelIdxList{i}));
%     center_xyz = s.Centroid(i,:);
%     %center_xyz = center_xyz([2 1 3]);
%     center_ind = get_index_from_center(center_xyz([2 1 3]), imgIn);
%     valid_center = ismember(center_ind,s.VoxelIdxList{i}); % center should inside region
% 
%     if sum(valid_center)==length(valid_center) && ...
%             length(s.VoxelIdxList{i}) >= minSz && tmp_level > min_level
%         cnt = cnt + 1;
%         synId2(s.VoxelIdxList{i}) = cnt;
%         cell_level(cnt,1) = tmp_level;
%         cell_level(cnt,2) = length(s.VoxelIdxList{i});
%     end
% end
% cell_level = cell_level(1:cnt,:);
% 
% synId_out = synId + synId2;
% % display
% overlay_im = zeros([H1,W1,3,D1], 'uint8');
% for i=1:D1
%     overlay_im(:,:,1,i) = uint8((synId(:,:,i)>0) * 100);
%     overlay_im(:,:,2,i) = uint8(imgIn_org(:,:,i)*2);
%     overlay_im(:,:,3,i) = uint8((synId2(:,:,i)>0) * 100);
% end
% svf = 'C:\Users\Congchao\Desktop\cell_detection_samples\crop_embryo_data_500x500x30x40\segmentation_from_synQuant\';
% files = dir(fullfile(svf,'*.tif'));
% tifwrite(overlay_im, fullfile(svf, ['overlay_synQuant_', num2str(numel(files)+1)]));
end


