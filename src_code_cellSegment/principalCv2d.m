function [eig_all, overlay_cl] = principalCv2d(vid, synId, sigma, fmap, save_folder)
% use principle curvature to segment connected cells
% INPUT:
% vid: original 3d image data
% synId: regions detected by synQuant
% sigma: smoothness scale, (sigma of gaussian, isotropic)
% save_folder: the folder to save results
% OUTPUT:
% eig_all: eigen value of all the foreground voxels (maximum eigen value)
% overlay_cl: overlapping the foreground with the detected valid gaps (gaps
% mean the voxels with positive eigen value)

% contact: ccwang@vt.edu, 02/01/2020
if nargin < 4
    fmap = [];
    save_folder = [];
end
if nargin < 5
    save_folder = [];
end

[h,w,z] = size(synId);
% get 2nd order derivative  % x direction is horiztional and y is vertical
Dxx = zeros(size(synId));
Dyy = zeros(size(synId));
Dxy = zeros(size(synId));
Dyx = zeros(size(synId));
sm_vid = zeros(h,w,z);
for i=1:z
    im = vid(:,:,i);
    im = imgaussfilt(im,sigma);
    sm_vid(:,:,i) = im;
    [lx, ly] = gradient(im);
    [lxx,lyx] = gradient(lx);
    [lxy, lyy] = gradient(ly);
    Dxx(:,:,i) = lxx;
    Dyy(:,:,i) = lyy;
    Dxy(:,:,i) = lxy;
    Dyx(:,:,i) = lyx;
end

% test each connected component
eig_all = zeros(size(synId));
eig1 = zeros(size(synId));
if isempty(fmap)
    fmap = imdilate(synId>0, strel('sphere', 5));
end
% fmap = scale_image(sm_vid,0,255)>50;%

fmap = ones(size(fmap)); % if we cal all the voxes
s = regionprops3(fmap, {'VoxelIdxList'});

for i=1:numel(s.VoxelIdxList)%[3 47 76]%
    %disp(i);
    vox = s.VoxelIdxList{i};
    xx = Dxx(vox); yy = Dyy(vox);
    xy = Dxy(vox); yx = Dyx(vox);
    
    C = zeros(numel(s.VoxelIdxList{i}),2);
    parfor j=1:numel(s.VoxelIdxList{i})
        MM = [xx(j) xy(j);yx(j) yy(j)];
        [~,Eval] = eig(MM);
        dEval = diag(Eval);
        c = sort(dEval,'descend');
        C(j,:) = c';
    end
    eig_all(vox) = C(:,1);
    if ~isempty(save_folder) && i<=50
        cur_id = find(synId(vox)>0,1);
        cur_id = synId(vox(cur_id));
        eig1(vox) = C(:,1);
        %save this region into single data in save_folder
        [yy,xx,zz] = ind2sub(size(synId),vox);
        ymin = max(1, min(yy)-3); ymax = min(h, max(yy)+3);
        xmin = max(1, min(xx)-3); xmax = min(w, max(xx)+3);
        zmin = max(1, min(zz)-3); zmax = min(z, max(zz)+3);
        cur_vid = vid(ymin:ymax, xmin:xmax, zmin:zmax);
        cur_synId = synId(ymin:ymax, xmin:xmax, zmin:zmax);
        cur_prCv = eig1(ymin:ymax, xmin:xmax, zmin:zmax);

        cur_vid = cur_vid ./ max(cur_vid(:));
        [c_h,c_w,c_z] = size(cur_vid);
        out_cell = zeros(c_h,c_w,3,c_z);
        for j=1:c_z
            out_cell(:,:,2,j) = cur_vid(:,:,j);
            out_cell(:,:,1,j) = 0.8*(cur_prCv(:,:,j)>0);
            out_cell(:,:,3,j) = 0.5*(cur_synId(:,:,j)>0);
        end
        tifwrite(out_cell, fullfile(save_folder, ['OneComp2dPrCv_',num2str(cur_id)]));
        eig1(vox) = 0;
    end
end

overlay_cl = uint8(zeros(h,w,3,z));
for i=1:z
    idMap = synId(:,:,i);
    idMap(idMap>0) = 125;
    b = labeloverlay(uint8(idMap), eig_all(:,:,i)>0);
    if size(b,3)==1
        b = cat(3,b,b,b);
    end
    overlay_cl(:,:,:,i) = b;
end
end