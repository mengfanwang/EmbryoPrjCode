function [eig_all, overlay_cl] = principalCv3d(vid, synId, sigma, fmap, save_folder)
% use principle curvature to segment connected cells

% get the hessian matrix
[Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, sm_vid] = Hessian3D(vid,sigma);
if nargin < 4
    fmap = [];
    save_folder = [];
end
if nargin < 5
    save_folder = [];
end
% if(sigma>1)
%     % Correct for scaling
%     c=(sigma^2);
%     Dxx = c*Dxx; Dxy = c*Dxy;
%     Dxz = c*Dxz; Dyy = c*Dyy;
%     Dyz = c*Dyz; Dzz = c*Dzz;
% end
[h,w,z] = size(synId);
% test each connected component
eig_all = zeros(size(synId));
eig1 = zeros(size(synId));

dir_y = zeros(size(synId));
dir_x = zeros(size(synId));
dir_z = zeros(size(synId));
% eig2 = zeros(size(synId));
% eig3 = zeros(size(synId));
valid_eig = 0.5;
%s = regionprops3(synId, {'VoxelIdxList'});
%fmap = imdilate(synId>0, strel('sphere', 5));
% if isempty(fmap)
%     fmap = imdilate(synId>0, strel('sphere', 5));
% end

fmap = ones(size(fmap)); % if we cal all the voxes
s = regionprops3(fmap, {'VoxelIdxList'});
for i=1:numel(s.VoxelIdxList)%[3 47 76]%
    %disp(i);
    vox = s.VoxelIdxList{i};
    xx = Dxx(vox); yy = Dyy(vox); zz = Dzz(vox);
    xy = Dxy(vox); xz = Dxz(vox); yz = Dyz(vox);
    
    C = zeros(numel(s.VoxelIdxList{i}),3);
    dir_xyz = zeros(numel(s.VoxelIdxList{i}),3);
    parfor j=1:numel(s.VoxelIdxList{i})
        MM = [xx(j), xy(j), xz(j);...
            xy(j), yy(j), yz(j);...
            xz(j), yz(j), zz(j)];
        [Evec,Eval] = eig(MM);
        dEval = diag(Eval);
        %[~,od] = sort(abs(dEval),'descend');
        %C(j,:) = dEval(od)';
        [c,od] = sort(dEval,'descend');
        C(j,:) = c';
        dir_xyz(j,:) = Evec(:, od(1))';
    end
    dir_x(vox) = dir_xyz(:,1);
    dir_y(vox) = dir_xyz(:,2);
    dir_z(vox) = dir_xyz(:,3);
    eig_all(vox) = C(:,1);
%     if ~isempty(save_folder) && i <= 50
%         cur_id = find(synId(vox)>0,1);
%         cur_id = synId(vox(cur_id));
%         eig1(vox) = C(:,1);
%     %     eig2(vox) = C(:,2);
%     %     eig3(vox) = C(:,3);
%         %save this region into single data in save_folder
%         [yy,xx,zz] = ind2sub(size(synId),vox);
%         ymin = max(1, min(yy)-3); ymax = min(h, max(yy)+3);
%         xmin = max(1, min(xx)-3); xmax = min(w, max(xx)+3);
%         zmin = max(1, min(zz)-3); zmax = min(z, max(zz)+3);
%         cur_vid = vid(ymin:ymax, xmin:xmax, zmin:zmax);
%         cur_synId = synId(ymin:ymax, xmin:xmax, zmin:zmax);
%         cur_prCv = eig1(ymin:ymax, xmin:xmax, zmin:zmax);
%         cur_vid = cur_vid ./ max(cur_vid(:));
%         [c_h,c_w,c_z] = size(cur_vid);
%         out_cell = zeros(c_h,c_w,3,c_z);
%         
% %         [c_x,c_y] = meshgrid(1:c_w,1:c_h);
%         for j=1:c_z
%             out_cell(:,:,2,j) = cur_vid(:,:,j);
%             %out_cell(:,:,1,j) = 0.8*(bwskel(cur_prCv(:,:,j)>valid_eig));
%             out_cell(:,:,1,j) = 0.8*(cur_prCv(:,:,j)>valid_eig);
%             out_cell(:,:,3,j) = 0.5*(cur_synId(:,:,j)>0);
%             
% 
% %             u11 = dir_x(ymin:ymax,xmin:xmax,zmin-1+j);
% %             u11(cur_prCv(:,:,j)<=valid_eig) = 0;
% %             v11 = dir_y(ymin:ymax,xmin:xmax,zmin-1+j);
% %             v11(cur_prCv(:,:,j)<=valid_eig) = 0;
% %             w11 = dir_z(ymin:ymax,xmin:xmax,zmin-1+j);     
% %             w11(cur_prCv(:,:,j)<=valid_eig) = 0;
% %             
% %             figure;
% %             imshow(out_cell(:,:,:,j));
% %             hold on
% %             quiver(c_x, c_y,v11,u11); 
% %             hold off;
% 
%         end
%         tifwrite(out_cell, fullfile(save_folder, ['OneComp3d_skel_PrCv_',num2str(cur_id)]));
%         eig1(vox) = 0;
%     end
end

overlay_cl = uint8(zeros(h,w,3,z));
% for i=1:z
%     idMap = synId(:,:,i);
%     idMap(idMap>0) = 125;
%     b = labeloverlay(uint8(idMap), eig_all(:,:,i)>valid_eig);
%     if size(b,3)==1
%         b = cat(3,b,b,b);
%     end
%     overlay_cl(:,:,:,i) = b;
% end
end