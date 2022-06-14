%% display the refine results with tracking purely based on overlapping 
disp_folder = 'synQuant_refine_res';
max_cell_num = max(cellfun(@(x) max(x(:)), refine_res));
cmap = jet(double(max_cell_num));
cmap = cmap(randperm(max_cell_num),:);

pre_ref_cmap = [];
remain_ref_cmap = cmap;
parents = [];
pre_syn_cmap = [];
remain_syn_cmap = cmap;

pre_syn_idmap = [];
pre_ref_idmap = [];
iou_thres = 1e-5;
for i=1:numel(tif_files)
    org_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    synId = id_mat{i};
    synId_ref = refine_res{i};
    
    %% overlayed synQuant results
    res_name = ['detResOv_',tif_files(i).name(1:end-4)];
    scaled_im = scale_image(org_im,0,2);
    [h,w,z] = size(scaled_im);
    im_cl = zeros(h,w,3,z);
    for j=1:z
        im_cl(:,:,1,j) = 0.5*(synId_ref(:,:,j)>0);
        im_cl(:,:,2,j) = scaled_im(:,:,j);
        im_cl(:,:,3,j) = 0.5*(synId(:,:,j)>0);
    end
    tifwrite(im_cl, fullfile(res_folder, disp_folder, res_name));
    
    %% refined results
    cell_ids = 1:max(synId_ref(:));
    [parents, ious] = ov2link(pre_ref_idmap, synId_ref, cell_ids);
    parents(ious<iou_thres) = 0;
    [pre_ref_cmap, remain_ref_cmap] = cl2parent(pre_ref_cmap, ...
        remain_ref_cmap, parents);

    res_name = ['Res_synRef_',tif_files(i).name(1:end-4)];
    tifwrite(label2rgb3d(synId_ref, pre_ref_cmap ,[0 0 0]), ...
        fullfile(res_folder, disp_folder, res_name));
     %% synQuant results
    cell_ids = 1:max(synId(:));
    [parents, ious] = ov2link(pre_syn_idmap, synId, cell_ids);
    parents(ious<iou_thres) = 0;
    [pre_syn_cmap, remain_syn_cmap] = cl2parent(pre_syn_cmap, ...
        remain_syn_cmap, parents);

    res_name = ['Res_syn_',tif_files(i).name(1:end-4)];
    tifwrite(label2rgb3d(synId, pre_syn_cmap ,[0 0 0]), ...
        fullfile(res_folder, disp_folder, res_name));

     %% principal curvature results
    res_name = ['Res_prCv_',tif_files(i).name(1:end-4)];
    tifwrite(eig_overlay{i}, fullfile(res_folder,disp_folder,res_name));
    
    pre_syn_idmap = synId;
    pre_ref_idmap = synId_ref;
end


%% tracking way 1: use distance the same as we tracked tips of processes
% build inputs
dets = cell(numel(tif_files),1);
for i=1:numel(tif_files)
    s = regionprops3(refine_res{i}, 'Centroid');
    dets{i} = [s.Centroid(:,2) s.Centroid(:,1) s.Centroid(:,3)];
end
movieInfo = tree2tracks(dets, false);% false means we do not use existing tracking results
xCoord = cell(numel(dets),1);
yCoord = cell(numel(dets),1);
zCoord = cell(numel(dets),1);
for i=1:numel(dets)
    yCoord{i} = dets{i}(:,1);
    xCoord{i} = dets{i}(:,2);
    zCoord{i} = dets{i}(:,3);
end
% main function
[movieInfo,movieInfoAll] = mcfTracking(movieInfo, xCoord,yCoord,zCoord);


% display the results one by one
video_folder = 'way_1';
movieInfo.xCoord = movieInfo.orgCoord(:,1);
movieInfo.yCoord = movieInfo.orgCoord(:,2);
movieInfo.zCoord = movieInfo.orgCoord(:,3);

embryo_vid = cell(numel(tif_files), 1);
for i=1:numel(tif_files)
    embryo_vid{i} = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
end
embryo_vid = cat(4, embryo_vid{:});
drawTracksStart1st(embryo_vid, dets,fullfile(res_folder,disp_folder,video_folder),...
    movieInfo);
% display the results frame by frame
file_names = cell(numel(tif_files),1);
for i=1:numel(tif_files)
    file_names{i} = ['Res_syn_',tif_files(i).name(1:end-4)];
end
display_cl_frame(movieInfo, refine_res, fullfile(res_folder, disp_folder, video_folder),...
    file_names);


%% summary statistics
mot = cell(numel(movieInfo.tracks),1);
for i=1:numel(movieInfo.tracks)
    cur_track = movieInfo.tracks{i};
    if length(cur_track) <=5
        continue;
    end
    mot{i} = movieInfo.orgCoord(cur_track(2:end),:) - movieInfo.orgCoord(cur_track(1:end-1),:);
end
mot_all=cat(1,mot{:});

figure;histogram(mot_all(:,1), 50);

dist_all = sqrt(mot_all(:,1).^2+mot_all(:,2).^2+mot_all(:,3).^2);

figure;histogram(dist_all, 50);