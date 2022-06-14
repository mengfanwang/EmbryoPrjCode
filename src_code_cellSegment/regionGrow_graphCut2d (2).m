function newIdMap = regionGrow_graphCut2d(idMap, eigMap,eigThres, vid, save_folder)
% region grow: shrink the gaps detected by principal curvature
% the idea is based on boyu's graph-cut (max-flow)
% INPUT:
% idMap: the label map of h*w*z
% eigMap: the principal curvature values of foreground in idMap
% eigThres: the threshold of principal curvature for eigMap
% OUTPUT:
% TODO:
%
% contact: ccwang@vt.edu

if nargin < 5
    save_folder = [];
end
newIdMap = zeros(size(idMap));
cost_design = 1;
connect = 8;
p_thres = 0.01; % pvalue threshold
app_str = '_pc2d_gc2d';
% deal with one component each time
s = regionprops3(idMap, {'VoxelIdxList'});
neiMap = zeros(3,3,3);
neiMap(:,:,2) = 1;
shift = 1;
regCnt = 0;

for i=1:50%:numel(s.VoxelIdxList)
    disp(i);
    yxz = s.VoxelIdxList{i};
    [eigComp, linerInd, edge_flag] = crop3D(eigMap, yxz, shift);
    if edge_flag % the region is on the edge
        continue;
    end
    idComp = crop3D(idMap, yxz, shift);
    idComp = idComp == i;
    eigCompMap = eigComp>eigThres; % gaps detected by principal curvature
    % !!! demean of principal curvature value
    eigComp = eigComp - eigThres;
    % if more than one component
    newIdComp = idComp;
    newIdComp(eigCompMap) = 0; % if we segment current region further
    [L, n] = bwlabeln(newIdComp, neiMap);
%     v_l = regionprops3(L,'VoxelIdxList');
%     com_len = cellfun(@length, v_l);
%     valid_comps = find(com_len>10);
    if n == 1 % still one component
        regCnt = regCnt + 1;
        newIdMap(yxz) = regCnt;
        continue;
    end
    newLabel = L;
    for j=1:n
        sMap = L==j;
        tMap = idComp<1 | (L~=j & L~=0);
        [dat_in, src, sink] = graphCut_build(eigComp,eigCompMap, sMap, tMap, connect, cost_design);
        G = digraph(dat_in(:,1),dat_in(:,2),dat_in(:,3));
        [~,~,cs,~] = maxflow(G, src, sink); % cs: fg, ct: bg
        cur_l = newLabel(cs(cs<numel(newLabel)));
        if length(unique(cur_l))>2
            keyboard;
        end
        newLabel(cs(cs<numel(newLabel))) = j;
    end
    % merge over-segmented regions
    vidComp = crop3D(vid, yxz, shift);
    inLabel = newLabel;
    [newLabel, regTest] = edgeTest(vidComp, inLabel, eigCompMap, connect, p_thres);

    if ~isempty(save_folder) % write image data into disk
        vidComp = crop3D(vid, yxz, shift);
        vidComp = scale_image(vidComp, 0, 1);
        [h,w,z] = size(L);
        combined1 = zeros(h,w,3,z);
        combined2 = zeros(h,w,3,z);
        for j=1:z
%             combined(:,:,1,j) = 0.8*(double(idComp(:,:,j)>0) - ...
%                 double(newLabel(:,:,j)>0));
            combined1(:,:,1,j) = 0.5*regTest(:,:,j)>0;
            combined1(:,:,2,j) = vidComp(:,:,j);
%             combined(:,:,3,j) = 0.5*(idComp(:,:,j)>0);
%             combined(:,:,1,j) = newLabel(:,:,j)>0;
%             combined(:,:,2,j) = L(:,:,j)>0;
            combined2(:,:,1,j) = 0.5*eigCompMap(:,:,j)>0;
            combined2(:,:,2,j) = vidComp(:,:,j);
        end
        %zzshow(combined)
        %tifwrite(combined, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_noDemean']));
        tifwrite(combined1, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_gap']));
        tifwrite(vidComp, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_org']));
        tifwrite(combined2, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_pc']));
        tmprgb = label2rgb3d(newLabel,'jet', [0 0 0], 'shuffle');
        tifwrite(tmprgb, ...
            fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_res']));
        tmprgb = label2rgb3d(inLabel,'jet', [0 0 0], 'shuffle');
        tifwrite(tmprgb, ...
            fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_input']));
    end
    valid_newLabel = newLabel(:)>0;
    newIdMap(linerInd(valid_newLabel)) = newLabel(valid_newLabel) + regCnt;
    
    regCnt = regCnt + j;
end

end