function [refine_res, embryo_vid, gt_mat, threshold_res, varMap, ...
    eigMaps] = data_scaling(sc_f, st_loc, sz_crop, org_refine_res, ...
    embryo_vid_org, gt_mat_org, org_threshold_res,...
    org_varMap, org_eigMaps)
% scaling the data and/or crop a sub region
if ~isempty(st_loc)
    h = sz_crop(1);
    w = sz_crop(2);
    z = sz_crop(3);
    end_loc = st_loc + sz_crop - 1;
else
    [h, w, z] = size(org_refine_res{1});
end
ds_sz = round([h/sc_f, w/sc_f, z]);

% st_loc(1:2) = round(st_loc(1:2)/sc_f);
% sz_crop(1:2) = round(sz_crop(1:2)/sc_f);


% embryo_vid = cell(numel(org_refine_res), 1);
% refine_res = cell(numel(org_refine_res), 1);
% threshold_res = cell(numel(org_refine_res), 1);
% eig_res_2d = cell(numel(org_refine_res), 1);
% eig_res_3d = cell(numel(org_refine_res), 1);
embryo_vid = embryo_vid_org;
refine_res = org_refine_res;
threshold_res = org_threshold_res;
eigMaps = org_eigMaps;

varMap = org_varMap;
for i = 1:numel(org_refine_res)
    if ~isempty(st_loc)
        refine_res{i} = refine_res{i}(st_loc(1):end_loc(1), ...
            st_loc(2):end_loc(2),st_loc(3):end_loc(3));
        refine_res{i} = rearrange_id(refine_res{i});
        embryo_vid{i} = embryo_vid{i}(st_loc(1):end_loc(1), ...
            st_loc(2):end_loc(2),st_loc(3):end_loc(3));
        threshold_res{i} = threshold_res{i}(st_loc(1):end_loc(1), ...
            st_loc(2):end_loc(2),st_loc(3):end_loc(3));
        
        eigMaps{i}{1} = eigMaps{i}{1}(st_loc(1):end_loc(1), ...
            st_loc(2):end_loc(2),st_loc(3):end_loc(3));
        eigMaps{i}{1} = eigMaps{i}{2}(st_loc(1):end_loc(1), ...
            st_loc(2):end_loc(2),st_loc(3):end_loc(3));
        
        varMap{i}{1,1} = varMap{i}{1,1}(st_loc(1):end_loc(1), ...
            st_loc(2):end_loc(2),st_loc(3):end_loc(3));
        varMap{i}{1,2} = varMap{i}{1,2}(st_loc(1):end_loc(1), ...
            st_loc(2):end_loc(2),st_loc(3):end_loc(3));
        %     else
        %         refine_res{i} = imresize3(org_refine_res{i},ds_sz,'method','nearest');
        %         embryo_vid{i} = imresize3(embryo_vid_org{i},ds_sz,'method','nearest');
        %         threshold_res{i} = imresize3(org_threshold_res{i},ds_sz,'method','nearest');
        %         eig_res_2d{i} = imresize3(org_eig_res_2d{i},ds_sz,'method','cubic');
        %         eig_res_3d{i} = imresize3(org_eig_res_3d{i},ds_sz,'method','cubic');
        %
        %         varMap{i}{1,1} = imresize3(varMap{i}{1,1},ds_sz,'method','cubic');
        %         varMap{i}{1,2} = imresize3(varMap{i}{1,2},ds_sz,'method','cubic');
    end
    refine_res{i} = imresize3(refine_res{i},ds_sz,'method','nearest');
    embryo_vid{i} = imresize3(embryo_vid{i},ds_sz,'method','nearest');
    threshold_res{i} = imresize3(threshold_res{i},ds_sz,'method','nearest');
    
    eigMaps{i}{1}(isnan(eigMaps{i}{1})) = 0;
    eigMaps{i}{2}(isnan(eigMaps{i}{2})) = 0;
    eigMaps{i}{1} = imresize3(eigMaps{i}{1},ds_sz,'method','cubic');
    eigMaps{i}{2} = imresize3(eigMaps{i}{2},ds_sz,'method','cubic');
    
    varMap{i}{1,1}(isnan(varMap{i}{1,1})) = 0;
    varMap{i}{1,2}(isnan(varMap{i}{1,2})) = 0;
    varMap{i}{1,1} = imresize3(varMap{i}{1,1},ds_sz,'method','cubic');
    varMap{i}{1,2} = imresize3(varMap{i}{1,2},ds_sz,'method','cubic');

end
gt_mat = gt_mat_org;
for i=1:numel(gt_mat)
    if ~isempty(st_loc)
        y_v = gt_mat{i}(:,2)>=st_loc(1) & gt_mat{i}(:,2)<=end_loc(1);
        x_v = gt_mat{i}(:,3)>=st_loc(2) & gt_mat{i}(:,3)<=end_loc(2);
        z_v = gt_mat{i}(:,4)>=st_loc(3) & gt_mat{i}(:,4)<=end_loc(3);
        v_ = y_v & x_v & z_v;
        [lmap, n] = bwlabeln(v_);
        if n==0
            gt_mat{i} = [];
            continue;
        end
        s = regionprops(lmap,'Area');
        [~, od] = max([s.Area]);
        gt_mat{i} = gt_mat{i}(lmap==od,:);
%         loc = find(v_==0,1);
%         if ~isempty(loc)
%             if loc==1
%                 gt_mat{i} = [];
%                 continue;
%             else
%                 gt_mat{i} = gt_mat{i}(1:loc-1,:);
%             end
%         end
        gt_mat{i}(:,2) = gt_mat{i}(:,2) - st_loc(1) + 1;
        gt_mat{i}(:,3) = gt_mat{i}(:,3) - st_loc(2) + 1;
        gt_mat{i}(:,4) = gt_mat{i}(:,4) - st_loc(3) + 1;
    end
    
    gt_mat{i}(:,2:4) = gt_mat{i}(:,2:4).*(ds_sz./[h,w,z]);

end
le = cellfun(@length, gt_mat);
gt_mat = gt_mat(le>0);
end