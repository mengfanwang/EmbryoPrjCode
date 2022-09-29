clc;clear;close all;
im = tifread('/work/Mengfan/Embryo/22-01-11/sameViewFusion_crop/1.tif');
load('/work/Mengfan/Embryo/22-01-11/sameViewDetection_crop/synQuant_refine_res_4d_v9.mat');
refine_res = refine_res{1};
load('/work/Mengfan/Embryo/22-01-11/sameViewDetection_crop/synQuant_res_2iter.mat');
id_map = id_mat_2nd{1};
id_list = unique(id_map(:,:,71:100));
id_list = id_list(2:end);
cell_list = [];
%%
miss_ind = getFN;
ss = 0;
s_list = [];
for ii = 1:size(miss_ind,1)
    if refine_res{1}(miss_ind(ii,2),miss_ind(ii,1),miss_ind(ii,3)) == 0
        ss = ss+1;
    end
end
disp(ss)
%%
miss_ind = getFP;
ss = 0;
s_list = [];
for ii = 1:size(miss_ind,1)
    if refine_res{1}(miss_ind(ii,2),miss_ind(ii,1),miss_ind(ii,3)) > 0
        ss = ss+1;
        s_list = [s_list ii];
    end
end
disp(ss)
%%
miss_ind = getMissingCell;
ss = 0;
s_list = [];
for ii = 1:size(miss_ind,1)
    if refine_res{1}(miss_ind(ii,2),miss_ind(ii,1),miss_ind(ii,3)) == 0
        ss = ss+1;
        s_list = [s_list ii];
    end
end
disp(ss)
%%
for ii = 1:size(miss_ind,1)
    cell_list = [cell_list id_map(miss_ind(ii,2),miss_ind(ii,1),miss_ind(ii,3))];
end
uncell_list = [];
for ii = 1:length(id_list)
    if ~ismember(id_list(ii), cell_list)
        uncell_list = [uncell_list id_list(ii)];
    end
end

s = regionprops3(id_map, {'VoxelIdxList','Volume','SurfaceArea', 'ConvexVolume', 'BoundingBox'});
%% volume
cell_volume = s.Volume(cell_list);
uncell_volume = s.Volume(uncell_list);
subplot(1,3,1);
histogram(cell_volume,0:50:600); hold on;
histogram(uncell_volume, 0:50:600);
title('Volume');

%% surface_volume ratio
cell_surface = s.SurfaceArea(cell_list);
uncell_surface = s.SurfaceArea(uncell_list);
subplot(1,3,2);
histogram(cell_surface./cell_volume,0:0.5:3.5); hold on;
histogram(uncell_surface./uncell_volume,0:0.5:3.5);
title('Surface volume ratio');

% %% fill ratio
% cell_convexVolume = s.ConvexVolume(cell_list);
% uncell_convexVolume= s.ConvexVolume(uncell_list);
% subplot(2,4,3);
% histogram(cell_volume./cell_convexVolume,0.2:0.1:1); hold on;
% histogram(uncell_volume./uncell_convexVolume,0.2:0.1:1);
% title('Fill ratio');
% 
% %% average intensity
% cell_intensity = arrayfun(@(x) mean(im(s.VoxelIdxList{x})), cell_list);
% uncell_intensity = arrayfun(@(x) mean(im(s.VoxelIdxList{x})), uncell_list);
% subplot(2,4,4);
% histogram(cell_intensity); hold on;
% histogram(uncell_intensity);
% title('Average intensity');

% %% constrast based on order-statistics
% cell_score = zeros(size(cell_list));
% uncell_score = zeros(size(uncell_list));
% tic;
% for ii = 1:length(cell_list)
%     cell_score(ii) = getZscore1(cell_list(ii), id_map, refine_res, im);
% end
% for ii = 1:length(uncell_list)
%     uncell_score(ii) = getZscore1(uncell_list(ii), id_map, refine_res, im);
% end
% toc
% subplot(1,3,3);
% histogram(cell_score); hold on;
% histogram(uncell_score);
% title('Z score');

%% faster way
cell_score = zeros(size(cell_list));
uncell_score = zeros(size(uncell_list));
tic;
radius = 5;
for ii = 1:length(cell_list)
    cell_id = cell_list(ii);
    ulf_x = max(ceil(s.BoundingBox(cell_id,2))-radius,1);
    ulf_y = max(ceil(s.BoundingBox(cell_id,1))-radius,1);
    ulf_z = max(ceil(s.BoundingBox(cell_id,3))-radius,1);
    bre_x = min(floor(s.BoundingBox(cell_id,2))+s.BoundingBox(cell_id,5)+radius,size(id_map,1));
    bre_y = min(floor(s.BoundingBox(cell_id,1))+s.BoundingBox(cell_id,4)+radius,size(id_map,2));
    bre_z = min(floor(s.BoundingBox(cell_id,3))+s.BoundingBox(cell_id,6)+radius,size(id_map,3));
    cell_score(ii) = getZscore2(cell_id, id_map(ulf_x:bre_x, ulf_y:bre_y, ulf_z:bre_z),...
        refine_res(ulf_x:bre_x, ulf_y:bre_y, ulf_z:bre_z), im(ulf_x:bre_x, ulf_y:bre_y, ulf_z:bre_z),radius);
end
for ii = 1:length(uncell_list)
    cell_id = uncell_list(ii);
    ulf_x = max(ceil(s.BoundingBox(cell_id,2))-radius,1);
    ulf_y = max(ceil(s.BoundingBox(cell_id,1))-radius,1);
    ulf_z = max(ceil(s.BoundingBox(cell_id,3))-radius,1);
    bre_x = min(floor(s.BoundingBox(cell_id,2))+s.BoundingBox(cell_id,5)+radius,size(id_map,1));
    bre_y = min(floor(s.BoundingBox(cell_id,1))+s.BoundingBox(cell_id,4)+radius,size(id_map,2));
    bre_z = min(floor(s.BoundingBox(cell_id,3))+s.BoundingBox(cell_id,6)+radius,size(id_map,3));
    uncell_score(ii) = getZscore2(cell_id, id_map(ulf_x:bre_x, ulf_y:bre_y, ulf_z:bre_z),...
        refine_res(ulf_x:bre_x, ulf_y:bre_y, ulf_z:bre_z), im(ulf_x:bre_x, ulf_y:bre_y, ulf_z:bre_z),radius);
end
toc
subplot(1,3,3);
histogram(cell_score); hold on;
histogram(uncell_score);
title('Z score');


function z_score = getZscore1(cell_id, id_map, refine_res, im)
    bnd_nei_radius = 5;
    sph = strel('sphere', bnd_nei_radius);
    se2 = strel(sph.Neighborhood(:,:,bnd_nei_radius+1));
    bnd_nei_radius = 3;
    sph = strel('sphere', bnd_nei_radius);
    se1 = strel(sph.Neighborhood(:,:,bnd_nei_radius+1));

    fg_locs = id_map == cell_id;
    bg_locs = imdilate(fg_locs,se2) & ~imdilate(fg_locs,se1)  & (refine_res == 0) & (id_map == 0);
    bg = mean(im(bg_locs));
    std = sqrt(bg);
    fg = mean(im(fg_locs));
    [mu, sigma] = ordStatApproxKsecWith0s_mat(fg/std, bg/std, []);
    z_score = (fg/std - bg/std - mu)/sigma;
end

function z_score = getZscore2(cell_id, id_map, refine_res, im, radius)
    bnd_nei_radius = radius;
    sph = strel('sphere', bnd_nei_radius);
    se2 = strel(sph.Neighborhood(:,:,bnd_nei_radius+1));
    bnd_nei_radius = ceil(radius/2);
    sph = strel('sphere', bnd_nei_radius);
    se1 = strel(sph.Neighborhood(:,:,bnd_nei_radius+1));

    fg_locs = id_map == cell_id;
    bg_locs = imdilate(fg_locs,se2) & ~imdilate(fg_locs,se1)  & (refine_res == 0) & (id_map == 0);
    bg = mean(im(bg_locs));
    std = sqrt(bg);
    fg = mean(im(fg_locs));
    [mu, sigma] = ordStatApproxKsecWith0s_mat(fg/std, bg/std, []);
    z_score = (fg/std - bg/std - mu)/sigma;
end

function miss_ind = getMissingCell
    roi = ReadImageJROI('/work/Mengfan/Embryo/22-01-11/sameViewDetection_crop/0097-0308-0399.roi');
    miss_ind = round([roi.mfCoordinates roi.vnSlices]);
end

function miss_ind = getFP
    roi = ReadImageJROI('/work/Mengfan/Embryo/22-01-11/sameViewDetection_crop2_test/sameViewDetection_crop2_test_inten_50/FP1.roi');
    miss_ind = round([roi.mfCoordinates roi.vnSlices]);
    roi = ReadImageJROI('/work/Mengfan/Embryo/22-01-11/sameViewDetection_crop2_test/sameViewDetection_crop2_test_inten_50/FP2.roi');
    miss_ind = [miss_ind; round([roi.mfCoordinates roi.vnSlices])];
    roi = ReadImageJROI('/work/Mengfan/Embryo/22-01-11/sameViewDetection_crop2_test/sameViewDetection_crop2_test_inten_50/FP3.roi');
    miss_ind = [miss_ind; round([roi.mfCoordinates roi.vnSlices])];
    miss_ind(miss_ind(:,3) > 40 | miss_ind(:,3) < 21,:) = [];
    b = load('/work/Mengfan/Embryo/22-01-11/sameViewDetection_crop2_test/sameViewDetection_crop2_test_inten_50/synQuant_refine_res_4d_v9.mat','refine_res');
    a = sub2ind_direct(size(b.refine_res{1}),miss_ind(:,2),miss_ind(:,1),miss_ind(:,3));
    a_new = zeros(size(a));
    a_list = [];
    for ii = 1:length(a)
        if ~ismember(b.refine_res{1}(a(ii)), a_list)
            a_list = [a_list b.refine_res{1}(a(ii))];
            a_new(ii) = 1;
        end
    end
    a = a(logical(a_new));
    miss_ind = zeros(length(a), 3);
    [miss_ind(:,2), miss_ind(:,1), miss_ind(:,3)] = ind2sub_direct(size(b.refine_res{1}), a);
end

function miss_ind = getFN
    roi = ReadImageJROI('/work/Mengfan/Embryo/22-01-11/sameViewDetection_crop3_test/FN.roi');
    miss_ind = round([roi.mfCoordinates roi.vnSlices]);
    miss_ind(miss_ind(:,3) > 40 | miss_ind(:,3) < 21,:) = [];
    b = load('/work/Mengfan/Embryo/22-01-11/sameViewDetection_crop3_test/synQuant_remain_FN.mat','b');
    a = sub2ind_direct(size(b.b),miss_ind(:,2),miss_ind(:,1),miss_ind(:,3));
    a_new = zeros(size(a));
    a_list = [];
    for ii = 1:length(a)
        if ~ismember(b.b(a(ii)), a_list)
            a_list = [a_list b.b(a(ii))];
            a_new(ii) = 1;
        end
    end
    a = a(logical(a_new));
    miss_ind = zeros(length(a), 3);
    [miss_ind(:,2), miss_ind(:,1), miss_ind(:,3)] = ind2sub_direct(size(b.b), a);
    roi = ReadImageJROI('/work/Mengfan/Embryo/22-01-11/sameViewDetection_crop2_test/sameViewDetection_crop2_test_inten_50/FN.roi');
    miss_ind2 = round([roi.mfCoordinates roi.vnSlices]);
    miss_ind2(miss_ind(:,3) > 40 | miss_ind(:,3) < 21,:) = [];
    miss_ind = [miss_ind; miss_ind2];
end