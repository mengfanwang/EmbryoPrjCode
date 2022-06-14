picked_track_id = [45];

it = numel(movieInfoAll);
% movieInfo = movieInfoAll{it};%
refine_res = refine_resAll{it};
zzshow(label2rgb3d(refine_res{1}))


%% selected_ids read from gt_label.xlsx
gt_id_map = zeros(size(refine_res{1}));
cnt = 0;
used_ids = [];
for i=1:max(selected_ids(:,1))
    locxyz = selected_ids(selected_ids(:,1)==i,2:4);
    if isempty(locxyz)
        continue;
    end
    cur_id_map = zeros(size(refine_res{1}));
    used_flag = false;
    for j=1:size(locxyz,1)
        cur_id = refine_res{1}(locxyz(j, 2), locxyz(j, 1), locxyz(j, 3));
        if cur_id == 0
            disp(i);
        end
        if sum(used_ids == cur_id)>0
            used_flag = true;
            break;
        end
        used_ids = cat(1, used_ids, cur_id);
        cur_id_map = cur_id_map | refine_res{1}==cur_id;
    end
    if used_flag
        continue;
    end
    if sum(i==[12 14 43 44 45])
        cur_id_map = imerode(cur_id_map, strel('disk',2));
    end
    if sum(i==[60 61 70 71])
        cur_id_map = imdilate(cur_id_map, strel('disk',1));
    end
    if i==78
        cur_id_map = imdilate(cur_id_map, strel('disk',2));
    end
    if i==24
        cur_id_map(:,:,1:2) = 0;
    end
    if cur_id==refine_res{1}(153,168,24)
        cur_id_map1 = cur_id_map;
        cur_id_map1(1:157,:,:) = 0;
        cur_id_map(158:end,:,:) = 0;
        %keyboard;
        cnt = cnt + 1;
        cur_id_map = imopen(cur_id_map, strel('disk',1));
        gt_id_map(cur_id_map) = cnt;
        
        cnt = cnt + 1;
        cur_id_map1 = imclose(cur_id_map1, strel('disk',1));
        gt_id_map(cur_id_map1) = cnt;
    else
        cnt = cnt + 1;
        cur_id_map = imopen(cur_id_map, strel('disk',1));
        gt_id_map(cur_id_map) = cnt;
    end
end
gt_id_map_org = gt_id_map;
cnt_org = cnt;
%% append other cells
gt_id_map = gt_id_map_org;
cnt = cnt_org;
% 1
cnt = cnt + 1;
cell = tifread('1.tif');
cur_id_map = imclose(cell>0,strel('disk', 3));
zzshow(cur_id_map);
gt_id_map(cur_id_map) = cnt;

% 2
cnt = cnt + 1;
cell = tifread('2.tif');
cur_id_map = imclose(cell>0,strel('disk', 3));
zzshow(cur_id_map);
gt_id_map(cur_id_map) = cnt;

% 3
cnt = cnt + 1;
cell = tifread('3.tif');
cur_id_map = imclose(cell>0,strel('disk', 3));
zzshow(cur_id_map);
gt_id_map(cur_id_map) = cnt;

% vid = scale_image(embryo_vid{1},0,1);
% cur_id_map = false(size(vid));
% bw = vid>(24/255);
% cur_id_map(183:197, 30:40, 1:3) = bw(183:197, 30:40, 1:3);
% tifwrite(double(cur_id_map), 'test');

% 4
cnt = cnt + 1;
cell = tifread('4.tif');
cur_id_map = imclose(cell>0,strel('disk', 3));
zzshow(cur_id_map);
gt_id_map(cur_id_map) = cnt;

% 5
cnt = cnt + 1;
vid = scale_image(embryo_vid{1},0,1);
cur_id_map = false(size(vid));
bw = vid>(12/255);
cur_id_map(125:136, 198:227, 1:2) = bw(125:136, 198:227, 1:2);
cur_id_map = imclose(cur_id_map,strel('disk', 3));
gt_id_map(cur_id_map) = cnt;
% 6
cnt = cnt + 1;
vid = scale_image(embryo_vid{1},0,1);
cur_id_map = false(size(vid));
bw = vid>(13/255);
cur_id_map(175:188, 133:155, 28:30) = bw(175:188, 133:155, 28:30);
cur_id_map = imclose(cur_id_map,strel('disk', 3));
gt_id_map(cur_id_map) = cnt;

% 7 append to existing one
vid = scale_image(embryo_vid{1},0,1);
cur_id_map = false(size(vid));
bw = vid>(24/255);
cur_id_map(145:155, 132:148, 1:2) = bw(145:155, 132:148, 1:2);
cur_id_map = imclose(cur_id_map,strel('disk', 3));
cur_id_map = imdilate(cur_id_map,strel('disk', 1));
gt_id_map(cur_id_map) = gt_id_map(148,137,4);


zzshow(label2rgb3d(gt_id_map));

zz = display_seg_res4PE(embryo_vid{1}, gt_id_map);
tifwrite(zz, 'gt_label');
%% 
total_slices = 0;
for i=1:size(gt_id_map,3)
    cur_slice = gt_id_map(:,:,i);
    total_slices = total_slices + length(unique(cur_slice(cur_slice>0)));
end
%% save data
fine_tune = (gt_id_map==54);
fine_tune(:,:,1:26) = 0;
gt_id_map(gt_id_map==54) = 0;
gt_id_map(fine_tune) = 54;

save('frame_one_with_gt.mat', 'gt_id_map', 'vid');

%% test the result of our framework
seg_std = zeros(numel(refine_resAll), 2);
iou = cell(numel(refine_resAll),1);
for i=1:numel(refine_resAll)
    [iou{i},c] = accuracy_measure_seg(gt_id_map, refine_resAll{i}{1});
    seg_std(i,:) = [mean(iou{i}), std(iou{i})];
%[iou,c] = accuracy_measure_seg(gt_id_map, refine_resAll{end}{1});
end
figure;plot(seg_std(:,1), 'LineWidth', 2);
ylabel('Jaccard similarity index (SEG)');
xlabel('# of interations');
seg = refine_resAll{1}{1};
for i=5:10
    display_bnd_map(vid(:,:,i)*1.5, seg(:,:,i))
end

max_proj_vid = max(vid,[],3);
max_proj_seg = max(seg,[],3);
s = regionprops(max_proj_seg, 'Area');
max_proj_seg_l = max_proj_seg;
max_proj_seg_l(ismember(max_proj_seg, find([s.Area]<50)))=0;
display_bnd_map(max_proj_vid*1.5, max_proj_seg_l);

gt_mask = label2rgb3d(gt_id_map);
bg = sum(gt_mask,3)==3;
bg = repmat(bg, 1,1,3);
gt_mask(bg) = 0;
tifwrite(gt_mask,'gt_mask');

zz = display_seg_res4PE(vid, gt_id_map);
tifwrite(zz, 'gt_label');


%% iterations
iou = cell(numel(refine_resAll),1);
acc_its = [];
for i=1:numel(refine_resAll)
    tmpM = movieInfoAll{i};
    if ~isfield(tmpM, 'validGapMaps')
        tmpM.validGapMaps = refine_resAll{i};
    end
    acc_ours = validate_link_acc(tmpM, [], gt_tracks);
    iou{i} = acc_ours;
    acc_its(i) = (sum(acc_ours(:,1))./ sum(acc_ours(:,2)));
end
acc_cell_its = cellfun(@(x) mean(x(:,3)), iou);
%save('C:\Users\Congchao\Google Drive\Reports\Phd_exams\Preliminary exam\test_data\tracking\movieInfo4PE_experiment.mat', 'movieInfoAll','threshold_resAll','-v7.3');
figure;plot(acc_its)
acc_cell_its = cellfun(@(x) mean(x(:,3)), iou);
hold on; plot(acc_cell_its);
ylabel('Linakge accuracy');
xlabel('# of interations');