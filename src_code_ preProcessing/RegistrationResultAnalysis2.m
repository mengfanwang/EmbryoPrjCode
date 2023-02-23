clc;clear;

folder_name = '/run/user/1001/gvfs/smb-share:server=rs0001.local,share=shared_yinan/tracking_50_149/';
link_table = readmatrix(fullfile(folder_name, 'mastodon_1127/MastodonTable-Link.csv'));
link_table = link_table(3:end, 4:6);
spot_table = readmatrix(fullfile(folder_name, 'mastodon_1127/MastodonTable-Spot.csv'));
spot_table = spot_table(3:end, [1 5 6 7]);

start_ind = [2266 1839 482];    % 1874 not in 
start_time = [1 1 1];
[tracks1, gt_loc1] = getTracks(spot_table, link_table, start_ind, start_time);

folder_name = '/run/user/1001/gvfs/smb-share:server=rs0001.local,share=shared_yinan/tracking_50_149/';
link_table = readmatrix(fullfile(folder_name, 'mastodon_1205/MastodonTable-Link.csv'));
link_table = link_table(3:end, 4:6);
spot_table = readmatrix(fullfile(folder_name, 'mastodon_1205/MastodonTable-Spot.csv'));
spot_table = spot_table(3:end, [1 5 6 7]);

start_ind = [1220 666 984 662];    % 1874 not in 
start_time = [1 1 1 1];
[tracks2, gt_loc2] = getTracks(spot_table, link_table, start_ind, start_time);
tracks2{10} = tracks2{10}(1:42,:);
gt_loc2{10} = gt_loc2{10}(1:41,:);
gt_loc2{10} = [gt_loc2{10};([1434 818 452]+1)./[2 2 5.86]];
tracks2 = tracks2(1:13);
gt_loc2 = gt_loc2(1:13);

% outlier: 415 913 953 1020
% annotation1:482  {11}(19) 415     true outlier, but reasonable prediction
%                                   [1114 567 524 48+49]
% annotation2:984  {21}(2)  913     true outlier,  
%                                   [1284 570 459 58+49]
% annotation2:984  {21}(42) 953     error              find by layer 4
% annotation2:662  {23}(23) 1020    error            

%%
tracks = [tracks1; tracks2];
gt_loc = [gt_loc1; gt_loc2];

preds_new = cell(size(tracks));
for ii = 1:length(preds_new)
    preds_new{ii} = zeros(size(tracks{ii},1), 3);
    preds_new{ii}(end,:) = gt_loc{ii}(end,:);
    ts = tracks{ii}(1,2);
    te = tracks{ii}(end,2);
    for tt = te-1:-1:ts
      load(['/work/Mengfan/Embryo/Registration/nonrigid_result_view11_l5_s1_linear/', num2str(tt+49), '.mat']);
        [x_bias, y_bias, z_bias] = moving_predict(phi_current_vec, preds_new{ii}(tt-ts+2,2), preds_new{ii}(tt-ts+2,1), preds_new{ii}(tt-ts+2,3),5);  
        preds_new{ii}(tt-ts+1,2) = preds_new{ii}(tt-ts+2,2) + x_bias;
        preds_new{ii}(tt-ts+1,1) = preds_new{ii}(tt-ts+2,1) + y_bias;
        preds_new{ii}(tt-ts+1,3) = preds_new{ii}(tt-ts+2,3) + z_bias;

    end
end

%plot
for ii = 1:length(tracks)
    if size(tracks{ii},1) > 50
    plot3(gt_loc{ii}(:,1), gt_loc{ii}(:,2), gt_loc{ii}(:,3), 'Color', [0 0.4470 0.7410], 'LineWidth',1); hold on;
    plot3(preds_new{ii}(:,1), preds_new{ii}(:,2), preds_new{ii}(:,3), 'Color', [0.8500 0.3250 0.0980], 'LineWidth',1); hold on;
    end
end
axis([0 960 0 960 0 160]);

%%
% tracks = tracks1;
% gt_loc = gt_loc1;
tracks = [tracks1; tracks2];
gt_loc = [gt_loc1; gt_loc2];


load('/work/public/sameViewFusion/sameViewDetection_050-149_11/tform_050-149_11_translation.mat');
xx = zeros(99,1);for ii = 1:99;xx(ii) = tform{ii}.T(4,2);end
yy = zeros(99,1);for ii = 1:99;yy(ii) = tform{ii}.T(4,1);end
zz = zeros(99,1);for ii = 1:99;zz(ii) = tform{ii}.T(4,3);end

% prediction
preds_new = cell(size(tracks));
preds_new5 = cell(size(tracks));
preds_rigid = cell(size(tracks));
for ii = 1:length(preds_new)
    preds_new{ii} = zeros(size(tracks{ii},1)-1, 3);
    preds_new5{ii} = zeros(size(tracks{ii},1)-1, 3);
    preds_rigid{ii} = zeros(size(tracks{ii},1)-1, 3);
    ts = tracks{ii}(1,2);
    te = tracks{ii}(end,2);
    for tt = ts:te-1
        load(['/work/Mengfan/Embryo/Registration/nonrigid_result_view11_l4_s1_linear/', num2str(tt+49), '.mat']);
        [preds_new{ii}(tt-ts+1,2), preds_new{ii}(tt-ts+1,1), preds_new{ii}(tt-ts+1,3)] = moving_predict(phi_current_vec, gt_loc{ii}(tt-ts+1,2), gt_loc{ii}(tt-ts+1,1), gt_loc{ii}(tt-ts+1,3),4);  
        load(['/work/Mengfan/Embryo/Registration/nonrigid_result_view11_l5_s1_linear/', num2str(tt+49), '.mat']);
        [preds_new5{ii}(tt-ts+1,2), preds_new5{ii}(tt-ts+1,1), preds_new5{ii}(tt-ts+1,3)] = moving_predict(phi_current_vec, gt_loc{ii}(tt-ts+1,2), gt_loc{ii}(tt-ts+1,1), gt_loc{ii}(tt-ts+1,3),5);  
        
        preds_rigid{ii}(tt-ts+1,1) = yy(tt);
        preds_rigid{ii}(tt-ts+1,2) = xx(tt);
        preds_rigid{ii}(tt-ts+1,3) = zz(tt);
    end
end

%%
distance_old_all = [];
distance_new_all = [];
distance_new5_all = [];
distance_rigid_all = [];
distance_temp_all = [];
for ii = 1:length(tracks)
    distance_old = sqrt(sum( (diff(gt_loc{ii})).^2 , 2));
    distance_new = sqrt(sum( (diff(gt_loc{ii}) + preds_new{ii}).^2 , 2));
    distance_new5 = sqrt(sum( (diff(gt_loc{ii}) + preds_new5{ii}).^2 , 2));
    distance_rigid = sqrt(sum( (diff(gt_loc{ii}) - preds_rigid{ii}).^2 , 2));
    
    d_temp = diff(gt_loc{ii}) + preds_new{ii};
%     distance_pred = sum( (diff(gt_loc{ii})).^2, 2);
    distance_temp = d_temp(:,2);
    

    distance_old_all = [distance_old_all; distance_old];
    distance_new_all = [distance_new_all; distance_new];
    distance_new5_all = [distance_new5_all; distance_new5];
    distance_rigid_all = [distance_rigid_all; distance_rigid];
    distance_temp_all = [distance_temp_all; distance_temp];

    ts = tracks{ii}(1,2);
    te = tracks{ii}(end,2);
%     plot(ts:te-1, distance_old - distance_new); hold on
end


%%
[xx, yy, zz] = moving_predict(phi_current_vec, 285.54, 643.52, 78.64, 5)


function [tracks, locations] = getTracks(spot_table, link_table, start_ind, start_time)

    traj_ind = 0;
    tracks = cell(length(start_ind)*4,1);
    while ~isempty(start_ind)
        prev_ind = start_ind(1);
        start_ind(1) = [];
        prev_time = start_time(1);
        start_time(1) = [];

        traj_ind = traj_ind + 1;
        tracks{traj_ind} = [prev_ind prev_time];
        while 1
            row_ind = find(ismember(link_table(:,1), prev_ind));
            curr_ind = link_table(row_ind, 2);
            if length(curr_ind) >= 2
                warning('Division find!');
                start_ind = [curr_ind' start_ind];
                start_time = [repmat(prev_time+1, 1, length(curr_ind)) start_time];
                break;
            elseif isempty(curr_ind)
                break;
            else
                prev_ind = curr_ind;
                prev_time = prev_time+1;
                tracks{traj_ind} = [tracks{traj_ind}; [prev_ind prev_time]];
            end
        end
    end
    tracks = tracks(1:traj_ind, 1);
    
    % get location
    locations = cell(traj_ind, 1);
    for ii = 1:traj_ind
        locations{ii} = spot_table(tracks{ii}(:,1)+1,2:4);
        locations{ii} = (locations{ii}+1)./[2 2 5.86];
    end

end

function [x_bias, y_bias, z_bias] = moving_predict(phi_current_vec, x_start, y_start, z_start, layer)
    [y_batch, x_batch, z_batch] = meshgrid(0:2^layer+1);
    grid_size = [2^layer 2^layer 2^layer];
    batch_size = round([960 960 160]./grid_size);
    resize_method = 'linear';

    y_batch = y_batch*batch_size(2) + 0.5 - batch_size(2)/2;
    x_batch = x_batch*batch_size(1) + 0.5 - batch_size(1)/2;
    z_batch = z_batch*batch_size(3) + 0.5 - batch_size(3)/2;
    x_pad = padarray(reshape(phi_current_vec(1:3:end-2), grid_size), [1 1 1], 'replicate'); 
    x_bias = interp3(y_batch, x_batch, z_batch, x_pad, y_start, x_start, z_start, resize_method);
    y_pad = padarray(reshape(phi_current_vec(2:3:end-1), grid_size), [1 1 1], 'replicate'); 
    y_bias = interp3(y_batch, x_batch, z_batch, y_pad, y_start, x_start, z_start, resize_method);
    z_pad = padarray(reshape(phi_current_vec(3:3:end), grid_size), [1 1 1], 'replicate'); 
    z_bias = interp3(y_batch, x_batch, z_batch, z_pad, y_start, x_start, z_start, resize_method);
end