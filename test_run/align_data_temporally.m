%% align the data temporally
if isunix
    addpath('/home/ccw/Dropbox/cc_ImHandle/');
else
    addpath('C:\Users\Congchao\Dropbox\cc_ImHandle\');
end
addpath('src_code_cellSegment');

%% align cropped data
data_folder = '/home/ccw/storage/Yinan_data/TM0-49/crop/';
save_folder = '/home/ccw/storage/Yinan_data/TM0-49/crop_temporal_aligned/';
tif_files = dir(fullfile(data_folder, '*.tif'));

embryo_vid_org = cell(numel(tif_files), 1); % original data
for i=1:numel(tif_files)
    embryo_vid_org{i} = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
end
tranlations = InitialDriftEstimate(embryo_vid_org);

tx = 0;
ty = 0;
tz = 0;
for i=2:numel(embryo_vid_org)
    tx = tx - tranlations(i, 1);
    ty = ty - tranlations(i, 2);
    tz = tz - tranlations(i, 3);
    tform = affine3d([1 0 0 0; 0 1 0 0; 0 0 1 0; tx, ty, tz, 1]);
    out = imwarp_same(embryo_vid_org{i}, tform);
    tifwrite(uint16(out), ...
        fullfile(save_folder, tif_files(i).name));
end

tifwrite(uint16(embryo_vid_org{1}), ...
        fullfile(save_folder, tif_files(1).name));

moving_reg = imregister(embryo_vid_org{9}, embryo_vid_org{8},...
'translation',optimizer,metric);

%% align the original data (one direction) : failed due to data size
data_folder = '/home/ccw/storage/Yinan_data/TM0-49/full_data/';
save_folder = '/home/ccw/storage/Yinan_data/TM0-49/full_data_temporal_aligned/';

[optimizer, metric] = imregconfig('Monomodal');
tif_files = dir(fullfile(data_folder, '*.tif'));
ref_im = tifread(fullfile(tif_files(1).folder, tif_files(1).name));
ds_sz = round([size(ref_im, 1)/2, size(ref_im, 2)/2, size(ref_im, 3)]);
ref_im = imresize3(ref_im,ds_sz,'method','nearest');
ndim = ndims(ref_im);
tranlations = zeros(numel(tif_files), ndim);
tx = 0;
ty = 0;
tz = 0;
for i=2:numel(tif_files)
    disp(i);
    mov_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    mov_im = imresize3(mov_im,ds_sz,'method','nearest');
    tic;
    tform = imregtform(mov_im, ref_im,...
        'rigid',optimizer,metric);
    toc;
    tranlations(i, :) = -tform.T(ndim+1, 1:ndim);
    tx = tx - tranlations(i, 1);
    ty = ty - tranlations(i, 2);
    tz = tz - tranlations(i, 3);
    tform = affine3d([1 0 0 0; 0 1 0 0; 0 0 1 0; tx, ty, tz, 1]);
    ref_im = imwarp_same(mov_im, tform);
    tifwrite(uint16(ref_im), ...
        fullfile(save_folder, tif_files(i).name));
    ref_im = mov_im;
end


%% align data using max-projected data (xy, xz, yz)
data_folder = '/home/ccw/storage/Yinan_data/TM0-49/full_data/';
save_folder = '/home/ccw/storage/Yinan_data/TM0-49/full_data_temporal_aligned/';

tif_files = dir(fullfile(data_folder, '*.tif'));
ref_im = tifread(fullfile(tif_files(1).folder, tif_files(1).name));
ndim = ndims(ref_im);
tranlations = zeros(numel(tif_files), ndim);
tx = 0;
ty = 0;
tz = 0;
for i=2:numel(tif_files)
    disp(i);
    mov_im = tifread(fullfile(tif_files(i).folder, tif_files(i).name));
    [drift_xyz, ~] = rigid_align_3d_by_max_proj(mov_im, ref_im);
    tranlations(i,:) = drift_xyz;
    tx = tx - tranlations(i, 1);
    ty = ty - tranlations(i, 2);
    tz = tz - tranlations(i, 3);
    disp([drift_xyz, tx, ty, tz]);
    
    tform = affine3d([1 0 0 0; 0 1 0 0; 0 0 1 0; tx, ty, tz, 1]);
    ref_im = imwarp_same(mov_im, tform);
    tifwrite(uint16(ref_im), ...
        fullfile(save_folder, tif_files(i).name));
    ref_im = mov_im;
end