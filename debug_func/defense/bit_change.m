from_path = '/home/ccw/storage/downsample_embryo_data_700x700x500_all';
to_path = '/home/ccw/storage/Mouse Embryogenesis/';
files = dir(fullfile(from_path, '*.tif'));

for i=1:2:numel(files)
    disp(i);
    data = tifread(fullfile(files(i).folder, files(i).name));
    data = uint8(255 * data/5000);
    tifwrite(data, fullfile(to_path, files(i).name));
end