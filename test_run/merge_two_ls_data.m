% To merge the data from two light-sheet microscopies. The reason behind is
% that these two are imaging the left/right right of the same depth in the
% sample.
mydir  = '/home/ccw/storage/Yinan_data/TM0-49/time_0_stackwise';
pairs = cell(4,1);

pairs{1} = {'/*s00_0.tif', '/*s04_0.tif'};
pairs{2} = {'/*s01_0.tif', '/*s05_0.tif'};
pairs{3} = {'/*s02_0.tif', '/*s06_0.tif'};
pairs{4} = {'/*s03_0.tif', '/*s07_0.tif'};

for i = 2:4
    tif_files = dir(fullfile(mydir, pairs{i}{1}));
    left = tifread(fullfile(tif_files(1).folder, tif_files(1).name));
    tif_files = dir(fullfile(mydir, pairs{i}{2}));
    right = tifread(fullfile(tif_files(1).folder, tif_files(1).name));

    out = uint16(double(left + right) / 2);

    tifwrite(out, fullfile(mydir, [tif_files(1).name(1:end-3), '_',...
        num2str(i), '+', num2str(i+4), '_0']));

end