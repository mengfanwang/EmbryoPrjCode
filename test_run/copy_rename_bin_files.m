bin_path = ['/home/ccw/Desktop', ...
    '/embryo_res_folder/downsample_crop_embryo_data_470x350x250x50'];
inputFiles = dir(fullfile(bin_path, '*.bin'));
outputFiles = dir(fullfile(bin_path, '*.tif'));
fileNames = { inputFiles.name };
outfileNames = {outputFiles.name};
for k = 1 : length(inputFiles )
  thisFileName = fileNames{k};
  % Prepare the input filename.
  inputFullFileName = fullfile(bin_path, thisFileName);
  % Prepare the output filename.
  for j = 1 : length(outfileNames)
      outputBaseFileName = sprintf('%s%s.bin', outfileNames{j}(1:end-4),...
          thisFileName(1:end-4));
      outputFullFileName = fullfile(bin_path, outputBaseFileName);
      % Do the copying and renaming all at once.
      copyfile(inputFullFileName, outputFullFileName);
  end
end