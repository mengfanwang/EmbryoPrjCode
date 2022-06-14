function out_ims = display_cl_frame(movieInfo, refine_res, savefolder, file_names)
% display the segmentation results frame by frame, the same cell are
% labeled with the same color
%INPUT:
% movieInfo: the structure containing particle2track information
% refine_res: maps consisting the labels of all regions
% savefolder:
% file_names: cells consisting the names to save the files, if it is a char
% variable, it is appended with number 1:framenumber
%OUTPUT: None

% contact: ccwang@vt.edu, 02/27/2020
num_fr = max(movieInfo.frames);
if nargin == 2
    savefolder = [];
    file_names = '';
else
    if isa(file_names, 'char')
        tmp_name = file_names;
        file_names = cell(num_fr, 1);
        for i=1:num_fr, file_names{i} = [tmp_name,'_',num2str(i)];end
    elseif ~iscell(file_names)
        error('file names should be cell or char data type');
    end
end
max_track_num = numel(movieInfo.tracks);
cmap = jet(double(max_track_num) + 1);
cmap = cmap(randperm(max_track_num),:);
cmap(end,:) = [1 1 1];
id_increment = 0;
out_ims = cell(num_fr, 1);
for i=1:num_fr
    fprintf('label %d frame with colors\n', i);
    synId_ref = refine_res{i};
    cell_ids = [1:max(synId_ref(:))] + id_increment;
    trackIds = movieInfo.particle2track(cell_ids);
    trackIds(isnan(trackIds)) = size(cmap,1);
    syn_ref_cmap = cmap(trackIds,:);
    out_ims{i} = label2rgb3d(synId_ref, syn_ref_cmap ,[0 0 0]);
    if ~isempty(savefolder)
        tifwrite(out_ims{i}, ...
            fullfile(savefolder, file_names{i}));
    end
    id_increment = id_increment + length(cell_ids);
end

end