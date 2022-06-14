function detections = cumulateXMLs(file_path)
% input: file path of all the mat files
% 
% all xmls has been re-save as .mat files
% structure name:GMEMfinalResult
% properties the same as those in xmls files

% centers of all the detections 
matfiles = dir(fullfile(file_path, '*.mat'));
detections = cell(numel(matfiles), 1);
for i=1:numel(matfiles)
    load(fullfile(matfiles(i).folder, matfiles(i).name));
    centers = extract(GMEMfinalResult, 'm');
    detections{i} = centers;
end

% matfiles = dir(fullfile(file_path, '*.mat'));
% detections = cell(numel(matfiles), 1);
% for i=1:numel(matfiles)
%     disp(i);
%     load(fullfile(matfiles(i).folder, matfiles(i).name));
%     if isempty(GMEMfinalResult)
%         continue;
%     end
%     %centers = extract(GMEMfinalResult, 'm');
%     centers = [GMEMfinalResult(:).m];
%     parents = [GMEMfinalResult(:).parent];
%     centers = reshape(centers, 3, [])';
%     detections{i} = [centers, parents'];
% end


end