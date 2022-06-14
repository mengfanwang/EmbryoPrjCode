function combined = writeRefineCell(vidComp, newLabel, regComp, regId, save_folder)
% write the region to a given folder
vidComp = scale_image(vidComp, 0, 1);
[h,w,z] = size(vidComp);
maxId = max(newLabel(:));
if maxId==0
    newLabel = regComp;
    maxId = 1;
    for ll = 1:maxId
        combined = zeros(h,w,3,z);
        for j=1:z
            combined(:,:,1,j) = 0.5*(newLabel(:,:,j)==ll); % detected region
            combined(:,:,2,j) = vidComp(:,:,j);
        end
        if ~isempty(save_folder)
            tifwrite(combined, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(regId),'_', num2str(ll),'_comp']));
        end
    end
else
    for ll = 1:maxId
        combined = zeros(h,w,3,z);
        for j=1:z
            %combined(:,:,1,j) = 0.5*regComp(:,:,j)>0;
            combined(:,:,3,j) = 0.5*(newLabel(:,:,j)==ll); % detected region
            combined(:,:,2,j) = vidComp(:,:,j);
        end
        if ~isempty(save_folder)
            tifwrite(combined, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(regId),'_', num2str(ll),'_comp']));
        end
    end
end
end