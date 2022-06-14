i = 4;
load('matlab.mat', 'newLabel');
grad_pc = [0.01, 1];
eig2d = principalCv2d(vidComp.^2, newLabel, 3); % use 2d principal curvature
eig3d = principalCv3d(vidComp.^2, newLabel, 3);

posEigMap = eig2d;
%posEigMap = eigComp;
if grad_pc(1) == 0
    scoreMap = eig3d;
else
    posEigMap(posEigMap<0) = 0;
    scoreMap = (scale_image(dMapComp, 0,grad_pc(1)) + ...
        scale_image(posEigMap, 0,grad_pc(2)))./norm(grad_pc);
end
[newLabel, regTest, seedRegion] = region_refine(newLabel, regTest, scoreMap,...
    eigComp~=0, 10, [1 2]);
inLabel = seedRegion; % seed label
valid_newLabel = newLabel(:)>0;
newIdMap(linerInd(valid_newLabel)) = newLabel(valid_newLabel) + regCnt;
regCnt = regCnt + n;

vidComp = crop3D(vid, yxz, shift);
vidComp = scale_image(vidComp, 0, 1);
scoreMap = scale_image(scoreMap, 0, 1);
[h,w,z] = size(regTest);
combined1 = zeros(h,w,3,z);
combined = zeros(h,w,3,z);
combined2 = zeros(h,w,3,z);
gradCombined = zeros(h,w,3,z);
for j=1:z
    %combined(:,:,1,j) = 0.8*(regTest(:,:,j)>0);
    combined(:,:,3,j) = 0.5*(newLabel(:,:,j)>0); % detected region
    combined(:,:,2,j) = vidComp(:,:,j);
    
    combined1(:,:,1,j) = 0.5*regTest(:,:,j)>0; % detected gap
    combined1(:,:,2,j) = vidComp(:,:,j);
    %             combined(:,:,2,j) = L(:,:,j)>0;
    combined2(:,:,1,j) = 0.1*eigCompMap(:,:,j)>0; % principal cur
    combined2(:,:,2,j) = vidComp(:,:,j);
    
    gradCombined(:,:,2,j) = vidComp(:,:,j);
    gradCombined(:,:,1,j) = scoreMap(:,:,j); % gradient
end
%zzshow(combined)
tifwrite(combined, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_comp']));
%         tifwrite(combined1, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_gap']));
%         tifwrite(vidComp, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_org']));

tifwrite(combined2, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_pc']));
tmprgb = label2rgb3d(newLabel,'jet', [0 0 0], 'shuffle');
tifwrite(tmprgb, ...
    fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_newLabel']));
tmprgb = label2rgb3d(inLabel,'jet', [0 0 0], 'shuffle');
tifwrite(tmprgb, ...
    fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_inLabel']));
tifwrite(scale_image(vidComp,0,1), ...
    fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_input']));
tifwrite(gradCombined, fullfile(save_folder, ['OneComp3d_gcut_PrCv_',num2str(i),'_grad']));
