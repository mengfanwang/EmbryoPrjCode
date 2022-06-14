function [drift_xyz, tranlations] = rigid_align_3d_by_max_proj(mov3d, fixed3d)

ndim = 3;
[optimizer, metric] = imregconfig('Monomodal');
tranlations = zeros(ndim, ndim-1);
for i=1:ndim % xz, yz, yx
    mov_proj = squeeze(max(mov3d, [], i));
    fixed_proj = squeeze(max(fixed3d, [], i));
    tform = imregtform(mov_proj, fixed_proj,...
        'translation',optimizer,metric);
    tranlations(i,:) = -tform.T(ndim, 1:ndim-1);
end
drift_xyz = zeros(1,ndim); % y, x, z
drift_xyz(1) = tranlations(3, 1);
drift_xyz(2) = tranlations(3, 2);
drift_xyz(3) = tranlations(1, 2);
end